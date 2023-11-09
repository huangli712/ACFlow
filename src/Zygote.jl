module Zygote

import ZygoteRules: @adjoint, @adjoint!, AContext, adjoint, _pullback, pullback,
  literal_getproperty, literal_getfield
using ChainRulesCore
using ChainRules: rrule, unthunk
using IRTools
using MacroTools, Requires
using IRTools: IR, Variable, Pipe, xcall, var, prewalk, postwalk,
  blocks, predecessors, successors, argument!, arguments, branches,
  insertafter!, finish, expand!, prune!, substitute!, substitute,
  block, block!, branch!, return!, stmt, meta
using IRTools.Inner: argnames!, update!
using IRTools: varargs!, inlineable!, pis!, slots!
import Base: copy!, tail, RefValue
using Base.Broadcast: AbstractArrayStyle, broadcasted
using Distributed: pmap
using ForwardDiff: Dual, value, partials

export gradient

const Numeric{T<:Number} = Union{T, AbstractArray{<:T}}

macro __splatnew__(T, args)
  esc(Expr(:splatnew, T, args))
end

@inline __new__(T, args...) = @__splatnew__(T, args)
@inline __splatnew__(T, args) = @__splatnew__(T, args)

literal_getindex(x, ::Val{i}) where i = getindex(x, i)
literal_indexed_iterate(x, ::Val{i}) where i = Base.indexed_iterate(x, i)
literal_indexed_iterate(x, ::Val{i}, state) where i = Base.indexed_iterate(x, i, state)

@inline tuple_va(N, xs) = xs
@inline tuple_va(N, x, xs...) = (x, tuple_va(N, xs...)...)
@inline tuple_va(::Val{N}, ::Nothing) where N = ntuple(_ -> nothing, Val(N))

iscall(x, m::Module, n::Symbol) = isexpr(x, :call) && x.args[1] == GlobalRef(m, n)

gradindex(x, i) = x[i]
gradindex(::Nothing, i) = nothing
xgetindex(x, i...) = xcall(Base, :getindex, x, i...)
xgradindex(x, i) = xcall(Zygote, :gradindex, x, i)

normalise!(ir) = ir |> IRTools.merge_returns!

function instrument_new!(ir, v, ex)
  isexpr(ex, :new) ? (ir[v] = xcall(Zygote, :__new__, ex.args...)) :
  isexpr(ex, :splatnew) ? (ir[v] = xcall(Zygote, :__splatnew__, ex.args...)) :
  ex
end

# Hack to work around fragile constant prop through overloaded functions
unwrapquote(x) = x
unwrapquote(x::QuoteNode) = x.value

is_getproperty(ex) = iscall(ex, Base, :getproperty)

# The initial premise of literal_getproperty was in some ways inherently flawed, because for
# getproperty it was intended that _pullback falls back to literal_getproperty, but we actually
# want the opposite to happen, since Zygote should fall back to recursing into the getproperty
# implementation by default. Users still want to define custom adjoints using only
# literal_getproperty, though. We can't really have mutually recursive definitions here, so we
# now always instrument getproperty as literal_getproperty, no matter whether the second
# argument is a literal or not.
function instrument_getproperty!(ir, v, ex)
  if is_getproperty(ex)
    obj, prop = ex.args[2], ex.args[3]
    if obj isa Module && prop isa QuoteNode && isconst(obj, unwrapquote(prop))
      # Metaprogramming can generate getproperty(::Module, ...) calls.
      # Like other types, these are type unstable without constprop.
      # However, literal_getproperty's heuristic is also not general enough for modules.
      # Thankfully, we can skip instrumenting these if they're const properties.
      ex
    elseif prop isa Union{QuoteNode,Integer}
      ir[v] = xcall(Zygote, :literal_getproperty, obj, Val(unwrapquote(prop)))
    else
      f = insert!(ir, v, :(Val($(prop))))
      ir[v] = xcall(Zygote, :literal_getproperty, obj, f)
    end
  else
    ex
  end
end

# Here, only instrumenting getfield with literals is fine, since users should never have to
# define custom adjoints for literal_getfield
function instrument_getfield!(ir, v, ex)
  if is_literal_getfield(ex)
    ir[v] = xcall(Zygote, :literal_getfield, ex.args[2], Val(unwrapquote(ex.args[3])))
  else
    ex
  end
end

is_literal_getfield(ex) =
  (iscall(ex, Core, :getfield) || iscall(ex, Base, :getfield)) &&
  ex.args[3] isa Union{QuoteNode,Integer}

is_literal_iterate(ex) =
  iscall(ex, Base, :indexed_iterate) && length(ex.args) >= 3 && ex.args[3] isa Union{Integer,QuoteNode}

is_literal_getindex(ex) =
  iscall(ex, Base, :getindex) && length(ex.args) == 3 && ex.args[3] isa Union{Integer,QuoteNode}

# TODO: is this always correct for user defined getindex methods?
function instrument_getindex!(ir, v, ex)
  if is_literal_getindex(ex)
    ir[v] = xcall(Zygote, :literal_getindex, ex.args[2], Val(unwrapquote(ex.args[3])))
  else
    ex
  end
end

function instrument_iterate!(ir, v, ex)
  if is_literal_iterate(ex)
    ir[v] = xcall(Zygote, :literal_indexed_iterate, ex.args[2],
                  Val(unwrapquote(ex.args[3])), ex.args[4:end]...)
  else
    ex
  end
end

function instrument_literals!(ir, v, ex)
  ex = instrument_getproperty!(ir, v, ex)
  ex = instrument_getfield!(ir, v, ex)
  ex = instrument_getindex!(ir, v, ex)
  ex = instrument_iterate!(ir, v, ex)
end

function instrument_global!(ir, v, ex)
  if istrackable(ex)
    ir[v] = xcall(Zygote, :unwrap, QuoteNode(ex), ex)
  else
    ir[v] = prewalk(ex) do x
      istrackable(x) || return x
      insert!(ir, v, xcall(Zygote, :unwrap, QuoteNode(x), x))
    end
  end
end

function istrackable(x)
  x isa GlobalRef && x.mod ∉ (Base, Core) || return false
  isconst(x.mod, x.name) || return true
  x = getfield(x.mod, x.name)
  !(x isa Type || sizeof(x) == 0)
end

function instrument(ir::IR)
  pr = Pipe(ir)
  for (v, st) in pr
    ex = st.expr
    if isexpr(ex, :foreigncall, :isdefined)
      continue
    elseif isexpr(ex, :enter, :leave)
      error("""try/catch is not supported.
            Refer to the Zygote documentation for fixes.
            https://fluxml.ai/Zygote.jl/latest/limitations
            """)
    elseif isexpr(ex, :(=))
      @assert ex.args[1] isa GlobalRef
      pr[v] = xcall(Zygote, :global_set, QuoteNode(ex.args[1]), ex.args[2])
    else
      ex = instrument_new!(pr, v, ex)
      ex = instrument_literals!(pr, v, ex)
      ex = instrument_global!(pr, v, ex)
    end
  end
  ir = finish(pr)
  # GlobalRefs can turn up in branch arguments
  for b in blocks(ir), br in branches(b), i in 1:length(arguments(br))
    (ref = arguments(br)[i]) isa GlobalRef || continue
    arguments(br)[i] = push!(b, xcall(Zygote, :unwrap, QuoteNode(ref), ref))
  end
  return ir
end

const BranchNumber = UInt8

function record_branches!(ir::IR)
  brs = Dict{Int,Variable}()
  for bb in blocks(ir)
    preds = predecessors(bb)
    length(preds) > 1 || continue
    brs[bb.id] = argument!(bb, BranchNumber(0), BranchNumber)
    i = length(arguments(bb))
    n = 0
    for aa in blocks(ir), br in branches(aa)
      br.block == bb.id && (arguments(br)[i] = BranchNumber(n += 1))
    end
  end
  return ir, brs
end

ignored_f(f) = f in (GlobalRef(Base, :not_int),
                     GlobalRef(Core.Intrinsics, :not_int),
                     GlobalRef(Core, :(===)),
                     GlobalRef(Core, :apply_type),
                     GlobalRef(Core, :typeof),
                     GlobalRef(Core, :throw),
                     GlobalRef(Base, :kwerr),
                     GlobalRef(Core, :kwfunc),
                     GlobalRef(Core, :isdefined))
ignored_f(ir, f) = ignored_f(f)
ignored_f(ir, f::Variable) = ignored_f(get(ir, f, nothing))

function ignored(ir, ex)
  isexpr(ex, :call) || return false
  f = ex.args[1]
  ignored_f(ir, f) && return true
  if f isa Variable && haskey(ir, f)
    f = ir[f].expr
  end
  if f == GlobalRef(Base, :getproperty) && length(ex.args) >= 3
    obj, prop = ex.args[2], ex.args[3]
    # Metaprogramming can generate getproperty(::Module, ...) calls.
    # These are type unstable without constprop, which transforming to _pullback breaks.
    # However, we can skip differentiating these if they're const properties.
    obj isa Module && prop isa QuoteNode && isconst(obj, unwrapquote(prop)) && return true
  end
  return false
end

function primal(ir::IR)
  pr = Pipe(ir)
  pbs = Dict{Variable,Variable}()
  argument!(pr, at = 1)
  cx = argument!(pr, Context, at = 2)
  for (v, st) in pr
    ex = st.expr
    if isexpr(ex, :call) && !ignored(ir, ex)
      yJ = insert!(pr, v, stmt(xcall(Zygote, :_pullback, cx, ex.args...),
                               line = ir[v].line))
      pr[v] = xgetindex(yJ, 1)
      J = insertafter!(pr, v, stmt(xgetindex(yJ, 2),
                                   line = ir[v].line))
      pbs[v] = substitute(pr, J)
    end
  end
  pr = finish(pr)
  pr, brs = record_branches!(pr)
  return pr, brs, pbs
end

struct Primal
  ir::IR
  pr::IR
  varargs::Union{Int,Nothing}
  branches::Dict{Int,Variable}
  pullbacks::Dict{Variable,Variable}
end

function Primal(ir::IR; varargs = nothing)
  ir = instrument(normalise!(ir))
  pr, brs, pbs = primal(ir)
  Primal(expand!(ir), pr, varargs, brs, pbs)
end

# Backwards Pass

struct Alpha
  id::Int
end

Base.show(io::IO, x::Alpha) = print(io, "@", x.id)

alpha(x) = x
alpha(x::Variable) = Alpha(x.id)
Variable(a::Alpha) = Variable(a.id)

sig(b::IRTools.Block) = unique([arg for br in branches(b) for arg in br.args if arg isa Variable])
sig(pr::Primal) = Dict(b.id => sig(b) for b in blocks(pr.ir))

function adjointcfg(pr::Primal)
  ir = empty(pr.ir)
  return!(ir, nothing)
  for b in blocks(pr.ir)[2:end]
    block!(ir)
    preds = predecessors(b)
    rb = block(ir, b.id)
    for i = 1:length(preds)
      cond = i == length(preds) ? nothing :
        push!(rb, xcall(Base, :(!==), alpha(pr.branches[b.id]), BranchNumber(i)))
      branch!(rb, preds[i].id, unless = cond)
    end
    if isempty(preds) || (!isempty(branches(b)) && branches(b)[end] == IRTools.unreachable)
      # If `b` is unreachable, then no context produced by the primal should end up branching to `rb`
      push!(rb, xcall(Core, :throw, "unreachable")) # `throw` is necessary for inference not to hit the `unreachable`
      branch!(rb, 0)
    end
  end
  sigs = sig(pr)
  for b in blocks(ir)[1:end-1], i = 1:length(sigs[b.id])
    argument!(b)
  end
  argument!(blocks(ir)[end])
  return ir, sigs
end

branchfor(ir, (from,to)) =
  get(filter(br -> br.block == to, branches(block(ir, from))), 1, nothing)

xaccum(ir) = nothing
xaccum(ir, x) = x
xaccum(ir, xs...) = push!(ir, xcall(Zygote, :accum, xs...))

function passthrough_expr(ex::Expr)
    # Metadata we want to preserve
    isexpr(ex, GlobalRef, :call, :isdefined, :inbounds, :meta, :loopinfo) && return true
    # ccalls and more that are safe to preserve/required for proper operation:
    # - jl_set_task_threadpoolid: added in 1.9 for @spawn
    isexpr(ex, :foreigncall) && unwrapquote(ex.args[1]) in (:jl_set_task_threadpoolid,) && return true
    return false
end

function adjoint(pr::Primal)
  #@show "adjoint"
  ir, sigs = adjointcfg(pr)
  for b in reverse(blocks(pr.ir))
    rb = block(ir, b.id)
    grads = Dict()
    grad(x, x̄) = push!(get!(grads, x, []), x̄)
    grad(x) = xaccum(rb, get(grads, x, [])...)
    # Backprop through (successor) branch arguments
    for i = 1:length(sigs[b.id])
      grad(sigs[b.id][i], arguments(rb)[i])
    end
    # Backprop through statements
    for v in reverse(keys(b))
      ex = b[v].expr
      if haskey(pr.pullbacks, v)
        g = push!(rb, stmt(Expr(:call, alpha(pr.pullbacks[v]), grad(v)),
                           line = b[v].line))
        for (i, x) in enumerate(ex.args)
          x isa Variable || continue
          grad(x, push!(rb, stmt(xgradindex(g, i),
                                 line = b[v].line)))
        end
      elseif ex isa Core.PiNode
        grads[ex.val] = grads[v]
      elseif isexpr(ex) && !passthrough_expr(ex)
        push!(rb, stmt(xcall(Base, :error, """
                             Can't differentiate $(ex.head) expression $ex.
                             You might want to check the Zygote limitations documentation.
                             https://fluxml.ai/Zygote.jl/latest/limitations
                             """),
                       line = b[v].line))
      else # A literal value
        continue
      end
    end
    if b.id > 1 # Backprop through (predecessor) branch arguments
      gs = grad.(arguments(b))
      for br in branches(rb)
        br.block == 0 && continue
        br′ = branchfor(pr.ir, br.block=>b.id)
        br′ === nothing && continue
        ins = br′.args
        for i = 1:length(br.args)
          ā = [gs[j] for j = 1:length(ins) if ins[j] == sigs[br.block][i]]
          br.args[i] = xaccum(rb, ā...)
        end
      end
    else # Backprop function arguments
      gs = [grad(arg) for arg = arguments(pr.ir)]
      Δ = push!(rb, pr.varargs === nothing ?
                      xcall(Zygote, :tuple, gs...) :
                      xcall(Zygote, :tuple_va, Val(pr.varargs), gs...))
      branches(rb)[1].args[1] = Δ
    end
  end
  return ir
end

struct Adjoint
  primal::IR
  adjoint::IR
end

function Adjoint(ir::IR; varargs = nothing, normalise = true)
  pr = Primal(ir, varargs = varargs)
  adj = adjoint(pr) |> prune!
  if normalise
    permute!(adj, length(adj.blocks):-1:1)
    adj = IRTools.domorder!(adj) |> IRTools.renumber
  end
  Adjoint(pr.pr, adj)
end

# Stacks

mutable struct Stack{T}
  idx::Int
  data::Vector{T}
end

Stack(data::Vector{T}) where T =
  Stack{T}(length(data), data)

function Base.pop!(stk::Stack)
  i = stk.idx
  stk.idx = i == 1 ? length(stk.data) : i-1
  @inbounds return stk.data[i]
end

function _push!(a::Vector{T}, x::T) where T
  Base._growend!(a, 1)
  @inbounds a[end] = x
  return
end

# Emit

xstack(T) = Expr(:call, Vector{T})

function alphauses(b)
  us = Set{Alpha}()
  postwalk(x -> x isa Alpha && push!(us, x), b)
  return us
end

xtuple(xs...) = xcall(:tuple, xs...)

concrete(T::DataType) = T
concrete(::Type{Type{T}}) where T = typeof(T)

runonce(b) = b.id in (1, length(b.ir.blocks))

function forward_stacks!(adj, F)
  #@show "forward"
  stks, recs = [], []
  pr = adj.primal
  for b in blocks(pr), α in alphauses(block(adj.adjoint, b.id))
    if runonce(b)
      push!(recs, Variable(α))
    else
      stk = pushfirst!(pr, xstack(Any))
      push!(recs, stk)
      push!(b, xcall(Zygote, :_push!, stk, Variable(α)))
    end
    push!(stks, (b.id, alpha(α)))
  end
  args = arguments(pr)[3:end]
  rec = push!(pr, xtuple(recs...))
  P = length(pr.blocks) == 1 ? Pullback{F} : Pullback{F,Any}
  # P = Pullback{F,Any} # reduce specialisation
  rec = push!(pr, Expr(:call, P, rec))
  ret = xtuple(pr.blocks[end].branches[end].args[1], rec)
  ret = push!(pr, ret)
  pr.blocks[end].branches[end].args[1] = ret
  return pr, stks
end

function reverse_stacks!(adj, stks)
  #@show "reverse"
  ir = adj.adjoint
  entry = blocks(ir)[end]
  self = argument!(entry, at = 1)
  t = pushfirst!(blocks(ir)[end], xcall(:getfield, self, QuoteNode(:t)))
  repl = Dict()
  runonce(b) = b.id in (1, length(ir.blocks))
  for b in blocks(ir)
    for (i, (b′, α)) in enumerate(stks)
      b.id == b′ || continue
      if runonce(b)
        val = insertafter!(ir, t, xcall(:getindex, t, i))
      else
        stk = push!(entry, xcall(:getindex, t, i))
        stk = push!(entry, xcall(Zygote, :Stack, stk))
        val = pushfirst!(b, xcall(:pop!, stk))
      end
      repl[α] = val
    end
  end
  return IRTools.prewalk!(x -> get(repl, x, x), ir)
end

function stacks!(adj, T)
  forw, stks = forward_stacks!(adj, T)
  back = reverse_stacks!(adj, stks)
  permute!(back, length(back.blocks):-1:1)
  IRTools.domorder!(back)
  return forw, back
end

varargs(m::Method, n) = m.isva ? n - m.nargs + 1 : nothing

function _generate_pullback_via_decomposition(T, world)
  (m = meta(T; world)) === nothing && return
  va = varargs(m.method, length(T.parameters))
  forw, back = stacks!(Adjoint(IR(m), varargs = va, normalise = false), T)
  m, forw, back
end

struct ZygoteRuleConfig{CTX<:AContext} <: RuleConfig{Union{HasReverseMode,NoForwardsMode}}
  context::CTX
end

_is_rrule_redispatcher(m::Method) = m.sig == Tuple{typeof(rrule), RuleConfig, Vararg}

"""
  has_chain_rrule(T)

For a type-tuple `T` e.g. `Tuple{typeof(f), Int, Float64}`, checks if there is a `rrule` defined for it.
Excluding the generic fallback.
The first return value is `true` if the `rrule` exists, `false` otherwise.
If it does not, then the second argument is a list of edges to attach to the CodeInfo for a generated function,
such that if a suitable rule is defined later, the generated function will recompile.
"""
function has_chain_rrule(T, world)
  config_T, arg_Ts = Iterators.peel(T.parameters)
  configured_rrule_m = meta(Tuple{typeof(rrule), config_T, arg_Ts...}; world)
  is_ambig = configured_rrule_m === nothing  # this means there was an ambiguity error, on configured_rrule

  if !is_ambig && _is_rrule_redispatcher(configured_rrule_m.method)
    # The config is not being used:
    # it is being redispatched without config, so we need the method it redispatches to
    rrule_m = meta(Tuple{typeof(rrule), arg_Ts...}; world)
    # Thus any no_rrule that might apply must also not have a config because if there was a
    # no_rrule with a config that applied then there would also be a rrule with config that applied
    no_rrule_m = meta(Tuple{typeof(ChainRulesCore.no_rrule), arg_Ts...}; world)
  else
    # Not being redispatched: it does have a config
    rrule_m = configured_rrule_m
    # Thus any no_rrule that might apply must also have a config because if it applied
    # it will be identical, and if it doesn't we don't care what it is.
    no_rrule_m = meta(Tuple{typeof(ChainRulesCore.no_rrule), config_T, arg_Ts...}; world)
  end

  is_ambig |= rrule_m === nothing  # this means there was an ambiguity error on unconfigured rrule

  # To understand why we only need to check if the sigs match between no_rrule_m and rrule_m
  # in order to decide if to use, one must consider the following facts:
  # - for every method in `no_rrule` there is a identical one in `rrule` that returns nothing
  # - this includes the general fallback `rrule(::Any...)=nothing`.
  # - a configured rrule/no_rrule is always more specific than a otherwise equivalent unconfigured rrule/no_rrule
  #  
  # Consider the following truth table, for what can occur:
  # rrule: fallback, no_rrule: fallback =>  matches => do not use rrule.
  # rrule: specific, no_rrule: fallback => !matches => do use rrule, as haven't opted out.
  # rrule: fallback, no_rrule: specific =>  IMPOSSIBLE, every no_rule is identical to some rrule
  # rrule: specific, no_rrule: specific =>  matches => do not use rrule as opted out
  # rrule: specific, no_rrule: general  => !matches => do use rrule as a more specific rrule takes preciedent over more general opted out
  # rrule: general , no_rrule: specific =>  IMPOSSIBLE, every no_rule us identical to some rrule so can't have a more general rrule being hit, as the specific one would hit first
  #
  # Note that the fallback cases are the same outcome as the general cases as fallback is just most general.
  # It can be seen that checking if it matches is the correct way to decide if we should use the rrule or not.

  if !is_ambig && matching_cr_sig(no_rrule_m, rrule_m)  # Not ambiguous, and opted-out.
    # Return instance for configured_rrule_m as that will be invalidated 
    # directly if configured rule added, or indirectly if unconfigured rule added
    # Do not need an edge for `no_rrule` as no addition of methods to that can cause this
    # decision to need to be revisited (only changes to `rrule`), since we are already not
    # using the rrule, so not using more rules wouldn't change anything.
    return false, configured_rrule_m.instance
  else
    # Either is ambiguous, and we should try to use it, and then error
    # or we are uses a rrule, no need to add any edges for `rrule`, as it will generate 
    # code with natural edges if a new method is defined there.
    # We also do not need an edge to `no_rrule`, as any time a method is added to `no_rrule`
    # a corresponding method is added to `rrule` (to return `nothing`), thus we will already
    # be revisiting this decision when a new opt-out is added.
    return true, nothing
  end
end

matching_cr_sig(t, s) = matching_cr_sig(t.method.sig, s.method.sig)
matching_cr_sig(::DataType, ::UnionAll) = false
matching_cr_sig(::UnionAll, ::DataType) = false
matching_cr_sig(t::Type, s::Type) = type_tuple_tail(t) == type_tuple_tail(s)

type_tuple_tail(d::DataType) = Tuple{d.parameters[2:end]...}
function type_tuple_tail(d::UnionAll)
    body = Base.unwrap_unionall(d)
    body_tt = type_tuple_tail(body)
    return Base.rewrap_unionall(body_tt, d)
end

"""
    is_kwfunc(sigt...)

Determines if `sigt` is the type signature of a kwfunction.
Each element of `sigt` should be a type.
Either the first 3 types are a kwfunc type, a NamedTuple and the matching base function type,
or the first argument is the base function type and it is not a kwfunction.
the remaining types in `sigt` are the types of the argument.

"""
is_kwfunc(k, ::Type{<:NamedTuple}, f, args...) = k===Core.kwftype(f)

"""
    wrap_chainrules_output(x)

Convert `x` from the differentials types ChainRules uses to the format Zygote uses internally.
"""
@inline wrap_chainrules_output(x) = x
@inline wrap_chainrules_output(x::AbstractThunk) = wrap_chainrules_output(unthunk(x))  # For now we are just not going to deal with thunks
@inline wrap_chainrules_output(x::Tuple) = map(wrap_chainrules_output, x)

"""
    wrap_chainrules_input(dx)

Convert `dx` from the format Zygote uses internally to differentials types ChainRules uses.
"""
@inline wrap_chainrules_input(dx) = dx
# For arrays, whitelist the safe ones, but always look inside Any[]:
@inline wrap_chainrules_input(dxs::AbstractArray{<:Number}) = dxs
@inline wrap_chainrules_input(dxs::AbstractArray) = map(wrap_chainrules_input, dxs)

"""
  _project(x, dx)

Uses `ChainRulesCore.ProjectTo` to standardise the gradient `dx` for type & shape.
Also handles some Zygote-specific corrections, such as `x::Array, dx::Tuple`.
Safe to apply to arbitrary input.
"""
@inline function _project(x, dx)
  wrap_chainrules_output(ProjectTo(x)(zygote2differential(dx, x)))
end

"""
  ZBack{F}(back) <: Function

Wrapper for a ChainRules pullback `back`, that causes it to follow Zygote conventions.
(A functor here is used rather than a closure to avoid boxing issues);
"""
struct ZBack{F} <: Function
  back::F
end
@inline (s::ZBack)(dy) = wrap_chainrules_output(s.back(wrap_chainrules_input(dy)))
# `nothing->nothing` can be deleted after https://github.com/FluxML/Zygote.jl/issues/603
# though it might be worth keeping as a performance optimization (benchmarking pending)
@inline (s::ZBack)(::Nothing) = nothing

"""
    chain_rrule(config, f, args...)

Returns a the (primal) value of `f(args...)` and a pullback, by invoking `ChainRulesCore.rrule(f, args...)`.
The pullback is appropriately wrapped up to follow Zygote conventions.
"""
@inline function chain_rrule(config, f, args...)
  y, back = rrule(config, f, args...)
  return y, ZBack(back)
end

"""
    zygote2differential(dx, primal)

Convert input `dx` from the Zygote format to the ChainRules differential types.
This is similar to `wrap_chainrules_input(dx)`, but because it gets `primal::T`,
it can turn `NamedTuple`s into `Tangent{T}(...)` not `Tangent{Any}(...)`.
"""
zygote2differential(x, primal) = z2d(x, primal)
z2d(dx, ::Any) = dx
z2d(dx::AbstractArray, primal::AbstractArray) = map(z2d, dx, primal)

# Internal container used to track accumulated gradients of mutable types (including params).
# Type param I ∈ (true, false) indicates whether implicit params are in use.
# By default, this should be false unless pullback(f, ::Params) is called.
mutable struct Context{I} <: AContext
  cache::Union{IdDict{Any,Any},Nothing}
end

Context() = Context{false}(nothing)
cache(cx::Context) = cx.cache === nothing ? (cx.cache = IdDict()) : cx.cache

struct Pullback{S,T}
  t::T
end

Pullback{S}(x) where S = Pullback{S,typeof(x)}(x)

struct CompileError
  T
  e
end

@inline function pullback(f, args...)
  y, back = _pullback(Context(), f, args...)
  y, Δ -> Base.tail(back(Δ))
end

"""
    gradient(f, args...)

Returns a tuple containing `∂f/∂x` for each argument `x`,
the derivative (for scalar `x`) or the gradient.

`f(args...)` must be a real number, see [`jacobian`](@ref) for array output.

See also [`withgradient`](@ref) to keep the value `f(args...)`,
and [`pullback`](@ref) for value and back-propagator.

```jldoctest; setup=:(using Zygote)
julia> gradient(*, 2.0, 3.0, 5.0)
(15.0, 10.0, 6.0)

julia> gradient(x -> sum(abs2,x), [7.0, 11.0, 13.0])
([14.0, 22.0, 26.0],)

julia> gradient([7, 11], 0, 1) do x, y, d
         p = size(x, d)
         sum(x.^p .+ y)
       end
([14.0, 22.0], 2.0, nothing)
```
"""
function gradient(f, args...)
  y, back = pullback(f, args...)
  grad = back(one(y))
  map(_project, args, grad)
end

# Interfaces

accum(x) = x
accum(x, y) =
  x === nothing ? y :
  y === nothing ? x :
  x + y

accum(x, y, zs...) = accum(accum(x, y), zs...)
accum(x::AbstractArray, ys::AbstractArray...) = accum.(x, ys...)

@generated function accum(x::NamedTuple, y::NamedTuple)
  #@show "haha"
  # assumes that y has no keys apart from those also in x
  fieldnames(y) ⊆ fieldnames(x) || throw(ArgumentError("$y keys must be a subset of $x keys"))

  grad(field) = field in fieldnames(y) ? :(y.$field) : :nothing
  Expr(:tuple, [:($f=accum(x.$f, $(grad(f)))) for f in fieldnames(x)]...)
end

function accum(x::RefValue, y::RefValue)
  #@show "baba"
  @assert x === y
  return x
end

# Core functions
accum_param(::Context{false}, _, Δ) = Δ
@generated function accum_param(cx::Context, x, Δ)
  isbitstype(x) && return :(Δ)
  quote
    if haskey(cache(cx), x)
      cache(cx)[x] = accum(cache(cx)[x],Δ)
      return
    else
      return Δ
    end
  end
end

function accum_global(cx::Context, ref, x̄)
  (x̄ === nothing || isconst(ref.mod, ref.name)) && return
  gs = cache(cx)
  gs[ref] = accum(get(gs, ref, nothing), x̄)
  return
end

unwrap(x) = x
@adjoint unwrap(x) = unwrap(x), x̄ -> (accum_param(__context__, x, x̄),)
unwrap(ref, x) = x
@adjoint unwrap(ref, x) = unwrap(x), function (x̄)
  accum_global(__context__, ref, x̄)
  (accum_param(__context__, x, x̄),)
end

function global_set(ref, val)
  @static if VERSION < v"1.9.0-DEV.265"
    ccall(:jl_set_global, Cvoid, (Any, Any, Any),
          ref.mod, ref.name, val)
  else
    setglobal!(ref.mod, ref.name, val)
  end
end

@adjoint! function global_set(ref, x)
  global_set(ref, x), function (x̄)
    gs = cache(__context__)
    x̄ = accum(get(gs, ref, nothing), x̄)
    gs[ref] = nothing
    return (nothing, x̄)
  end
end

@adjoint tuple(xs...) = xs, identity

unapply(t, xs) = _unapply(t, xs)[1]
_unapply(t, xs) = first(xs), tail(xs)
_unapply(t::Tuple{}, xs) = (), xs

function _unapply(t::Tuple, xs)
  t1, xs1 = _unapply(first(t), xs)
  t2, xs2 = _unapply(tail(t), xs1)
  (t1, t2...), xs2
end

@adjoint! function Core._apply_iterate(::typeof(iterate), f, args...)
  y, back = Core._apply(_pullback, (__context__, f), args...)
  st = map(x -> map(_->nothing, x), args)
  y, function (Δ)
    Δ = back(Δ)
    Δ === nothing ? nothing :
      (nothing, first(Δ), unapply(st, Base.tail(Δ))...)
  end
end

@generated nt_nothing(x) = Expr(:tuple, [:($f=nothing) for f in fieldnames(x)]...)
@generated pair(::Val{k}, v, _=nothing) where k = :($k = v,)
@generated pair(::Val{k}, v, ::NamedTuple{keys}) where {k,keys} = k isa Int ? :($(getfield(keys, k)) = v,) : :($k = v,)

@adjoint function literal_getfield(x, ::Val{f}) where f
  val = getfield(x, f)
  function back(Δ)
    accum_param(__context__, val, Δ) === nothing && return
    if isimmutable(x)
      dx = (; nt_nothing(x)..., pair(Val(f), Δ, x)...)
      (_project(x, dx), nothing)
    else
      dx = grad_mut(__context__, x)
      dx[] = (; dx[]..., pair(Val(f), accum(getfield(dx[], f), Δ))...)
      return (dx,nothing)
    end
  end
  unwrap(val), back
end

_pullback(cx::AContext, ::typeof(getfield), x, field_name::Symbol) =
  _pullback(cx, literal_getfield, x, Val(field_name))

grad_mut(x) = Ref{Any}(nt_nothing(x))
grad_mut(cx::Context, x) = _get!(() -> grad_mut(x), cache(cx), x)

# needed for reverse-over-reverse pending rrule for Base.get!
function _get!(default::Base.Callable, ch, x)
  if haskey(ch, x)
    ch[x]
  else
    ch[x] = default()
  end
end

# Dictionaries

grad_mut(d::AbstractDict) = Dict()

@adjoint function getindex(d::AbstractDict, k)
  val = d[k]
  function dict_getindex_pullback(Δ)
    accum_param(__context__, val, Δ) === nothing && return
    grad = grad_mut(__context__, d)
    grad[k] = accum(get(grad, k, nothing), Δ)
    return (grad, nothing)
  end
  val, dict_getindex_pullback
end

for (mapfunc,∇mapfunc) in [(:map,:∇map),(:pmap,:∇pmap)]
  @eval function $∇mapfunc(cx, f::F, args::Vararg{Any, N}) where {F, N}
    ys_and_backs = $mapfunc((args...) -> _pullback(cx, f, args...), args...)
    ys = map(first, ys_and_backs)
    function map_back(Δ)
        Δarg = $mapfunc(((_,pb), δ) -> last(pb(δ)), ys_and_backs, Δ)
        (nothing, Δarg)
    end
    return ys, map_back
  end

  @eval @adjoint function $mapfunc(f, args::Union{AbstractArray,Tuple}...)
    $∇mapfunc(__context__, f, args...)
  end
end

accum_sum(xs::AbstractArray{<:Number}; dims = :) = sum(xs, dims = dims)
accum_sum(xs::Number; dims = :) = xs

function unbroadcast(x::AbstractArray, x̄)
  _project(x, x̄)
end

function unbroadcast(x::Number, x̄)
  accum_sum(x̄)
end

@adjoint broadcasted(::typeof(+), xs::Numeric...) =
  broadcast(+, xs...), ȳ -> (nothing, map(x -> unbroadcast(x, ȳ), xs)...)

@adjoint broadcasted(::typeof(-), x::Numeric, y::Numeric) = x .- y,
  Δ -> (nothing, unbroadcast(x, Δ), -unbroadcast(y, Δ))

@adjoint broadcasted(::typeof(*), x::Numeric, y::Numeric) = x.*y,
   Δ -> (nothing, unbroadcast(x, Δ .* conj.(y)), unbroadcast(y, Δ .* conj.(x)))

@adjoint function broadcasted(::typeof(/), x::Numeric, y::Numeric)
  res = x ./ y
  res, Δ -> (nothing, unbroadcast(x, Δ ./ conj.(y)), unbroadcast(y, .-Δ .* conj.(res ./ y)))
end

# Avoid hitting special cases for `Adjoint` etc.
@adjoint broadcasted(::AbstractArrayStyle, f::F, args...) where {F} = broadcast_forward(f, args...)
@adjoint! (b::typeof(broadcast))(f, args...) = _pullback(__context__, broadcasted, f, args...)

# We do this because it ensures type stability so it compiles nicely on the gpu
# The val is needed for some type stability
@inline dual(x, i, ::Val{N}) where {N} = x
@inline dual(x::Real, i, ::Val{N}) where {N} = Dual(x, ntuple(==(i), N))

function dualize(args::Vararg{Any, N}) where {N}
    ds = map(args, ntuple(identity,N)) do x, i
        return dual(x, i, Val(N))
      end
      return ds
end

@inline function dual_function(f::F) where F
    function (args::Vararg{Any,N}) where N
      ds = dualize(args...)
      return f(ds...)
    end
  end

@inline function broadcast_forward(f, args::Vararg{Any,N}) where N
  out = dual_function(f).(args...)
  T = eltype(out)
  T <: Union{Dual, Complex{<:Dual}} || return (out, _ -> nothing)
  _broadcast_forward(T, out, args...)
end

# Real input and real output pullback
@inline function _broadcast_forward(::Type{<:Dual}, out, args::Vararg{Any, N}) where {N}
  valN = Val(N)
  y = broadcast(x -> value(x), out)
  function bc_fwd_back(ȳ)
    dargs = ntuple(valN) do i
      unbroadcast(args[i], broadcast((y1, o1) -> y1 * partials(o1,i), ȳ, out))
    end
    (nothing, nothing, dargs...) # nothings for broadcasted & f
  end
  return y, bc_fwd_back
end

# we need to define this late, so that the genfuncs see lib.jl
# Move using statements out of this file to help with sysimage building
ignore_sig(T) = all(T -> T <: Type, T.parameters)

function edge!(m::IRTools.Meta, edge::Core.MethodInstance)
  m.code.edges === nothing && (m.code.edges = Core.MethodInstance[])
  push!(m.code.edges, edge)
  return
end

function _generate_pullback(ctx, world, f, args...)
    cr_T = Tuple{ZygoteRuleConfig{ctx}, f, args...}
    chain_rrule_f = :chain_rrule

  hascr, cr_edge = has_chain_rrule(cr_T, world)
  hascr && return :($chain_rrule_f(ZygoteRuleConfig(ctx), f, args...))

  # No ChainRule, going to have to work it out.
  T = Tuple{f,args...}
  ignore_sig(T) && return :(f(args...), Pullback{$T}(()))

  g = try
    _generate_pullback_via_decomposition(T, world)
  catch e
    if VERSION < v"1.8"
      # work around Julia bug
      rethrow(CompileError(T,e))
    end
    return :(throw($(CompileError(T,e))))
  end
  g === nothing && return :(f(args...), Pullback{$T}((f,)))
  meta, forw, _ = g
  argnames!(meta, Symbol("#self#"), :ctx, :f, :args)
  forw = varargs!(meta, forw, 3)
  # IRTools.verify(forw)
  forw = slots!(pis!(inlineable!(forw)))
  # be ready to swap to using chainrule if one is declared
  cr_edge != nothing && edge!(meta, cr_edge)
  return update!(meta.code, forw)
end

function _generate_callable_pullback(j::Type{<:Pullback{T}}, world, Δ) where T
  ignore_sig(T) && return :nothing
  g = try
    _generate_pullback_via_decomposition(T, world)
  catch e
    if VERSION < v"1.8"
      # work around Julia bug
      rethrow(CompileError(T,e))
    end
    return :(throw($(CompileError(T,e))))
  end
  if g === nothing
    Δ == Nothing && return :nothing
    return :(error("Non-differentiable function $(repr(j.t[1]))"))
  end
  meta, _, back = g
  argnames!(meta, Symbol("#self#"), :Δ)
  # IRTools.verify(back)
  back = slots!(inlineable!(back))
  return update!(meta.code, back)
end

@generated function _pullback(ctx::AContext, f, args...)
  _generate_pullback(ctx, nothing, f, args...)
end

@generated function (j::Pullback)(Δ)
  _generate_callable_pullback(j, nothing, Δ)
end

end # module
