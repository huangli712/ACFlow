# Summary of internal tests

*Please run the gendata.jl script in each folder to generate input data at first.*

## Standard tests

The following twelve tests are used to test the basic features of the ACFlow toolkit (for `MaxEnt`, `StochAC`, `StochSK`, and `StochOM` solvers only).

* **basic/A01**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink or bryan algorithm)
    * StochAC solver
    * StochSK solver
    * StochOM solver
    * Flat model
    * Linear mesh
    * Fixed error bar

* **basic/A02**
    * Bosonic (symm) + Matsubara
    * MaxEnt solver (chi2kink or bryan algorithm)
    * StochAC solver
    * StochSK solver
    * StochOM solver
    * Flat model
    * Half-Lorentzian mesh
    * Fixed error bar

* **basic/A03**
    * Fermionic + Matsubara
    * Realistic self-energy function
    * MaxEnt solver (chi2kink algorithm)
    * StochOM solver
    * Two-Lorentzians model
    * Linear mesh
    * Script mode

* **basic/A04**
    * Fermionic + Matsubara
    * Realistic Green's function
    * MaxEnt solver (chi2kink algorithm)
    * StochOM solver
    * Gaussian model
    * Lorentzian mesh
    * Script mode
    * Fixed error bar

* **basic/A05**
    * Fermionic + Matsubara (Auxiliary Green's function)
    * Realistic self-energy function
    * MaxEnt solver (chi2kink algorithm + preblur)
    * StochOM solver
    * Gaussian model
    * Tangent mesh
    * Script mode

* **basic/A06**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm + preblur)
    * StochOM solver
    * Gaussian model
    * Lorentzian mesh
    * Fixed error bar

* **basic/A07**
    * Bosonic (symm) + Matsubara
    * MaxEnt solver (chi2kink algorithm + preblur)
    * StochOM solver
    * Rise-And-Decay model
    * Half-Lorentzian mesh
    * Fixed error bar

* **basic/A08**
    * Fermionic + Matsubara
    * Matrix-valued Green's functions
    * MaxEnt solver (chi2kink algorithm)
    * Gaussian model + File model
    * Linear mesh
    * Script mode
    * Fixed error bar
    * MaxEntAux algorithm + Direct algorithm

* **basic/A09**
    * Fermionic + Imaginary time
    * MaxEnt solver (bryan algorithm)
    * StochAC solver (with constraints)
    * StochSK solver (with constraints)
    * Two-Gaussians model + Flat model
    * Linear mesh
    * Fixed error bar

* **basic/A10**
    * Bosonic + Imaginary time
    * Current-current correlation function
    * MaxEnt solver (chi2kink algorithm)
    * StochAC solver (with constraints)
    * StochSK solver (with constraints)
    * Rise-And-Decay model + Flat model
    * Half-Lorentzian mesh
    * Fixed error bar

* **basic/A11**
    * Fermionic + Matsubara (incomplete data)
    * MaxEnt solver (chi2kink or bryan algorithm)
    * StochAC solver
    * StochSK solver
    * StochOM solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * A clone of **A01**

* **basic/A12**
    * Bosonic (symm) + Matsubara (incomplete data)
    * MaxEnt solver (chi2kink or bryan algorithm)
    * StochAC solver
    * StochSK solver
    * StochOM solver
    * Flat model
    * Half-Lorentzian mesh
    * Fixed error bar
    * A clone of **A02**

## Experimental tests 1

The following tests are designed to test the newly developed `NevanAC` solver (based on the `Nevanlinna analytical continuation` method) only.

* **nac/N01**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * NevanAC solver + Hardy optimization
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Two gaussian peaks

* **nac/N02**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * NevanAC solver + Hardy optimization
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Two gaussian peaks + big gap

* **nac/N03**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * NevanAC solver + no Hardy optimization
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Three delta-like peaks

## Experimental tests 2

The following tests are designed to test the newly developed `StochPX` solver (based on the `stochastic pole expansion` method) only. Note that the `X` series are for the fermionic systems, the `Y` series are for the bosonic systems, and the `Z` series are for the matrix-valued Green's functions.

* **pole/X01**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Two gaussian peaks + big gap

* **pole/X02**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Broad spectrum, two peaks
    * A clone of **A01**

* **pole/X03**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * With sharp quasiparticle peak

* **pole/X04**
    * Fermionic + Matsubara
    * Green's function
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver (with constraints)
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Sharp gap edges
    * A clone of **T03**

* **pole/X05**
    * Fermionic + Matsubara
    * Green's function + Self-energy function
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model + Gaussian model
    * Linear mesh + Tangent mesh
    * Standard mode
    * Script mode
    * A clone of **T01**

* **pole/X06**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Single off-centered delta peak

* **pole/X07**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Four off-centered delta peaks

* **pole/X08**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver (with constraints)
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Four off-centered delta peaks

* **pole/X09**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Eight off-centered delta peaks

* **pole/Y01**
    * Bosonic + Matsubara (imaginary part = 0)
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Two gaussian peaks + big gap

* **pole/Y02**
    * Bosonic (symm) + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Two gaussian peaks + big gap

* **pole/Y03**
    * Bosonic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Two gaussian peaks + big gap

* **pole/Y04**
    * Bosonic + Matsubara (imaginary part = 0)
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Four off-centered delta peaks

* **pole/Y05**
    * Bosonic (symm) + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Four off-centered delta peaks

* **pole/Y06**
    * Bosonic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Four off-centered delta peaks

* **pole/Y07**
    * Bosonic (symm) + Matsubara
    * Current-current correlation function
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Half-Lorentzian mesh
    * Fixed error bar

* **pole/Y08**
    * Bosonic (symm) + Matsubara
    * Lindhard function on square lattice
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * linear mesh
    * Fixed error bar
    * Script mode

* **pole/Y09**
    * Bosonic (symm) + Matsubara
    * Green's function
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver (with constraints)
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Sharp gap edges + broad tail

* **pole/Z01**
    * Fermionic + Matsubara
    * Matrix-valued Green's functions
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver + self-adaptive mesh
    * Gaussian model + File model
    * Linear mesh
    * Script mode
    * Fixed error bar
    * Two gaussian peaks
    * MaxEntAux algorithm + Direct algorithm

* **pole/Z02**
    * Fermionic + Matsubara
    * Matrix-valued Green's functions
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver + self-adaptive mesh
    * Gaussian model + File model
    * Linear mesh
    * Script mode
    * Fixed error bar
    * Four delta peaks
    * MaxEntAux algorithm + Direct algorithm

## Experimental tests 3

The following tests are designed to test the newly developed `BarRat` solver (based on the `Barycentric rational function approximation` method) only.

* **rfa/R01**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * BarRat solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Two lorentzian peaks + one gaussian peak

* **rfa/R02**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * BarRat solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Two off-centered delta peaks

* **rfa/R03**
    * Fermionic + Matsubara
    * Matrix-valued Green's functions (offdiagonal part)
    * MaxEnt solver (chi2kink algorithm)
    * BarRat solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Three gaussian peaks

* **rfa/R04**
    * Fermionic + Matsubara
    * Matrix-valued Green's functions (offdiagonal part)
    * MaxEnt solver (chi2kink algorithm)
    * BarRat solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Six off-centered delta peaks

* **rfa/R05**
    * Bosonic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * BarRat solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Two gaussian peaks

## LQCD tests

The following tests are designed to test the newly developed `StochPX` solver (based on the `stochastic pole expansion` method) for analytical continuation of the lattice QCD datasets.

* **lqcd/Q01**
    * Bosonic (symm) + Matsubara
    * Breit-Wigner model
    * MaxEnt solver
    * StochPX solver
    * Single broad peak

* **lqcd/Q02**
    * Bosonic (symm) + Matsubara
    * Breit-Wigner model
    * MaxEnt solver
    * StochPX solver
    * Two broad peaks

* **lqcd/Q03**
    * Bosonic (symm) + Matsubara
    * Hadron spectral function
    * MaxEnt solver
    * StochPX solver
    * Resonance peak + continuum function

* **lqcd/Q04**
    * Bosonic (symm) + Matsubara
    * Hadron spectral function
    * MaxEnt solver
    * StochPX solver
    * Resonance peak only

* **lqcd/Q05**
    * Bosonic (symm) + Matsubara
    * Hadron spectral function
    * MaxEnt solver
    * StochPX solver
    * Continuum function only

* **lqcd/Q06**
    * Bosonic (symm) + Matsubara
    * Gaussian mixture model
    * MaxEnt solver
    * StochPX solver + self-adaptive mesh
    * Three separated peaks

* **lqcd/Q07**
    * Bosonic (symm) + Matsubara
    * Bottomonium spectral functions
    * MaxEnt solver
    * StochPX solver + self-adaptive mesh
    * Very large gap + multiple peaks at high energy

* **lqcd/Q08**
    * Bosonic (symm) + Matsubara
    * Charmonium spectral functions (``\eta_c`` channel)
    * MaxEnt solver
    * StochPX solver
    * Very broad spectrum + multiple peaks at high energy
