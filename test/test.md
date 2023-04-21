# Summary of internal tests

*Please run the gendata.jl script in each folder to generate input data at first.*

## Standard tests

The following ten tests are used to test the basic features of the ACFlow toolkit.

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
    * MaxEnt solver (chi2kink algorithm)
    * StochOM solver
    * Two-Lorentzians model
    * Linear mesh
    * Script mode

* **basic/A04**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochOM solver
    * Gaussian model
    * Lorentzian mesh
    * Script mode
    * Fixed error bar

* **basic/A05**
    * Fermionic + Matsubara (Auxiliary Green's function)
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
    * Fermionic + Matsubara (Matrix-valued functions)
    * MaxEnt solver (chi2kink algorithm)
    * Gaussian model + File model
    * Linear mesh
    * Script mode
    * Fixed error bar
    * MaxEntAux algorithm

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
    * MaxEnt solver (chi2kink algorithm)
    * StochAC solver (with constraints)
    * StochSK solver (with constraints)
    * Rise-And-Decay model + Flat model
    * Half-Lorentzian mesh
    * Fixed error bar

## Experimental tests

The following tests are designed to test the newly developed StochPX solver (based on the stochastic pole expansion). Note that the `X` series are for the fermionic systems, the `Y` series are for the bosonic systems, and the `Q` series are for the lattice QCD data.

* **future/X01**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Two gaussian peaks + big gap

* **future/X02**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Broad spectrum, two peaks
    * A clone of **A01**

* **future/X03**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochAC solver
    * StochSK solver
    * StochOM solver
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * With sharp quasiparticle peak

* **future/X04**
    * Fermionic + Matsubara
    * Green's function
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver (with constraints)
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Sharp gap edges
    * A clone of **T03**

* **future/X05**
    * Fermionic + Matsubara
    * Green's function + Self-energy function
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model + Gaussian model
    * Linear mesh + Tangent mesh
    * Standard mode
    * Script mode
    * A clone of **T01**

* **future/X06**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Single off-centered delta peak

* **future/X07**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochAC solver
    * StochSK solver
    * StochOM solver
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Four off-centered delta peaks

* **future/X08**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver (with constraints)
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Four off-centered delta peaks

* **future/X09**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Eight off-centered delta peaks

* **future/Y01**
    * Bosonic + Matsubara (imaginary part = 0)
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Two gaussian peaks + big gap

* **future/Y02**
    * Bosonic (symm) + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Two gaussian peaks + big gap

* **future/Y03**
    * Bosonic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Two gaussian peaks + big gap

* **future/Y04**
    * Bosonic + Matsubara (imaginary part = 0)
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Four off-centered delta peaks

* **future/Y05**
    * Bosonic (symm) + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Four off-centered delta peaks

* **future/Y06**
    * Bosonic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Four off-centered delta peaks

* **future/Y07**
    * Bosonic (symm) + Matsubara
    * Current-current correlation function
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Half-Lorentzian mesh
    * Fixed error bar

* **future/Y08**
    * Bosonic (symm) + Matsubara
    * Lindhard function on square lattice
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * linear mesh
    * Fixed error bar
    * Script mode

* **future/Y09**
    * Bosonic (symm) + Matsubara
    * Green's function
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver (with constraints)
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Sharp gap edges + broad tail

* **Q01**
* **Q02**
* **Q03**
