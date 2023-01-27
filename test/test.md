# Summary of internal tests

*Please run the gendata.jl script in each folder to generate input data at first.*

## Standard tests

The following ten tests are used to test the basic features of the ACFlow toolkit.

* **A01**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink or bryan algorithm)
    * StochAC solver
    * StochSK solver
    * StochOM solver
    * Flat model
    * Linear mesh
    * Fixed error bar

* **A02**
    * Bosonic (symm) + Matsubara
    * MaxEnt solver (chi2kink or bryan algorithm)
    * StochAC solver
    * StochSK solver
    * StochOM solver
    * Flat model
    * Half-Lorentzian mesh
    * Fixed error bar

* **A03**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochOM solver
    * Two-Lorentzians model
    * Linear mesh
    * Script mode

* **A04**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochOM solver
    * Gaussian model
    * Lorentzian mesh
    * Script mode
    * Fixed error bar

* **A05**
    * Fermionic + Matsubara (Auxiliary Green's function)
    * MaxEnt solver (chi2kink algorithm + preblur)
    * StochOM solver
    * Gaussian model
    * Tangent mesh
    * Script mode

* **A06**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm + preblur)
    * StochOM solver
    * Gaussian model
    * Lorentzian mesh
    * Fixed error bar

* **A07**
    * Bosonic (symm) + Matsubara
    * MaxEnt solver (chi2kink algorithm + preblur)
    * StochOM solver
    * Rise-And-Decay model
    * Half-Lorentzian mesh
    * Fixed error bar

* **A08**
    * Fermionic + Matsubara (Matrix-valued functions)
    * MaxEnt solver (chi2kink algorithm)
    * Gaussian model + File model
    * Linear mesh
    * Script mode
    * Fixed error bar

* **A09**
    * Fermionic + Imaginary time
    * MaxEnt solver (bryan algorithm)
    * StochAC solver (with constraints)
    * StochSK solver (with constraints)
    * Two-Gaussians model + Flat model
    * Linear mesh
    * Fixed error bar

* **A10**
    * Bosonic + Imaginary time
    * MaxEnt solver (chi2kink algorithm)
    * StochAC solver (with constraints)
    * StochSK solver (with constraints)
    * Rise-And-Decay model + Flat model
    * Half-Lorentzian mesh
    * Fixed error bar

## Experimental tests

The following tests are designed to test the newly developed StochPX solver (based on the stochastic pole expansion).

* **X01**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Two gaussian peaks + big gap

* **X02**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * A clone of **A01**
    * Broad spectrum, two peaks

* **X03**
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

* **X04**
    * Fermionic + Matsubara + Imaginary time
    * Green's function
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver (with constraints)
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Standard mode
    * A clone of **T03**
    * Sharp gap edges

* **X05**
    * Fermionic + Matsubara
    * Green's function + Self-energy function
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model + Gaussian model
    * Linear mesh + Tangent mesh
    * Standard mode
    * Script mode
    * A clone of **T01**

* **X06**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Single off-centered delta peak

* **X07**
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

* **X08**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver (with constraints)
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Four off-centered delta peaks

* **X09**
    * Fermionic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Eight off-centered delta peaks

* **Y01**
    * Bosonic + Matsubara (imaginary part = 0)
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Two gaussian peaks + big gap

* **Y02**
    * Bosonic (symm) + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Two gaussian peaks + big gap

* **Y03**
    * Bosonic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Two gaussian peaks + big gap

* **Y04**
    * Bosonic + Matsubara (imaginary part = 0)
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Four off-centered delta peaks

* **Y05**
    * Bosonic (symm) + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Four off-centered delta peaks

* **Y06**
    * Bosonic + Matsubara
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Linear mesh
    * Fixed error bar
    * Four off-centered delta peaks

* **Y07**
    * Bosonic (symm) + Matsubara
    * Current-current correlation function
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * Half-Lorentzian mesh
    * Fixed error bar

* **Y08**
    * Bosonic (symm) + Matsubara
    * Lindhard function on square lattice
    * MaxEnt solver (chi2kink algorithm)
    * StochPX solver
    * Flat model
    * linear mesh
    * Fixed error bar
    * Script mode
