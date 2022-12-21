# Summary of internal tests

*Please run the gendata.jl script in each folder to generate input data at first.*

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
