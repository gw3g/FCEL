# Fully Coherent Energy Loss

All FCEL particulars are bundled in the header file ``fcel.h``, which includes the function ``channel(int i)`` the partonic 2->2 subprocess to be selected. This will prepare the appropriate SU(3) irreps of the final state and their colour probabilities as global quantities. The argument of ``channel(..)``  corresponds with the QCD processes

1. `i=1` g, g -> g, g
2. `i=2` q, g -> q, g
3. `i=3` g, q -> q, g
4. `i=4` g, g -> q, qbar

Here we consider two applications: i) LHC physics, ii) atmospheric neutrinos.

## nuclear mod. factor @ LHC ``pA``

Compute R<sub>pA</sub> for light [[1](https://arxiv.org/abs/2003.06337),[2](https://arxiv.org/abs/2003.01987)] and heavy mesons production, due solely to FCEL. The LHC parameters are the collisional energy, the nucleon number, and the strong coupling. Respectively, their default parameters are
```
double SQRTS = 8.16;  // TeV
double A = 208.;      // Pb
double alpha_s = 0.5; // coupling
```
Data is saved under ``pA/out/``, where file names indicate the partonic subprocess ``{ab;cd}`` (for a,b -> c,d) and the fixed transverse momentum (pT) or rapidity (y).

* ``R_reps(double pT)`` generates a table of R(y)
* ``R_scan_y(double pT)`` generates a table of R(y)
* ``R_scan_pT(double y)`` generates a table of R(pT)

Uncertainty bands are generated by the Hessian method, which is applied to the parameter list:
```
#define PARAMS 5  // q0, xi, z, n,  m_Q
double dS[PARAMS] = {.02,.25,.2,1.,0.5};
double  S[PARAMS] = {.07,.50,.8,4.,1.5};
```


## atmospheric neutrino flux ``nu/``

Conventional and prompt ν's are produced by cosmic rays (CRs) colliding with Air nuclei (<A>=14.5). Using the approximate Z-moment solution to the cascade equation, the flux on the earths surface can be expressed by a convolution of the charm cross-section with the CR flux. The cross-section dσ/dx<sub>F</sub> is defined in ``sigma.h``, as a functions of x<sub>F</sub> and the proton energy E<sub>p</sub>. Three CR fluxes are included in ``crflux.h``:

* knee-spectrum
* H3a
* Global-spline fit

