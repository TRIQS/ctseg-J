---
title: 'CTSEG: A segment picture quantum impurity solver based on TRIQS'
tags:
  - Quantum Monte Carlo
  - impurity solver
authors:
  - name: Nikita Kavokine
    orcid: 0000-0002-8037-7996
    affiliation: "1, 2" 
  - name: Hao Lu
    orcid: 0009-0000-4581-9544
    affiliation: 2
  - name: Thomas Ayral
    orcid: 0000-0003-0960-4065
    affiliation: 3
  - name: Michel Ferrero
    orcid: 0000-0003-1882-2881
    affiliation: "4, 5"
  - name: Nils Wentzell
    orcid: 0000-0003-3613-007X
    affiliation: 1
  - name: Olivier Parcollet
    orcid: 0000-0002-0389-2660
    affiliation: 1
affiliations:
 - name: Center for Computational Quantum Physics, Flatiron Institute, 162 5th Avenue, NY 10010, New York, USA
   index: 1
 - name: Max Planck Institute for Polymer Research, Ackermannweg 10, 55128 Mainz, Germany
   index: 2
 - name: Eviden Quantum Laboratory, 78340 Les Clayes-sous-Bois, France
   index: 3
 - name: CPHT, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris, 91128 Palaiseau, France
   index: 4
 - name: Collège de France, 11 place Marcelin Berthelot, 75005 Paris, France
   index: 5
date: 25 May 2024
bibliography: paper.bib
---

# Summary

Electron-electron interactions are a key determinant of the electrical and 
optical response properties in solid-state materials. Yet, the many-electron 
problem is outstandingly difficult, and can be tackled analytically only when 
interactions are weak, or in rare exactly solvable cases. One of the popular 
numerical schemes for addressing strongly interacting systems is dynamical 
mean-field theory (DMFT) [@georges1996]. In DMFT, the many-body problem formulated on a
crystal lattice is self-consistently mapped onto a single atom (or impurity)
immersed in an effective environment (or bath). The remaining computational 
task is then the solution of the impurity problem, which can be carried out 
through various quantum Monte Carlo algorithms [@gull2011]. Here, we present an implementation
of the continous time hybridization expansion algorithm in the segment picture (`CTSEG`)
based on `TRIQS`, a comprehensive library for the numerical investigation of interacting 
quantum systems [@TRIQS2015]. 

# Statement of need

The Monte Carlo algorithms for quantum impurity problems are 
based on stochastically exploring the terms in the perturbative expansion of the solution 
around an exactly solvable limit. Continuous time hybridization expansion algorithms -- chief of which `CTHYB` -- involve expanding around the limit of an isolated atom [@gull2011]. 
Currently, there exist implementations of `CTHYB` within three different libraries: `ALPS` [@ALPS2018], `w2dynamics` [@w2dynamics2019] and `TRIQS` [@CTHYB2016].

However, a simpler and potentially faster version of the `CTHYB` algorithm, 
called `CTSEG`, can be used under the restriction of density-density
interactions on the impurity [@werner2006a;@werner2006b;@haule2007]. `CTSEG` can be further generalized to allow for time-dependent[@werner2007;@werner2010] and 
spin-spin interactions [@otsuki2013;@steiner2015]: see `@werner2016` for a review. To our knowledge, there exists so far one published implementation of `CTSEG` based on ALPS [@ALPS-CTSEG], but it does not allow for spin-spin interactions. 

Our `CTSEG` solver is about twice as fast as `TRIQS-CTHYB` for a single orbital problem, and has
better scaling with the number of orbitals (40 times faster in our 5 orbital test case, see Fig. 1a). 
`CTSEG` has already allowed us to obtain the first numerically-exact solution of the 
quantum Heisenberg spin glass [@kavokine2024]. 

![**a**. Running time comparison between the TRIQS implementations of CTSEG and CTHYB. The test system is a multi-orbital impurity at half-filling and inverse temperature $\beta = 20$. The Coulomb repulsion is $U = 2$ for two electrons on the same orbital and $U' = 1$ for two electrons on different orbitals. The hybridization is diagonal and identical for all orbitals: $\Delta(\omega) = 1/(\omega - 0.3)$. **b**. Spin-spin correlation function $\chi(\tau) = \langle \mathbf{S}(\tau) \cdot \mathbf{S}(0) \rangle$ of the $t-J-U$ model studied by Dumitrescu et al., obtained using CTSEG at inverse temperature $\beta = 300$ and different values of doping $p$. At long times $\chi(\tau) \sim 1/\tau^{\theta}$, with $\theta = 1$ at the QCP. Inset: exponent $\theta$ as a function of doping $p$. The QCP is located at $p \approx 0.16$.](figure_JOSS.pdf){width=80%}

# Example of use

As a further illustration of our solver's performance, we apply it to the fully connected $t-J-U$ model
studied by @dumitrescu2022. At half-filling, the model forms a spin glass phase, which melts into 
a metal at a doping-induced quantum critical point (QCP). Dumitrescu et al. 
obtained solutions at inverse temperatures up to $\beta = 65$, limited by the fermionic sign problem 
of their interaction expansion solver. The hybridization expansion carried out by `CTSEG` is 
sign-problem-free for the $t-J-U$ model, allowing us to reach $\beta = 300$ and to obtain a 
more accurate localization of the QCP at doping $p \approx 0.16$ (Fig. 1b).

# Acknowledgements

We thank Alexander Hampel for testing the code and providing valuable feedback.
The Flatiron Institute is a division of the Simons Foundation. 

# References
