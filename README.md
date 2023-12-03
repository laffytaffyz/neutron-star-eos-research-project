# Resolving the Constraints Imposed by Chiral Effective Field Theory and Perturbative Quantum Chromodynamics on the Neutron Star Equation of State
### Tiffany Zhang<sup>1,2</sup>, James M. Lattimer<sup>2</sup>
#### <sup>1</sup>Great Neck South High School, <sup>2</sup>Stony Brook University

**Abstract** <br>
The structure and composition of neutron stars (NS) are supported by observational data for outer regions but are not fully understood for deeper layers, especially the inner core. This knowledge gap in the inner core is reflected in the equation of state (EOS), the function between pressure and energy density that governs all NS properties. Chiral effective field theory (χEFT) accounts for the less dense outer regions, and perturbative quantum chromodynamics (pQCD) accounts for the extreme NS densities, beyond the maximal densities in the inner core. However, no observational data supports pQCD constraining the lower-density EOS, so other possible high-density limits were sought. Moreover, the constraints on the inner-core EOS have been dependent on the specific theoretical approach. Hence, this study aims to define inner-core EOSs with other high-density limits in a model-agnostic way. Varying the high-density limit and interpolating with the low-density χEFT limit based on fundamental thermodynamic constraints, the generated EOSs inconsistent with observational NS mass and radii data were removed. The novel EOSs identified in this study produced an accurate and precise mean radius 12.238±0.541 km for the theorized NS maximum mass. These EOSs yielded similar properties to χEFT extrapolations and pQCD–χEFT interpolations, resolving between these two opposing models while challenging the existence of deconfined quark cores. This finding suggests the possibility of the EOS blending the two models, which further theoretical and empirical analysis can validate to determine the inner core’s composition and to establish the neutron star equation of state.
<br><br>
[View paper here](https://docs.google.com/document/d/1v-nBOz4MStbWZRciXdfEbUac13yA1kjHQwrRJMa4FG0/edit?usp=sharing)

## File Descriptions
- **[Research notebook](/research%20notebook.ipynb)**: processes all the data and generates all of the figures
- **[Data generation program](/data_generation.py)**: generates all the EOS data
- **TOVsolver files (credit: [Anton Motornenko](https://github.com/amotornenko/TOVsolver))**: solves the Tolman-Oppenheimer-Volkoff equations given the neutron star equation of state
  - [TOV.py](/tov.py): class for solving for TOV equations
  - [constants.py](/constants.py): class containing physical constants
- **[Data folder](/data/)**: all theoretical and observational data used in this study (see [data read me](/data/README.md))


