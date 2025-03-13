# Programming bistability in geometrically perturbed mechanical metamaterials

This repository provides a versatile design framework for designing and optimizing bistable kirigami metamaterials with target shapes and mechanical properties. 

The framework allows users to design custom kirigami patterns by prescribing shape-change lattice vectors in two stable states ($\boldsymbol{\ell}_1^R, \boldsymbol{\ell}_2^R$ in the first stable state and $\boldsymbol{\ell}_1^D, \boldsymbol{\ell}_2^D$ in the second), enabling the realization of homogeneous shape-morphing kirigami metamaterials. 

_<ins>**Please use the code to design your own Kirigami pattern!** </ins>_


<img src = "https://github.com/Hougeichao/Programming-bistability-in-geometrically-perturbed-mechanical-metamaterials/blob/main/Figure/Github_Fig1.png" width="1000"/>

We highlight the richness of the design space by producing kirigami patterns that exhibit a variety of axial and shearing shape changes.

<img src = "https://github.com/Hougeichao/Programming-bistability-in-geometrically-perturbed-mechanical-metamaterials/blob/main/Figure/Design-exploration.png" width="1000"/>

## Paper
This repository contains the code developed for the paper below. Feel free to try this code to create your own planar Kirigami pattern. 

>Peng, Y., Niloy, I., Kam, M., Celli, P., & Plucinsky, P. (2024). "[Programming bistability in geometrically perturbed mechanical metamaterials](https://doi-org.libproxy2.usc.edu/10.1103/PhysRevApplied.22.014073)". Physical Review Applied, 22(1), 014073.



## Getting Started

### Prerequisites
* MATLAB (tested on R2023a and later versions)
* Optimization Toolbox (for `FMINCON`)

### Usage 

1. Main Script:

- Run `Main.m` to generate homogeneous bistable kirigami designs with prescribed Bravais lattice vectors.
  - **Figure 1 & 2**: Display the unit cell in the two stable states.
  - **Figure 3 & 4**: Illustrate the 2-by-2 designs in the reference and deformed states, respectively.

2. Post-Processing:

- Run `Postresult.m` to generate energy-related results.

  - **Figure 5**: Depicts the relationship between elastic energy and the characteristic stretch $\lambda$
  - **Figure 6**: Illustrates the relationship between the energy and the actuation angle $\xi$.


### Key Functions

* `Main.m`: The main MATLAB code to implement optimization.
* `Para.m`: Parameterizes the kirigami geometry by the design variables and plots the figures in two states by setting p = 1 (with prescribed Bravais lattice vectors).
* `Obj.m`: Defines the objective function used for the simulation such as the target energy barrier,and target stiffnesses at two stable states.
* `Constraint.m`: Specifies the inequality constraints used for the optimization.
* `Energy_spr.m` and `EnergyBarrier.m`: Construct the functions for the elastic energy and energy barrier respectively.  
* `Rot.m`: Represents the 2D rotation tensor.
* `ccc.m`: Is used for the clearance.

### Customization and Optimization

* Feel free to experiment with different shape-change lattice vectors and energy barriers to customize your kirigami pattern.
* If the optimization does not yield a desirable pattern, consider the following:
  * Adjusting fmincon parameters (e.g., MaxFunctionEvaluations, MaxIterations, ConstraintTolerance) to refine the optimization process.
  * Modifying the prescribed Bravais lattice vectors to explore different configurations.
  * Tuning the objective function to better match your design goals.
  * Changing the initial point.
* **Note**: Some extreme shape changes, while mathematically valid, may not yield physically realizable kirigami patterns. Users are encouraged to explore a wide range of configurations.

### How to cite
If you use this code in your work, please cite the following [paper](https://doi-org.libproxy2.usc.edu/10.1103/PhysRevApplied.22.014073)

```bibtex
@article{peng2024programming,
  title={Programming bistability in geometrically perturbed mechanical metamaterials},
  author={Peng, Yingchao and Niloy, Imtiar and Kam, Megan and Celli, Paolo and Plucinsky, Paul},
  journal={Physical Review Applied},
  volume={22},
  number={1},
  pages={014073},
  year={2024},
  publisher={APS}
}
```
