# Programming bistability in geometrically perturbed mechanical metamaterials

<img src = "https://github.com/Hougeichao/Programming-bistability-in-geometrically-perturbed-mechanical-metamaterials/blob/main/Figure/Design-exploration.png" />

If you use this code in your work, please cite the following paper:

>Peng, Y., Niloy, I., Kam, M., Celli, P., & Plucinsky, P. (2024). "[Programming bistability in geometrically perturbed mechanical metamaterials](https://doi-org.libproxy2.usc.edu/10.1103/PhysRevApplied.22.014073)". Physical Review Applied, 22(1), 014073.

===============================================================

How to use:

Please run `Main.m` to generate homogeneous bistable kirigami designs with prescribed Bravais Lattice vectors. Figure 1 and Figure 2 show the unit cell in two stable states. Figure 3 and Figure 4 illustrate the 2-by-2 designs in the reference state and deformed state respectively. 

Please run `Postresult.m` to generate energy results. Figure 5 depicts the relationship between the elastic energy and the characteristic stretch $\lambda$. Figure 6 illustrates the relationship between the energy and the actuation angle $\xi$.

Main programs:

* `Main.m` is the main MATLAB code to implement optimization.
* `Para.m` parameterizes the kirigami geometry by the design variables and plots the figures in two states by setting p = 1 (with prescribed Bravais lattice vectors).
* `Obj.m` contains the objective function used for the simulation such as the target energy barrier,and target stiffnesses at two stable states.
* `Constraint.m` contains the inequality constraints used for the optimization.
* `Energy_spr.m` and `EnergyBarrier.m` construct the functions for the elastic energy and energy barrier respectively.  
* `Rot.m` represents the 2D rotation tensor.
* `ccc.m` is used for the clearance.

Remarks:
