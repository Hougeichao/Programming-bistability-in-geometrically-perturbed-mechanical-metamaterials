### Code for optimization part of 'Programming bistability in geometrically perturbed mechanical meta materials' ###

Please run 'Main.m' to generate homogeneous bistable kirigami designs with prescribed Bravais Lattice vectors. Figure 1 and Figure 2 show the unit cell in two stable states. Figure 3 and Figure 4 illustrate the 2-by-2 designs in the reference plane and deformed plane respectively. 

Run **Postresult.m** to generate energy results: Figure 5 and Figure 6. Figure 5 depicts the relationship between the elastic energy and the characteristic stretch \lambda. Figure 6 illustrates the relationship between the energy and the actuation angle \xi.

**Main.m** is the main MATLAB code to implement optimization.

**Para.m** contains Matlab code for parameterizing all the vectors in kirigami design by the design variables. The theoretical stiffnesses are calculated in this function as well. Also, it can be used to create the geometry of the patterns in two states by giving the eight design variables and Bravais Lattice vectors and defining the input parameter p = 1. The coordinates of one unit cell in two states are stored as the output of Para function (Coor_ref and Coor_def)

**Obj.m** contains the objective function used for the simulation such as target energy barrier, target stiffnesses at two stable states, etc...

**Constraint.m** contains the inequality constraints used for the optimization.

**Energy_spr.m** and **EnergyBarrier.m** construct the functions for the elastic energy and energy barrier respectively. 

**Rot.m** represents the rotation tensor. **ccc.m** is used for the clearance.
