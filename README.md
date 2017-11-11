# Finite Element Mesh
The project consists of the source code of a 2D geometry of an
electrical machine. The [gmsh] (http://gmsh.info/) .geo file can be used to mesh the domains of the machine geometry. The **gmsh** generated mesh can be used for the numerical analysis, such as the finite element method.

## Surface mounted permanent magnet machine
+ `./50kW_PMM/pmsm_machine.geo:` consists of the source code which generates one quarter
of the cross-section of a surface mounted permanent magnet (SMPM) machine. The air-gap region is divided
into three layers; the middle layer could be utilized for the rotation.
