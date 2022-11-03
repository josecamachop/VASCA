
# Variable-Selection ASCA (VASCA)

This repository contains the code for reproducing the simulation results of the paper Camacho J, Vitale R, Morales-Jimenez D and Gómez-Llorente C. Variable-Selection ANOVA Simultaneous Component Analysis. Submitted to Bioinformatics. 2022 

This repository requires the MEDA toolbox v1.3 at https://github.com/josecamachop/MEDA-Toolbox/tree/v1.3

Contact person: José Camacho (josecamacho@ugr.es)

Last modification of this document: 3/Nov/22


## Organization

The code is organized in a number of Matlab scripts:

- The scripts named "RunExample#" perform the complete computation of the corresponding simulated example in the paper.
- The script named "RunWheat" performs the motivating example with ASCA.

The present repository is organized as follows:

- README.md					This document.

- RunExample1.m				Contains the code for the example with non-significant relationships.
- RunExample1b.m			Contains the code for the example with non-significant relationships in two factors with several levels.
- RunExample2.m				Contains the code for the example with significant one-to-one relationships.
- RunExample2b.m			Contains the code for the example with significant one-to-one relationships and a small bias.
- RunExample3.m				Contains the code for the example with a multivariate relationship.
- RunExample4.m				Contains the code for the example with a multivariate relationship in two factors with several levels.
- RunWheat.m				Contains the code for the ASCA motivating example.

- html/						Contains the web reports obtained with command "publish".
- Fig/						Contains the figures of the paper.

- example1.mat				Contains the workspace after the computation of Example 1 for simuleMV as the simulation engine.
- example1_randn.mat		Contains the workspace after the computation of Example 1 for randn as the simulation engine.
- example1b.mat				Contains the workspace after the computation of Example 1b.
- example2.mat				Contains the workspace after the computation of Example 2.
- example2b.mat				Contains the workspace after the computation of Example 2b.
- example3.mat				Contains the workspace after the computation of Example 3.
- example4.mat				Contains the workspace after the computation of Example 4.
- example5.mat				Contains the workspace after the computation of Example 5.
- wheat						Contains the wheat data.

- parglm_genes				Contains our interpretation of the code of ASCA-genes.

