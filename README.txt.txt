This repository includes R code for reproducing Figure 2 and 3, Table 1 to 3, Web Table 1 to 6, and Web Figure 1 in the Biometrics paper "Planning Three-Level Cluster Randomized Trials to Detect Treatment Effect Heterogeneity" (under review).

For questions or comments about the code, please contact Zizhong Tian at <zqt5121@psu.edu>.

I. Supporting Files: These supporting files are sourced in the corresponding main files that reproduce the simulation tables and illustrative plots in the main manuscript as well as the supplementary web appendix.

Folder: functions
1) L1_sim_functions.R = functions about sample size/power formula and data generation under the scenario of level-1 (individual-level) randomization;
2) L2_sim_functions.R = functions about sample size/power formula and data generation under the scenario of level-2 (subcluster-level) randomization;
3) L3_sim_functions.R = functions about sample size/power formula and data generation under the scenario of level-3 (cluster-level) randomization;
4) L1_plot_functions.R = functions used for data application under the scenario of level-1 (individual-level) randomization;
5) L2_plot_functions.R = functions used for data application under the scenario of level-2 (individual-level) randomization;
6) L3_plot_functions.R = functions used for data application under the scenario of level-3 (individual-level) randomization;

II. Main Files: These main files are used to reproduce the simulation results and illustrative plots in the main manuscript as well as the web appendix.

Folder: simulations
7) L1_HTE.R = reproduce the simulation table for testing HTE under level-1 randomization;
8) L1_OTE.R = reproduce the simulation table for testing ATE under level-1 randomization;
9) L1_OTE_unadj.R = reproduce the simulation table for testing unadjusted ATE under level-1 randomization;
10) L2_HTE.R = reproduce the simulation table for testing HTE under level-2 randomization;
11) L2_OTE.R = reproduce the simulation table for testing ATE under level-2 randomization;
12) L2_OTE_unadj.R = reproduce the simulation table for testing unadjusted ATE under level-1 randomization;
13) L3_HTE.R = reproduce the simulation table for testing HTE under level-3 randomization;
14) L3_OTE.R = reproduce the simulation table for testing ATE under level-3 randomization;
15) L3_OTE_unadj.R = reproduce the simulation table for testing unadjusted ATE under level-1 randomization;
16) xtable.R = generate the LaTeX editing version of all table contents;
Folder: plots
17) Remark_1_Level2.R = reproduce the illustrative variance plot regarding Remark 1 (level-2 randomization);
18) Remark_1_Level3.R = reproduce the illustrative variance plot regarding Remark 1 (level-3 randomization);
19) HALI_example_BaselineScore.R = reproduce the power analysis plot in the data application section;

III. Software 

Analyses were conducted with R, version 4.0.3 (https://www.r-project.org/)
The calculations used R packages nlme (version 3.1-143).

IV. R commands for the installation of R packages 

install.packages(c("nlme")) 