# RGD_CT_Work
This repository contains the analysis code and related data files developed for the Color Transparency experiment which ran in fall 2023 by the CLAS12 Collaboration of Jefferson Labs


The code file contains all analysis code I used.

The main series of analysis data is the rgd_ct*.cc code. They take list files at the input and produce lengthy root files of histograms and graphs
The latest verison I used is rgd_ct2_6.cc and this code has comments added to help with navigation. These codes are ment to be run on ifarm. The codes with tree in the name produce root tree files.
The latest code I produced are the ctest and they contain some of the fiducial graphs. The cleaest data file is the one call bcolloquial, as these were the files I used to produce the graphs for my colloquial talk in August of 2025.

the code fitsaps_5.cc is the latest version of the fitting code I would use to apply the Briet-Winger fits and get the Rho0 yield of each Q2 bin. This code takes the root files produced by rgd_ct*.cc and produces a root file of the fitted graphs.
the latest version of the data is Invarient_mass_plots_for_Mathieu.cc

the code graph_cusn_better.cc is the latest version of the code for calculating nuclear transparency and producting the T vs Q2 graph. This code takes the yields produced by fitsaps_5.cc.
The latest data file is Nuclear_transparency_plots_for_Mathieu.cc
