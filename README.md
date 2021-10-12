Resting State Network Analysis with Asymptotical Surprise and PACO
================
**Repo Manager:** Marisa C. Ross, PhD <br />
**Last updated:** October 12, 2021

## Overview

This project uses the PGC-ENIGMA PTSD resting-state dataset and resolution-limit-free community detection via PACO to investigate differences in functional network topology in PTSD. 

**How to Use This Repo:**

1.  All necessary functions for network analysis are called from the [master.m](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/blob/main/master.m) script. This is the primary script that individual researchers can use to enter paths to and specifications for their datasets. Examples of text files for subject numbers and demographic/clinical information are contained [here](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/blob/main/sub_ids_all_col2) and [here](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/blob/main/clinical_data_all_col2).  

2. Functions for four main steps for network analysis along with helper functions are contained in the [network_functions/](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/tree/main/network_functions) directory. The steps for network analysis should be executed in the following order:
    1. [find_common_rois_250_function](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/blob/main/network_functions/find_common_ROIs_250_function.m)
    2. [stack_sub_matrices_function_combat_enigma](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/blob/main/network_functions/stack_sub_matrices_function_combat_enigma.m)
    3. [calculate_group_level_community_structure_function_surprise_250.m](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/blob/main/network_functions/calculate_group_level_community_structure_function_surprise_250.m)
    4. [calculate_pre_defined_network_within_subject.m](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/blob/main/network_functions/calculate_pre_defined_network_within_subject.m)

3.  The user's system MUST have the Partitioning Cost Optimization (PACO) from Carlo Nicolini installed. Please see [the PACO GitHub](https://github.com/CarloNicolini/paco) for instructions and downloads. 

4.  The vast majority of this pipeline runs through Matlab, so familiarity with Matlab is essential. This pipeline also requires the [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/) in Matlab. 


## Documents and Reports

Before beginning analysis, please review the following articles which are ESSENTIAL to understanding the execution of this pipeline.  

-   ***[Ross et al. 2021 Resolution limit-free community detection reveals unique patterns of resting-state network connectivity in posttraumatic stress disorder: A PGC-ENIGMA PTSD Consortium investigation ](https://www.medrxiv.org/content/10.1101/2021.06.24.21259102v1)*** *This is the pre-print that describes the original investigation with the ENGIMA PTSD resting-state set, led by Marisa C. Ross, PhD*
-   ***[Nicolini, Bordier, & Bifone 2017 Community detection in weighted brain connectivity networks beyond the resolution limit:](https://www.sciencedirect.com/science/article/pii/S1053811916306449)*** *This is the primary resource for description of PACO and its use in functional brain networks*
-   ***[Yu et al. 2018 Statistical harmonization corrects site effects in functional connectivity measu rements from multi-site fMRI data](https://onlinelibrary.wiley.com/doi/epdf/10.1002/hbm.24241)*** *Key paper on ComBat harmonization for multi-site projects*

## Methodology

1) Preprocessed .nii files from HalfPipe resting-state pipeline were first [resmoothed](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/blob/main/resmooth_and_rescale) with an 8mm FWHM kernel for quality control of spatial normalization

2) Resting State Network Analysis with ComBat and Asymptotical Surprise Overview
   - Troubleshooting may be necessary to determine the best atlas for group-level analyses. You want to maximize spatial coverage amongst all participants but            also have ROIs small enough to be meaningful. The atlas used in the preprint was the 250-ROI parcellation from [Craddock et al. 2012](https://pubmed.ncbi.nlm.nih.gov/21769991/)
   - Use the master script to call all of the functions required to do network detection using Asymptotical Surprise for resting-state scans of all subjects. This        master script uses a 250-ROI atlas and 10000 PACO iterations to define group-level network structure, then applies the group-level network structure as an          atlas to individuals to calculate graph theory metrics of interest
   - The steps for network analysis should be executed in the following order:
    1. [find_common_rois_250_function](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/blob/main/network_functions/find_common_ROIs_250_function.m) 
        - Finds ROIs that are common across all subjects
        - Do trouble-shooting on spatial coverage and atlas size at this step
    2. [stack_sub_matrices_function_combat_enigma](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/blob/main/network_functions/stack_sub_matrices_function_combat_enigma.m)
        - stack connectivity matrices together and generate group-level correlation matrix 
        - use [ComBat](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/tree/main/ComBat_Functions) to harmonize across scanners
    3. [calculate_max_connected_threshold](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/blob/main/network_functions/calculate_max_connected_threshold.m)
        - calculate the correct percolation threshold to generate a large, sparse graph where no nodes are unconnected from group-level correlation matrix for PACO           to use
    4. [calculate_group_level_community_structure_function_surprise_250.m](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/blob/main/network_functions/calculate_group_level_community_structure_function_surprise_250.m)
        -  Use AS to define a group-level spatial map based on the group correlation matrix 
        -  This script calls PACO to do network detection. 
        -  Saves group-level correlation matrix, maximum value of surprise, number of unique modules 
        -  Want to MAXIMIZE surprise (higher S == better partitioning). Use this step to do more trouble-shooting/jackknife analysis on included subjects to                    maximize surprise
    5. [calculate_pre_defined_network_within_subject.m](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/blob/main/network_functions/calculate_pre_defined_network_within_subject.m)
        -  Use the group mask (network structure) to calculate network measures for each subject. 
        -  Saves individual spatial maps and graph theory indices (positive participation cofficient, clustering coefficient, etc.) for every subject.
        -  User can make changes to this function if other graph theory measures are desired. Anything in the Brain Connectivity Toolbox should be compatible 
3) Define network-specific connectivity metrics
    - Now that a group-level spatial map has been generated, use it to define networks of interest. Select the modules you want to investigate (the original               analysis selected the DMN and CEN, based on similarity to canonical topology and the boundaries defined in the [atlas from Yeo et al. 2011]                  (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3174820/)). 
    - [within_network_connectivity_and_pos_par.m](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/blob/main/within_network_connectivity_and_pos_par.m) provides some examples of how one might go         about calculating network-specific graph theory metrics (e.g. participation coefficient of the DMN)
4) Within- and between- group statistics and node-level analysis
    - [rs_group_stats.m](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/blob/main/rs_group_stats.m) provides example Linear Mixed Effects Models and       plots for within- and between-group comparisons of graph theory measures.
    - This script also calls [group_level_node_analysis.m](https://github.com/marisacross/pgc_engima_ptsd_networks_analysis/blob/main/network_functions/group_level_node_analysis.m) to determine critical values for each graph theory index. These values can be used to investigate whether individual nodes are driving network- or whole-brain findings. 
    - Options allow for either regression or partial correlation

## Table of Contents

-   **[ComBat_Functions:](top_level_folder/)** *contains all functions necessary to use ComBat Harmonization*
-   **[network_functions:](top_level_folder/)** *contains necessary functions for network analysis*
-   **[master.m:](top_level_file.ext)** *example master script to call all community detection functions and save results*
-   **[rs_group_stats.m:](top_level_file.ext)** *example LMEs and node-level analysis*

## References

Key references include:
- Ross, M.C. and the PGC-ENIGMA PTSD Consortium. Resolution limit-free community detection reveals unique patterns of resting-state network connectivity in posttraumatic stress disorder: A PGC-ENIGMA PTSD Consortium investigation. *medRxiv*, https://doi.org/10.1101/2021.06.24.21259102
- Bordier, C., Nicolini, C., Forcellini, G. & Bifone, A. Disrupted modular organization of primary sensory brain areas in schizophrenia. *Neuroimage Clin 18*, 682–693 (2018)
- Craddock, R. C., James, G. A., Holtzheimer, P. E., Hu, X. P. & Mayberg, H. S. A whole brain fMRI atlas generated via spatially constrained spectral clustering. *Hum Brain Mapp 33*, 1914–1928 (2012)
- Yu, M. et al. Statistical harmonization corrects site effects in functional connectivity measurements from multi-site fMRI data. *Hum Brain Mapp 39*, 4213–4227 (2018)
- Nicolini, C., Bordier, C. & Bifone, A. Community detection in weighted brain connectivity networks beyond the resolution limit. *NeuroImage 146*, 28–39 (2017)
- Aldecoa, R. & Marín, I. Surprise maximization reveals the community structure of complex networks. *Sci Rep 3*, 1060 (2013)
