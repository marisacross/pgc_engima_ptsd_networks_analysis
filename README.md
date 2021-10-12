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


## Documents and Reports

Before beginning analysis, please review the following articles which are ESSENTIAL to understanding the execution of this pipeline.  

-   ***[Ross et al. 2021 Resolution limit-free community detection reveals unique patterns of resting-state network connectivity in posttraumatic stress disorder: A PGC-ENIGMA PTSD Consortium investigation ](https://www.medrxiv.org/content/10.1101/2021.06.24.21259102v1)*** *This is pre-print that describes the original investigation with the ENGIMA PTSD resting-state set, led by Marisa C. Ross, PhD*
-   ***[Nicolini, Bordier, & Bifone 2017 Community detection in weighted brain connectivity networks beyond the resolution limit:](https://www.sciencedirect.com/science/article/pii/S1053811916306449)*** *This is the primary resource for description of PACO and its use in functional brain networks *

## Methodology

Briefly describe the most essential parts of your methodology, making
sure to include the names of datasets that were used and a link or
reference to any external datasets. Give an in-depth description of any
code or scripts you need to run the project, including the order they
should be run in and any unusual dependencies, quirks, or features you
think others should know about to reproduce your project.

## Table of Contents

A linked table of contents for all the files on the repo (you *do not*
need to include the README and .gitignore files in this TOC)

Use the following format:

-   **[top\_level\_folder:](top_level_folder/)** *a description of this
    folder (if necessary)*
    -   **[second\_level\_folder:](second_level_folder/)**
        -   **[file.ext:](second_level_folder/file.ext)** *a description
            of this file (required)*
-   **[top\_level\_file.ext:](top_level_file.ext)** *a description of
    this file (required)*

## References

No one likes plagiarizers, cite your work!
