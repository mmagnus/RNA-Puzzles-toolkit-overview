# RNA-Puzzles Toolkit Overview

## Overview of the RNA-Puzzle Toolkit in Jupyter:

See you see examples of usages of the Toolkit with Jupyter: 

[See The RNA-Puzzles Tookit Overview](https://github.com/mmagnus/RNA-Puzzles-toolkit-overview/blob/master/rna-puzzle-toolkit.ipynb)

![fig](figures/fig1.png)

**Figure 1.  Scheme of the RNA-Puzzles Toolkit.**
The Toolkit implements various tools for processing RNA structure files, and it is composed of three parts: tools for formatting and analyzing structural models. (A) The user can start an analysis with one or multiple models. (B) The first step is to standardize the formatting of analyzed structural models. (C) When the structures are standardized, the user can run an analysis for a single model, such as Clash Score evaluation or base pair identification using various methods. For multiple models, various comparison methods are implemented. The tools can be accessed via the command-line interface or Jupyter. The tools are provided as a Docker image that can be easily used. (D) Using the Toolkit, a normalized dataset of all submissions was generated and is available for download.

See also https://github.com/RNA-Puzzles/RNA-Puzzles-toolkit-overview

## RNA structure formatting tools:

  - [RNA-Puzzles Normalization](https://github.com/RNA-Puzzles/BasicAssessMetrics)
  - [rna-tools](https://github.com/RNA-Puzzles/rna-tools)

## RNA structure comparison tools:

  - [RMSD](https://github.com/RNA-Puzzles/BasicAssessMetrics)
  - [Deformation Profile](https://github.com/RNA-Puzzles/DeformationProfile)
  - [Interaction Network Fidelity](https://github.com/RNA-Puzzles/BasicAssessMetrics)
  - [Mean of Circular Quantities](https://github.com/tzok/mcq4structures)
  - [P value](https://github.com/RNA-Puzzles/BasicAssessMetrics)
  - [LCS-TA](https://github.com/tzok/mcq4structures)
  - [rnaqua](https://github.com/mantczak/rnaqua)
  - [Clash score](http://molprobity.biochem.duke.edu/)

## Datasets from RNA-Puzzles:

  - [Raw submission](https://github.com/RNA-Puzzles/RNA-Puzzles_dataset/tree/master/raw)
  - [Processed For Assessment](https://github.com/RNA-Puzzles/RNA-Puzzles_dataset/tree/master/normalized)
  - [Normalized submissions](https://github.com/mmagnus/RNA-Puzzles-Normalized-submissions)
  
## Docker Image for the RNA-Puzzles Toolkit:

- Docker Image to run the Toolkit (https://github.com/RNA-Puzzles/toolkit-ansible-role)

<!--
- rna-tools (https://github.com/RNA-Puzzles/rna-tools)
- RNA_assessment (https://github.com/RNA-Puzzles/RNA_assessment)
- mcq (https://github.com/RNA-Puzzles/mcq)
- rnaqua (https://github.com/RNA-Puzzles/rnaqua)
-->
