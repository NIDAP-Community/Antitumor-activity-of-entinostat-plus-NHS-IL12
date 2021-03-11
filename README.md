# Antitumor-activity-of-entinostat-plus-NHS-IL12
Code to reproduce figures and tables in paper

This code was written in the NIH Integrated Data Analysis Portal (NIDAP), a user interface developed on the Foundry Platform (Palantir Technologies).

Code Authors: Margaret Cam, Thomas Joshua Meyer, Christian Sidak, Matthew Angel, Richard Finney, Jing Bian.

## Basic Workflow

Quality_Control.R > Color_by_Genes.R > ModScore.R > DEG_l2p.R

Note: unused R objects are deleted after each section/step to reduce memory load.

## Quality Control

Step 1: Filter & Generate Quality Control Plots from h5 files

Step 2: Generate Post-filter Data

Step 3: PCA & Normalization

Step 4: Combine and Renormalize

## Color by Genes

Produces figures 3b, 5a, 5e, 5g, 7e, 7m, S5b, S5c

## ModScore

Step 5: ModScore and Cell Classification (produces figure 5b)
Note: This code replaces gene expression information of irrelevant genes (e.g. Vamp4 and Vash2) with information from negative transcriptional markers (e.g. Cd4_neg, Sell_low). See Supplementary Table 3: Marker genes used for cell type identity by scRNAseq.

Step 6: Dotplot and Contingency Table (produces figures 3a, S3d, table S4)

## DEG_l2p

Step 7: DEG Analysis

Step 8: List to Pathway Visualisation (produces figures 4a, 5f, 6b, 6d, 7d, S7k, and tables S5, S6, S8, S9, S10, S11)

Note: List to Pathway Visualisation (l2p) requires l2p_0.0-1 (see attached tar.gz file). Github link: https://github.com/CCBR/l2p
