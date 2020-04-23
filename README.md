# cMAp_Terzi


# This is data from connectivety map
## Part1: RNA seq
It was performed with R on Data from EGF-Tgfa Nx mice.
The significant changes are stored in *tgfa_fixed2_parent_ShamvsNx_res_mouse.csv*

## Part2: cMAP
Data where downloaded from **GSE92742**, and analyzed with python/cmapPy. file is a file that can be modified by Jupyter notebook (Anaconda) or Google Colab

## Part3: Results
Saved in PowerPoint presentation format


## To compare with cMAP site
cMAP site is able to test list of upregulated and downregulated genes against their database. The results are **very different** that what the analysis shows.
A list of scores are saved in "clue_cMAP_HA1E_vs_Terzi.txt"



## Files:

**cMAP_code_corrected.R** : code for reading RNA-Seq counts, analyzing them and changing gene names into human analoge.

**cMAP_Terzi_Lisa.ipynb** : In python, code to analyse L1000 data and RNA-Seq. Convert them into (-1,0,1) and get correlation between the RNA-Seq and the change induced by the drugs.

**cMAP_Terzi_Lisa_noDMSO_normalization.ipynb** : In python, same code apart from removing normalization per plate.

**tgfa_fixed2_parent_ShamvsNx_res_mouse.csv** : File of differentially expressed genes between Tgfa p/mut. It is done twice, once for sham Tgfa-p/mut and once for Nx Tgfa-p/mut. The results are two merged DESeq2 tables, genes were traformed to analog human names. Includes: mean expression, log2 fold change, logFC SE, Wald-test score,  p-value and p-adjusted value. *output* of R code and *input* of python.

**cMAP_stats.pptx** : PowerPoint of the results.

**clue_cMAP_HA1E_vs_Terzi.txt** Table from the connectivity map site, reveals score of drugs vs top 150+150 genes in RNA-Seq experiment.
