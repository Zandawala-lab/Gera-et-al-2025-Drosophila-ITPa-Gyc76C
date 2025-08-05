# Gera-et-al-2025-Drosophila-ITPa-Gyc76C

This is the R code to replicate connectomic analysis in Figure 13 (and Fig 13 S1 & S2) from: <br>
Gera et al. (2025) Anti-diuretic hormone ITP signals via a guanylate cyclase receptor to modulate systemic homeostasis in *Drosophila* <br>
https://elifesciences.org/reviewed-preprints/97043


## Getting Started
Download the repository. Instructions are also found in the README_FolderStructure document.
<br>
<br>
## Download Files Needed From Codex: https://codex.flywire.ai/api/download?dataset=fafb <br>
If files have **not** been downloaded, the code will prompt you to do so. <br>

Files needed from Codex download are: <br>
1. Synapse Table: fafb_v783_princeton_synapse_table.csv
2. Classification / Hierarchical Annotations: classification.csv
3. Neurotransmitter Type Predictions: neurotransmitters.csv <br>
    * *Note*: we **renamed** the filename to be **neurotransmitters** (not neurons)
4. Cell Types: consolidated_cell_types.csv

Place the files in the 'input' directory: ITP connectome/**input**/
<br>

## Files & Folders Provided in this Repository

**input**
1. version.csv - v783
2. ITP_v783.csv <br>
    * Supplementary Table 3: https://elifesciences.org/reviewed-preprints/97043
3. NSC_v783.csv
    *  Available from McKim et al. (2024) NSC connectome paper
    *  Provided in Supplementary Table 3 - Supporting Information section: https://elifesciences.org/reviewed-preprints/102684#d1e2716
4. brainmesh.obj

  * If you download the repo from github, this folder should already exist with the files above provided <br>
  * The code checks to see if input/output folders exist, and will make them if they do not. However, without the input files provided (or needed from codex              download), the rest of the script will not run.
<br>

**output** <br>
  * Empty directory to start
  * Files & figures will be saved here when Figure_13.R is run <br>
<br>

**scripts** 
  * Figure_13.R is the main script to replicate analyses and figures <br>
<br>

## Other Files Provided
1. ITP connectome.Rproj
    * Open this file first for project settings and to run Figure_13.R
2. README_FolderStructure.docx
    * Contains instructions and images of folder structure for reference

<br>


## Questions:
tmckim[at]unr.edu




