# Gera-et-al-2025-Drosophila-ITPa-Gyc76C

This is the R code to replicate connectomic & scRNA analyses in Figure 13 (and Supplemental Figures 1-3) from: <br>
Gera et al. (2025) Anti-diuretic hormone ITP signals via a guanylate cyclase receptor to modulate systemic homeostasis in *Drosophila* <br>
https://doi.org/10.7554/eLife.97043


## Getting Started
Download this repository<br>
<br>

There are 2 folders, one for each analysis type:
1. Connectome analyses -- ITP connectome folder
2. scRNA analyses -- ITP scRNA analysis folder<br>
<br>


Running the R scripts:<br>
1. Check that you have downloaded any files needed
2. Navigate to the folder of interest (ITP connectome or ITP scRNA analysis)
3. Open the ITP connectome.Rproj (or ITP scRNA analysis.Rproj)
     * Opens new Rstudio window with associated project settings (paths, etc.)
     * Located in the main folder
4. Open the Figure 13.R script (or Figure 13 supplement 3.R) and run code
     * Found in the **scripts** subfolder<br>
<br>


Instructions are also found in the README_FolderStructure documents in each folder<br>

   
<br>

## 1. Connectome Analysis (Figure 13 - Panels B - D & Supplement 1) <br>

### 1.1. Download Files Needed From Codex: https://codex.flywire.ai/api/download?dataset=fafb <br>
If files have **not** been downloaded, the code will prompt you to do so. <br>

Files needed from Codex download are: <br>
1. Synapse Table: fafb_v783_princeton_synapse_table.csv
2. Classification / Hierarchical Annotations: classification.csv
3. Neurotransmitter Type Predictions: neurotransmitters.csv <br>
    * *Note*: we **renamed** the filename to be **neurotransmitters** (not neurons)
4. Cell Types: consolidated_cell_types.csv

Place the files in the 'input' directory: ITP connectome/**input**/
<br>
<br>

### 1.2. Files & Folders Provided in: ITP connectome

**input**
1. version.csv - v783
2. ITP_v783.csv <br>
    * Supplementary Table 3: https://doi.org/10.7554/eLife.97043
3. NSC_v783.csv
    *  Available from McKim et al. (2024) NSC connectome paper
    *  Provided in Supplementary Table 3 - Supporting Information section: https://doi.org/10.7554/eLife.102684
4. brainmesh.obj<br>
<br>

* If you download the repo from github, this folder should already exist with the files above provided <br>
* The code checks to see if input/output folders exist, and will make them if they do not. However, without the input files provided (or needed from codex download), the rest of the script will not run.
<br>

**output** <br>
  * Empty directory to start
  * Files & figures will be saved here when Figure_13.R is run <br>
<br>

**scripts** 
  * Figure_13.R is the main script to replicate analyses and figures in the paper and supplement <br>
<br>

### 1.3. Other Files Provided
1. ITP connectome.Rproj
    * Open this file first for project settings and to run Figure_13.R
2. README_FolderStructure.docx
    * Contains instructions and images of folder structure for reference

<br>

## 2. scRNA Analysis (Figure 13 - Panels G-J & Supplement 3) <br>

### 2.1. Optional: Download Files Needed From SCope For the Brain & VNC
Brain: https://scope.aertslab.org/#/Davie_et_al_Cell_2018/Davie_et_al_Cell_2018%2FAerts_Fly_AdultBrain_Filtered_57k.loom/gene <br>
VNC: https://scope.aertslab.org/#/Davie_et_al_Cell_2018/Davie_et_al_Cell_2018%2FGoodwin_Fly_AdultVNC_elife54074.loom/gene <br>

Place the files in the 'input' directory: ITP scRNA analysis/**input**/
<br>

### Note: Running the script to reproduce the analysis does not require that you download these large files
* The subset of data needed is provided in the **input** directory to replicate analyses and figures


### 2.2. Files & Folders Provided in: ITP scRNA analysis

**input**
1. Neuropeptides.rda
2. Neuropeptide_receptors.rda
3. Monoamine_receptors.rda
4. Neurotransmitter_receptors.rda
5. brainITPcells.rda
6. VNCITPcells.rda

  * If you download the repo from github, this folder should already exist with the files above provided <br>

<br>

**output** <br>
  * Empty directory to start
  * PDF figures will be saved here after running Figure 13.R or Figure 13 supplement 3.R <br>
<br>

**scripts** 
  * Figure 13.R is the main script to replicate analyses and figures in the paper <br>
  * Figure 13 supplement 3.R is the script to replicate analyses and figures in the supplement
<br>

### 2.3. Other Files Provided
1. ITP scRNA analysis.Rproj
    * Open this file first for project settings and to run Figure 13.R and Figure 13 supplement 3.R
2. README_FolderStructure.docx
    * Contains instructions and images of folder structure for reference

<br>


## Questions:
mzandawala[at]unr.edu <br>
tmckim[at]unr.edu




