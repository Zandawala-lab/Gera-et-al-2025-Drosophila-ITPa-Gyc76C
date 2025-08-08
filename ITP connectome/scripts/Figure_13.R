# --- Figure 13 ----------------------------------------------------------------
# ------------------------------------------------------------------------------
# This is the code for:
# Figure 13: Panels B-D (ITP input/output)

# ------------------------------------------------------------------------------
# Files to download from Codex before running this script:
# ------------------------------------------------------------------------------
# All csv files were downloaded from Codex on 07/27/2025
# https://codex.flywire.ai/api/download?dataset=fafb

# If files have *not* been downloaded, the code below will prompt you to do so.

# Files needed from Codex download are:
# 1.) Synapse Table: fafb_v783_princeton_synapse_table.csv
# 2.) Classification / Hierarchical Annotations: classification.csv
# 3.) Neurotransmitter Type Predictions: neurotransmitters.csv
      # *Note*: rename the filename to be: neurotransmitters (not neurons)
# 4.) Cell Types: consolidated_cell_types.csv

# Put files in the 'input' directory: ITP connectome/input/

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Load packages
# ------------------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(fafbseg)  # https://github.com/natverse/fafbseg
library(natverse) # https://github.com/natverse/natverse

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------
options(scipen=999)

PATH_input = "./input/"
PATH_output = "./output/"

# Check if folders exist, if not, make them
if(!dir.exists(file.path(PATH_input))){
  dir.create(file.path(PATH_input))
}
if(!dir.exists(file.path(PATH_output))){
  dir.create(file.path(PATH_output))
}

# Uncomment to check input_files for what exists
# input_files = list.files(path = PATH_input, full.names = FALSE, recursive = FALSE)

# Read version file (provided) & get # to use in output filenames - v783
v = read_delim(paste0(PATH_input,"version.csv"),
               col_types  =  cols(version  =  col_character()),delim  =  ";")
v = v$version[1]

# Folder name for saving figs + check if exists
fig_folder_name = "Figure_13"
if(!dir.exists(file.path(PATH_output, fig_folder_name))){
  dir.create(file.path(PATH_output, fig_folder_name))
}

# ------------------------------------------------------------------------------
# Options for saving plots & csv files
# ------------------------------------------------------------------------------

# When running script for the first time, set both to TRUE:
write_plots = TRUE           # TRUE - save/replicate figure plots
                             # FALSE - plots not saved outside of R
write_csv = TRUE             # TRUE - save processed data associated w/figures
                             # FALSE - data not saved outside of R

# ------------------------------------------------------------------------------
# Set colors
# ------------------------------------------------------------------------------

# Synapse location plots
l_NSC_ITP = "#838383" 
LNv_5th = "#00FFFF"  
LNd_ITP = "#ffe200"

syn_col <- c("l_NSC_ITP" = l_NSC_ITP, "5th_LNv" = LNv_5th, "LNd_ITP" = LNd_ITP)

# Super class plots
super_class_colors <- c("endocrine" = "#5f50a1", # other
                        "undefined" = "white",
                        "optic" = "black", # 999999
                        "descending" = "#c3905f",  
                        "visual_centrifugal" = "#f46d43",
                        "central" = "#3384b8", 
                        "visual_projection" = "#c8d684", 
                        "sensory" = "#b73545", 
                        "ascending" = "#a9d5a3")
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# Files downloaded from Codex:--------------------------------------------------

# Synapses
if (!paste0("fafb_v783_princeton_synapse_table.csv") %in% list.files(PATH_input)) {
  stop("please go to https://codex.flywire.ai/api/download and download the 
       synapses table for the current version and save it in './input'.")
}
synapses = read_delim(paste0(PATH_input,"fafb_v783_princeton_synapse_table.csv"),
                      delim = ",",
                      escape_double = FALSE,
                      col_types = cols(pre_root_id_720575940 = col_character(),
                                      post_root_id_720575940 = col_character()),
                      trim_ws = TRUE)

# Classification
if (!paste0("classification.csv") %in% list.files(PATH_input)) {
  stop("please go to https://codex.flywire.ai/api/download and download the 
       classification file for the current version and save it in './input'.")
}
classification <- read_delim(paste0(PATH_input, "classification.csv"),
                             delim = ",",
                             escape_double = FALSE,
                             col_types = cols(root_id = col_character(), 
                                              flow = col_character()),
                             trim_ws = TRUE)

# Neurotransmitters
if (!("neurotransmitters.csv") %in% list.files(PATH_input)) {
  stop("please go to https://codex.flywire.ai/api/download and download the 
       neurotransmitters file for the current version and save it in './input'.")
}
neurotransmitters <- read_delim(paste0(PATH_input, "neurotransmitters.csv"),
                             delim = ",",
                             escape_double = FALSE,
                             col_types = cols(root_id = col_character()),
                             trim_ws = TRUE)
# Cell types
if (!("consolidated_cell_types.csv") %in% list.files(PATH_input)) {
  stop("please go to https://codex.flywire.ai/api/download and download the 
       cell types file for the current version and save it in './input'.")
}
cell_types <- read_delim(paste0(PATH_input, "consolidated_cell_types.csv"),
                         delim = ",",
                         escape_double = FALSE,
                         col_types = cols(root_id = col_character(),
                                          primary_type = col_character(),
                                          `additional_type(s)` = col_character()),
                         trim_ws = TRUE)

# Files provided:---------------------------------------------------------------

# ITP ids identified from Codex
# Provided in current paper (Gera et al.)
# Supplementary Table 3: https://doi.org/10.7554/eLife.97043
ITP_ids = read_delim(paste0(PATH_input,"ITP_v",v,".csv"),
                     delim = ",",
                     col_types = cols(ITP_id = col_character()))

# IDs & names from McKim et al. NSC connectome to link ITP & NSC connections
# Provided in Supplementary Table 3 - Supporting Information section:
# https://doi.org/10.7554/eLife.102684
NSC = read_delim(paste0(PATH_input,"NSC_v",v,".csv"),
                 delim = ",",
                 col_types = cols(NSC_id = col_character()))

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Data processing, filtering, joining, & organization, etc.
# ------------------------------------------------------------------------------

# Adjust col names in dfs: -----------------------------------------------------
synapses <- synapses %>%
  rename(pre_root_id = pre_root_id_720575940,
         post_root_id = post_root_id_720575940)

cell_types <- cell_types %>%
  rename(cell_type = primary_type,
         cell_type_additional = `additional_type(s)`)

# Filter synapses: -------------------------------------------------------------

# Filter out autapses (currently not shown in Codex)
synapses_no_autapses <- synapses %>%
  filter(pre_root_id != post_root_id)

# IDs in synapses only contain last 9 digits of full id, add first 9
# Needed for consistent joining based on ids
synapses_no_autapses <- synapses_no_autapses %>%
  mutate(pre_root_id = paste0("720575940", as.character(pre_root_id)),
         post_root_id  = paste0("720575940", as.character(post_root_id)))

# ITP ids are pre, so find outputs
synapses_output_ITP <- synapses_no_autapses %>%
  semi_join(ITP_ids, by = c("pre_root_id" = "ITP_id"))

# ITP ids are post, so find inputs
synapses_input_ITP <- synapses_no_autapses %>%
  semi_join(ITP_ids, by = c("post_root_id" = "ITP_id"))

# The individual *INPUT* & *OUTPUT* synapse dfs are processed below
# Remove the large synapses dfs to prevent memory/lag issues
rm(synapses)
rm(synapses_no_autapses)

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Figure 13 - Panel B (left) - Plot all *INPUT* synapses of ITP
# Note: This is all data - >=5 synapse threshold is not yet used here
# ------------------------------------------------------------------------------

# filename: Figure_13B_ITP_input_all_synapses_v783
panel_name = "B"
panel_type = "_ITP_input_all_synapses_"
filename_input_syn = paste0(fig_folder_name, panel_name, panel_type)

# Open a new 3D plot
open3d()
# Set high resolution (4K) for plot
par3d(windowRect = c(0, 0, 3840, 2160))  

# Iterate over unique ITP names to plot synapses
for (i in rev(unique(ITP_ids$ITP_name))) {
  # Get ids for current ITP_name
  ids <- ITP_ids %>% filter(ITP_name == i) %>% pull(ITP_id)
  # Filter synapses for current ids
  synapses_filtered <- synapses_input_ITP %>% filter(post_root_id %in% ids)
  # Plot pre-synaptic positions
  plot3d(synapses_filtered$post_x, synapses_filtered$post_y, synapses_filtered$post_z,
         col = syn_col[i], size = 0.5, type = "s", add = TRUE)
}

# Add surface model
brainmesh <- readOBJ(paste0(PATH_input,"brainmesh.obj"))
# Plot
plot3d(brainmesh, add = TRUE, alpha = 0.1, col = "grey")
# Adjust view
view3d(userMatrix = rotationMatrix(90 * pi / 90, 1, 0, 0), zoom = 0.5)  

if(write_plots){
  # Export as html
  p<-rglwidget(webgl=TRUE, width = 1920, height = 1080)
  htmltools::save_html(p, file.path(PATH_output, fig_folder_name, 
                                    paste0(filename_input_syn, "v", v, ".html")))
  
  # Export as png
  png_filename <- file.path(PATH_output, fig_folder_name, 
                            paste0(filename_input_syn, "v", v, ".png"))
  rgl.snapshot(png_filename)
}

# Close 3D plot
close3d()
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Figure 13- Panel B (right) - Plot all *OUTPUT* synapses of ITP
# Note: This is all data - >=5 synapse threshold is not yet used here
# ------------------------------------------------------------------------------

# filename: Figure_13B_ITP_output_all_synapses_v783
panel_name = "B"
panel_type = "_ITP_output_all_synapses_"
filename_output_syn = paste0(fig_folder_name, panel_name, panel_type)

# Open a new 3D plot
open3d()
# Set high resolution (4K) for plot
par3d(windowRect = c(0, 0, 3840, 2160)) 

# Iterate over unique ITP names to plot synapses
for (i in rev(unique(ITP_ids$ITP_name))) {
  # Get ids for current ITP_name
  ids <- ITP_ids %>% filter(ITP_name == i) %>% pull(ITP_id)
  # Filter synapses for current ids
  synapses_filtered <- synapses_output_ITP %>% filter(pre_root_id %in% ids)
  # Plot pre-synaptic positions
  plot3d(synapses_filtered$pre_x, synapses_filtered$pre_y, synapses_filtered$pre_z,
         col = syn_col[i], size = 0.5, type = "s", add = TRUE)
}

# Add surface model
brainmesh <- readOBJ(paste0(PATH_input,"brainmesh.obj"))
# Plot
plot3d(brainmesh, add = TRUE, alpha = 0.1, col = "grey")
# Adjust view
view3d(userMatrix = rotationMatrix(90 * pi / 90, 1, 0, 0), zoom = 0.5)  

if(write_plots){
  # Export as html
  p<-rglwidget(webgl=TRUE, width = 1920, height = 1080)
  htmltools::save_html(p, file.path(PATH_output, fig_folder_name, 
                                    paste0(filename_output_syn, "v", v, ".html")))
  
  # Export as png
  png_filename <- file.path(PATH_output, fig_folder_name, 
                            paste0(filename_output_syn, "v", v, ".png"))
  rgl.snapshot(png_filename)
}

# Close 3D plot
close3d()
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# *INPUT* - Connectivity data for figure panels
# Data processing, joining, & organization, etc. (>=5 synapse threshold applied)
# ------------------------------------------------------------------------------

# Summarize synapse data :------------------------------------------------------

# Summarize synapses w/neuropil - individual #'s for every row (every connection)
ITP_input_sum <- synapses_input_ITP %>%
  group_by(pre_root_id, post_root_id, neuropil) %>%
  summarise(n_synapses = n())

# Summarize synapses w/o neuropil - total #'s per *unique* id combo
ITP_input_total_sum <- synapses_input_ITP %>%
  group_by(pre_root_id, post_root_id) %>%
  summarise(total_synapses = n())

# Join these^
ITP_input_sum <- ITP_input_sum %>%
  left_join(ITP_input_total_sum, by = c("pre_root_id", "post_root_id"))

# Join  w/info of interest to ITP inputs:---------------------------------------

# Join classification with neurotransmitters:-----------------------------------
classification <- left_join(classification, neurotransmitters, by = "root_id")

# Join cell types to classification :-------------------------------------------
# Includes neurotransmitters from above join^
classification <- left_join(classification, cell_types, by = "root_id")
classification <- classification %>%
  relocate(cell_type, .after = sub_class) %>%
  relocate(cell_type_additional, .after = cell_type)

# Prep data for join below
classification_join <- classification
# These are *INPUTS* to ITP, so classification info goes with pre_root_id
# Add pre_ to all col names
names(classification_join) <- paste0("pre_", names(classification_join))

# Join & apply threshold of >= 5 synapses:-------------------------------
ITP_input <- left_join(ITP_input_sum[ITP_input_sum$total_synapses >= 5, ], 
                       classification_join, by = "pre_root_id")

# Code check: unique(ITP_input$pre_root_id) # 106 - matches Codex
# autapses were filtered out in our data above, and do not appear in Codex

# Join dfs and rename cols - add info for ITP to ITP connections
ITP_join <- ITP_ids
colnames(ITP_join) <- c("ITP_name_post", "post_root_id", "ITP_hemisphere_post")
ITP_input <- left_join(ITP_input, ITP_join, by = "post_root_id")
colnames(ITP_join) <- c("ITP_name_pre", "pre_root_id", "ITP_hemisphere_pre")
ITP_input <- left_join(ITP_input, ITP_join, by = "pre_root_id")
ITP_input <- ITP_input %>%
  relocate(ITP_name_pre, ITP_hemisphere_pre, ITP_name_post, ITP_hemisphere_post, 
           .after = total_synapses) 

# Save thresholded data used for input figures: ITP_input_synthresh5_v783.csv
if(write_csv){
  write.csv(ITP_input, paste0(PATH_output, "ITP_input_synthresh5_v", v, ".csv"))
}

# Data & figures after this for *INPUT* include the >=5 synapse threshold ------

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Figure Panel C - *INPUT* Proportions
# ------------------------------------------------------------------------------

# Prepare & organize data:------------------------------------------------------

# Summarize data by ITP types
ITP_input_total <- ITP_input %>%
  group_by(ITP_name_post) %>%
  summarise(n_synapses_total = sum(total_synapses, na.rm = TRUE),
            n_pre_partners_total = length(unique(pre_root_id)))

# Summarize groups by ITP types and super_class
ITP_input_sum_grouped <- ITP_input %>%
  group_by(pre_super_class, ITP_name_post) %>%
  summarise(n_synapses_sum = sum(total_synapses, na.rm = TRUE),
            avrg_synapses = mean(total_synapses, na.rm = TRUE),
            n_pre_partners = length(unique(pre_root_id)),
            n_post_partners = length(unique(post_root_id)))

# Join summarized data from above^ for calcs
ITP_input_sum_grouped <- left_join(ITP_input_sum_grouped, ITP_input_total, by = "ITP_name_post")
ITP_input_sum_grouped$perc_of_input <- ITP_input_sum_grouped$n_synapses_sum / ITP_input_sum_grouped$n_synapses_total

# Simplify to only cols needed for plot
ITP_input_prop <- ITP_input_sum_grouped %>%
  group_by(ITP_name_post, pre_super_class) %>%
  summarise(perc_of_input = sum(perc_of_input))

# Set order for x-axis
ITP_input_prop$ITP_name_post <- factor(ITP_input_prop$ITP_name_post,
                                       levels = c("5th_LNv", "LNd_ITP", "l_NSC_ITP"))

# Save data from figure
# filename: Figure_13C_ITP_input_synapses_prop.csv
panel_name = "C"
panel_type = "_ITP_input_synapses_prop_"
filename_input_prop = paste0(fig_folder_name, panel_name, panel_type)
if(write_csv){
  write.csv(ITP_input_prop, file.path(PATH_output, fig_folder_name, 
                                      paste0(filename_input_prop, "v", v, ".csv")), 
            row.names = FALSE)
}

# Plot *INPUT* synapse proportions:---------------------------------------------

# filename: Figure_13C_ITP_PropInputSynapses_v783.pdf
panel_name = "C"
panel_type = "_ITP_PropInputSynapses_"
filename_input_prop_syn = paste0(fig_folder_name, panel_name, panel_type)

# Address sub- and superscript characters in x-axis names
custom_labels <- c(expression("5th-LN"[v]), expression("LN"["d"]^"ITP"),"l-NSC^ITP")

p <- ggplot(ITP_input_prop, aes(x = ITP_name_post, y = perc_of_input, fill = pre_super_class)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  facet_grid(scales = "free_x", space = "free_x") +
  scale_fill_manual(values = super_class_colors,
                    guide = guide_legend(nrow = 1)) +
  scale_x_discrete(labels = parse(text = custom_labels)) +
  ylab("Proportion of input synapses") +
  xlab("") +
  theme(panel.background = element_rect(fill = NA, color = NA),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0),
        axis.text.y = element_text(angle = 0),
        axis.ticks.x = element_blank())
print(p)

# Save plot
if(write_plots){
  ggsave(file.path(PATH_output, fig_folder_name, paste0(filename_input_prop_syn, "v", v,".pdf")), 
         plot = p, width = 12, height = 10, units = "cm")
}

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Figure Panel D - *INPUT* - Numbers of neurons based on super class (schematic)
# ------------------------------------------------------------------------------

# filename: Figure_13D_ITP_input_by_superclass.csv
panel_name = "D"
panel_type = "_ITP_input_by_superclass_"
filename_input_superclass = paste0(fig_folder_name, panel_name, panel_type)

ITP_input_super_class <- ITP_input %>%
  group_by(ITP_name_post, pre_super_class) %>%
  summarise(n_pre_partners_total = length(unique(pre_root_id))) %>%
  arrange(ITP_name_post, desc(n_pre_partners_total))

# Save data for figure
if(write_csv){
  write.csv(ITP_input_super_class, file.path(PATH_output, fig_folder_name, 
                                             paste0(filename_input_superclass,"v", 
                                                    v, ".csv")), row.names = FALSE)
}

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure Panel D - *INPUT* IDs for neuroglancer neurons colored by super class
# ------------------------------------------------------------------------------

# filename_all: Figure_13D_synthresh5_all_input_to_5th_LNv_from_central_v783.csv
# filename_top: Figure_13D_top10_input_to_5th_LNv_from_central_v783.csv
panel_name = "D"
output_filename = "_top10_input_to_"
output_filename_all = "_synthresh5_all_input_to_"

# Folder names for saving figs & check if exists
folder_name_top = "Top_inputs_ITP"
if(!dir.exists(file.path(PATH_output, fig_folder_name, folder_name_top))){
  dir.create(file.path(PATH_output, fig_folder_name, folder_name_top))
}
# Folder name for saving figs & check if exists
folder_name_all = "ITP_all_inputs"
if(!dir.exists(file.path(PATH_output, fig_folder_name, folder_name_all))){
  dir.create(file.path(PATH_output, fig_folder_name, folder_name_all))
}

# Subset of cols to get in loop (don't need all orig cols)
col_to_keep = c("pre_root_id","pre_cell_type","pre_cell_type_additional", "pre_class",
               "pre_super_class", "ITP_name_pre","ITP_name_post","ITP_hemisphere_post",
               "post_root_id", "neuropil", "n_synapses", "total_synapses")

# Outer loop over 3 types of ITP_name_post: LNd_ITP, 5th_LNv, l_NSC_ITP
for (itp_type in unique(ITP_input$ITP_name_post)) {
  # Subset input data to get only rows that match current itp_type in ITP_name_post col & all cols
  ITP_input_subset <- ITP_input[ITP_input$ITP_name_post == itp_type, ]
  # Inner loop to get relevant input (pre) super classes based on data for current itp_type
  for (i in unique(ITP_input_subset$pre_super_class)) {
    if (!is.na(i)) {                                    
      # Filter data for current ITP ids
      ITP_subset_superclass <- unique(ITP_input_subset[ITP_input_subset$pre_super_class == i,][col_to_keep])
      # Name and save file
      if(write_csv){
        filename <- paste0(fig_folder_name, panel_name, output_filename_all, itp_type,
                           "_from_", i, sep = "")
        write.csv(ITP_subset_superclass, file.path(PATH_output, fig_folder_name,
                                                   folder_name_all, 
                                                   paste0(filename,"_v", v,".csv")))
      }
      # Assign values to pre_cell_type from pre_class or pre_root_id if they are NA
      ITP_subset_superclass <- ITP_subset_superclass %>%
        mutate(pre_cell_type = ifelse(
            is.na(pre_cell_type),
            ifelse(is.na(pre_class), pre_root_id, pre_class),
            pre_cell_type))
      
      # Top 10 pre_cell_types by # of synapses
      top_tmp = ITP_subset_superclass %>%
        group_by(pre_cell_type) %>%
        summarize(total_syn = sum(n_synapses))
      
      # Get the top10 names of pre_cell_types to filter the data based on
      top10_tmp = top_tmp[order(-top_tmp$total_syn),]$pre_cell_type[1:10]
      ITP_subset_superclass_top = ITP_subset_superclass[ITP_subset_superclass$pre_cell_type %in% top10_tmp,]
      
      # Create a df with pre_cell_type and their respective order
      order_df = data.frame(pre_cell_type = top10_tmp, order = 1:10)
      # Merge order information into ITP_subset_superclass_top
      ITP_subset_superclass_top = merge(ITP_subset_superclass_top, order_df, by = "pre_cell_type")
      
      # Get previous # of neurons from ITP_input_super_class data
      value_check <- ITP_input_super_class %>%
        filter(ITP_name_post == itp_type, pre_super_class == i) %>%
        pull(n_pre_partners_total) 
      # Check that # neurons here matches above
      if (length(unique(ITP_subset_superclass$pre_root_id)) != value_check) {
        stop(paste0(
          "Numbers do not match for ", itp_type, " ", i, 
          ". Go back to above code & check *ITP_input_super_class*. This is in Figure",
          " Panel D - *INPUT* - Numbers of neurons based on super class (schematic)."
        ))
      }
      # Name and save file
      if(write_csv){
        filename <- paste0(fig_folder_name, panel_name, output_filename, itp_type,
                           "_from_", i, sep = "")
        write.csv(ITP_subset_superclass_top,file.path(PATH_output, fig_folder_name,
                                                      folder_name_top, 
                                                      paste0(filename,"_v", v,".csv")))
      }
    }
  }
}

# pre_root_ids from each file based on ITP type & super class 
# (e.g.Figure_13D_top10_input_to_LNd_ITP_from_central_v783.csv) are put into
# neuroglancer (https://edit.flywire.ai/) for visualization & screenshots

#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# *OUTPUT* - Connectivity data for figure panels
# Data processing, joining, & organization, etc. (>=5 synapse threshold applied)
# ------------------------------------------------------------------------------

# Summarize synapse data :------------------------------------------------------

# Summarize synapses w/neuropil - individual #'s for every row (every connection)
ITP_output_sum <- synapses_output_ITP %>%
  group_by(pre_root_id, post_root_id, neuropil) %>%
  summarise(n_synapses = n())

# Summarize synapses w/o neuropil - total #'s per *unique* id combo
ITP_output_total_sum <- synapses_output_ITP %>%
  group_by(pre_root_id, post_root_id) %>%
  summarise(total_synapses = n())

# Join these^ dfs
ITP_output_sum <- ITP_output_sum %>%
  left_join(ITP_output_total_sum, by = c("pre_root_id", "post_root_id"))

# Join dfs w/info of interest to ITP outputs:-----------------------------------

# Assumes classification, neurotransmitters, and cell types already joined above

# Prep data for join below
classification_join <- classification
# These are *OUTPUTS* from ITP, so classification info goes with prost_root_id
# Add post_ to all col names
names(classification_join) <- paste0("post_", names(classification_join))

# Join & apply threshold of >= 5 synapses:-------------------------------
ITP_output <- left_join(ITP_output_sum[ITP_output_sum$total_synapses >= 5, ], 
                        classification_join, by = "post_root_id")

# Code check: unique(ITP_output$post_root_id) # 301 - matches Codex
# autapses were filtered out in our data above, and do not appear in Codex

# Join dfs and rename cols - add info for ITP to ITP connections
ITP_join <- ITP_ids
colnames(ITP_join) <- c("ITP_name_post", "post_root_id", "ITP_hemisphere_post")
ITP_output <- left_join(ITP_output, ITP_join, by = "post_root_id")
colnames(ITP_join) <- c("ITP_name_pre", "pre_root_id", "ITP_hemisphere_pre")
ITP_output <- left_join(ITP_output, ITP_join, by = "pre_root_id")

ITP_output <- ITP_output %>%
  relocate(ITP_name_pre, ITP_hemisphere_pre, ITP_name_post, ITP_hemisphere_post, 
           .after = total_synapses) 

# Save thresholded data used for input figures: ITP_output_synthresh5_v783.csv
if(write_csv){
  write.csv(ITP_output, paste0(PATH_output, "ITP_output_synthresh5_v", v, ".csv"))
}

# Data & figures after this for *OUTPUT* include the >=5 synapse threshold -----


#-------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure Panel D - *OUTPUT*  - Numbers of neurons based on super class (schematic)
# ------------------------------------------------------------------------------

# filename: Figure_13D_ITP_output_by_superclass_v783.csv
panel_name = "D"
panel_type = "_ITP_output_by_superclass_"
filename_output_superclass = paste0(fig_folder_name, panel_name, panel_type)

ITP_output_super_class <- ITP_output %>%
  group_by(ITP_name_pre, post_super_class) %>%
  summarise(n_post_partners_total = length(unique(post_root_id))) %>%
  arrange(ITP_name_pre, desc(n_post_partners_total))

# Save data for figure
if(write_csv){
  write.csv(ITP_output_super_class, file.path(PATH_output, fig_folder_name, 
                                              paste0(filename_output_superclass, "v", 
                                                     v, ".csv")), row.names = FALSE)
}

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure Panel D - *OUTPUT*  IDs for neuroglancer neurons colored by super class
# ------------------------------------------------------------------------------

# filename_all: Figure_13D_synthresh5_all_output_from_5th_LNv_to_ascending_v783.csv
# filename_top: Figure_13D_top10_output_from_5th_LNv_to_ascending_v783.csv
panel_name = "D"
output_filename_top = "_top10_output_from_"
output_filename_all = "_synthresh5_all_output_from_"

# Folder names for saving figs & check if exists
folder_name_top = "Top_outputs_ITP"
if(!dir.exists(file.path(PATH_output, fig_folder_name, folder_name_top))){
  dir.create(file.path(PATH_output, fig_folder_name, folder_name_top))
}
# Folder names for saving figs & check if exists
folder_name_all = "ITP_all_outputs"
if(!dir.exists(file.path(PATH_output, fig_folder_name, folder_name_all))){
  dir.create(file.path(PATH_output, fig_folder_name, folder_name_all))
}

# Subset of cols to get in loop (don't need all orig cols)
col_to_keep =c("pre_root_id","post_cell_type","post_cell_type_additional", "post_class",
               "post_super_class", "ITP_name_pre","ITP_name_post","ITP_hemisphere_post",
               "post_root_id", "neuropil", "n_synapses", "total_synapses")

# Outer loop over 3 types of ITP_name_pre: LNd_ITP, 5th_LNv, l_NSC_ITP
for (itp_type in unique(ITP_output$ITP_name_pre)) {
  # Subset output data to get only rows that match current itp_type in ITP_name_post col & all cols
  ITP_output_subset <- ITP_output[ITP_output$ITP_name_pre == itp_type, ]
  # Inner loop to get relevant output (post) super classes based on data for current itp_type
  for (i in unique(ITP_output_subset$post_super_class)) {
    if (!is.na(i)) {
      # Filter data for current ITP ids
      ITP_subset_superclass <- unique(ITP_output_subset[ITP_output_subset$post_super_class == i,][col_to_keep])
      # Name and save file
      if(write_csv){
        filename <- paste0(fig_folder_name, panel_name, output_filename_all, itp_type, 
                           "_to_", i, sep = "")
        write.csv(ITP_subset_superclass, file.path(PATH_output, fig_folder_name, 
                                                   folder_name_all, 
                                                   paste0(filename,"_v", v,".csv")))
      }
      # Assign values to post_cell_type from class or post_root_id if they are NA
      ITP_subset_superclass <- ITP_subset_superclass %>%
        mutate(post_cell_type = ifelse(
            is.na(post_cell_type),
            ifelse(is.na(post_class), post_root_id, post_class),
            post_cell_type))
      
      # Top 10 cell_types by # of synapses
      top_tmp = ITP_subset_superclass %>%
        group_by(post_cell_type) %>%
        summarize(total_syn = sum(n_synapses))
      
      # Get the top10 names of post_cell_types to filter the data based on
      top10_tmp = top_tmp[order(-top_tmp$total_syn),]$post_cell_type[1:10]
      ITP_subset_superclass_top = ITP_subset_superclass[ITP_subset_superclass$post_cell_type %in% top10_tmp,]
      
      # Create a df with post_cell_type and their respective order
      order_df = data.frame(post_cell_type = top10_tmp, order = 1:10)
      # Merge order information into ITP_subset_superclass_top
      ITP_subset_superclass_top = merge(ITP_subset_superclass_top, order_df, by = "post_cell_type")
      
      # Get previous # of neurons from ITP_output_super_class data
      value_check <- ITP_output_super_class %>%
        filter(ITP_name_pre == itp_type, post_super_class == i) %>%
        pull(n_post_partners_total) 
      # Check that # neurons here matches above
      if (length(unique(ITP_subset_superclass$post_root_id)) != value_check) {
        stop(paste0(
          "Numbers do not match for ", itp_type, " ", i, 
          ". Go back to above code & check *ITP_output_super_class*. This is in Figure",
          " Panel D - *OUTPUT* - Numbers of neurons based on super class (schematic)."
        ))
      }
      # Name and save file
      if(write_csv){
        filename <- paste0(fig_folder_name, panel_name, output_filename_top, itp_type,
                           "_to_", i, sep = "")
        write.csv(ITP_subset_superclass_top,file.path(PATH_output, fig_folder_name,
                                                      folder_name_top, 
                                                      paste0(filename,"_v", v,".csv")))
      }
    }
  }
}

# post_root_ids from each file based on ITP type & super class 
# (e.g.Figure_13D_top10_output_from_5th_LNv_to_central_v783.csv) are put into
# neuroglancer (https://edit.flywire.ai/) for visualization & screenshots

# ------------------------------------------------------------------------------
