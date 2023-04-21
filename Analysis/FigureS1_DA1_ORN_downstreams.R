# Taisz, Galili et al., (2023)
# Figure S1F, H, I
# find all downstream cell types of DA1 (Or67d) ORNs that are not ALLNs in the hemibrain dataset
# ORN cell type and hemisphere annotation from Schlegel, Bates et al., (2021)
# see https://github.com/natverse/hemibrainr on how to access hemibrain data
# ----------------------------------------------------------------------------------------------
library(natverse)
library(ggplot2)
library(dplyr)

# set working directory to repository
setwd("/Users/GitHub/2023_cVA_Taisz_Galili")

# helper function around neuprintr::neuprint_connection_table to aggregate connectivity by cell type
# arguments:
# bodyids: a hemibrain cell type or bodyids (usually members of a cell type, but can be a subset, or a combination of cell types)
# prepost: whether to look upstream (PRE) or downstream (POST) of the query neurons
# returns a data frame with the names of connected cell type and their synaptic weights with the query. In the case upstream partners, also relative inputs are calculated.

neuprint_conn_by_type = function (bodyids, prepost, ...) {
  if (is.character(bodyids)) {
    bodyids = neuprintr::neuprint_search(bodyids, field = "type")$bodyid
  }
  conn_table = neuprintr::neuprint_connection_table(bodyids, prepost = prepost, all_segments = TRUE)
  meta = neuprintr::neuprint_get_meta(conn_table$partner)
  conn_table$type = meta$type
  conn_table_by_type = stats::aggregate(weight ~ type, data = conn_table, FUN = sum) # drops NAs
  if (prepost == "PRE") {
    conn_table_by_type = dplyr::arrange(conn_table_by_type, dplyr::desc(weight))
    conn_table_by_type
  } else {
    conn_table$post = meta$post
    # calculate all inputs (make sure that you only count them once, not every time a connection between your quary and the partner cell type is made)
    conn_table_sub = conn_table[, c(2, 5, 6)]
    all_input_df = conn_table_sub[!duplicated(conn_table_sub),]
    all_inputs = stats::aggregate(post ~ type, data = all_input_df, FUN = sum) # drops NAs
    all_inputs = all_inputs[match(conn_table_by_type$type, all_inputs$type),]
    
    if (all.equal(all_inputs$type, conn_table_by_type$type)) {
      conn_table_by_type$all_inputs = all_inputs$post
      conn_table_by_type$percent_input = (conn_table_by_type$weight / conn_table_by_type$all_inputs)*100
      conn_table_by_type = dplyr::arrange(conn_table_by_type, dplyr::desc(weight))
      conn_table_by_type
    } else {
      stop("Couldn't count relative inputs")
    }
  }
}

# load ORN ids with hemisphere information

ORN_df = read.csv("Data/DA1_ORN_side.csv")
ORN_R = ORN_df$bodyid[ORN_df$side == "R"]
ORN_L = ORN_df$bodyid[ORN_df$side == "L"]

ORN = ORN_df$bodyid
ORN_ds = neuprint_conn_by_type(ORN, prepost = "POST")

# there are a few ORNs where left-right information isn't clear
# these will be included when counting connectivity from all ORNs
# but not included when comparing ipsilateral vs. contralateral ORN input

ORN_NOside_ds = neuprint_conn_by_type(ORN[!(ORN %in% c(ORN_R, ORN_L))], prepost = "POST")
ORN_NOside_ds_PN = ORN_NOside_ds[(!grepl("LN", ORN_NOside_ds$type)) & (!grepl("ORN", ORN_NOside_ds$type)),]
ORN_NOside_us = neuprint_conn_by_type(ORN[!(ORN %in% c(ORN_R, ORN_L))], prepost = "PRE")
ORN_NOside_us_PN = ORN_NOside_us[(!grepl("LN", ORN_NOside_us$type)) & (!grepl("ORN", ORN_NOside_us$type)),]

# initially, anything that doesn't have LN (local neuron) in its hemibrain cell type
# will be referred to as PN (projection neuron)


# we want to identify every non-local neuron type (PN) that gets DA1 ORN ( = cVA) input
# we remove any cell type with very low number of absolute inputs (10 or fewer synapses from the 142 DA1 ORNs overall)
ORN_ds_PN = ORN_ds[(!grepl("LN", ORN_ds$type)) & (!grepl("ORN", ORN_ds$type)) & ORN_ds$weight > 10,]

# we then calculate the proportion of DA1 ORN inputs for each of these cell types to see how selective they are to cVA

# find the inputs of these PN cell types
DA1_PN_conn_l = lapply(as.list(ORN_ds_PN$type), neuprint_connection_table, prepost = "PRE")
# make empty data frame
DA1_input_ratio = as.data.frame(matrix(NA, nrow = length(DA1_PN_conn_l), ncol =7))
colnames(DA1_input_ratio) = c("PN_type", "ORN_weight", "DA1_ORN_ratio", "other_ORN_w", "DA1_PN_ratio", "other_PN_w", "nonORN_inputs")
DA1_input_ratio[, 1] = ORN_ds_PN$type
DA1_input_ratio[, 2] = ORN_ds_PN$weight
# fill data frame with values
for (i in 1:length(DA1_PN_conn_l)) {
  x = DA1_PN_conn_l[[i]]
  x$type = neuprint_get_neuron_names(x$partner)
  x_full_ORN_in = sum(x$weight[grep("ORN", x$type)])
  x_DA1_ORN_in = sum(x$weight[grep("ORN_DA1", x$type)])
  
  x_full_PN_in = sum(x$weight[grep("PN", x$type)])
  x_DA1_PN_in = sum(x$weight[grep("DA1_lPN", x$type)])
  
  x_nonORN_inputs = sum(x$weight) - x_full_ORN_in
  
  DA1_input_ratio[i, 3] = x_DA1_ORN_in/x_full_ORN_in
  DA1_input_ratio[i, 5] = x_DA1_PN_in/x_full_PN_in
  
  DA1_input_ratio[i, 4] = x_full_ORN_in - x_DA1_ORN_in
  DA1_input_ratio[i, 6] = x_full_PN_in - x_DA1_PN_in
  
  DA1_input_ratio[i, 7] = x_nonORN_inputs
}



# in columns "DA1_ORN_ratio" we calculated the proportion of DA1 ORN input from all ORN input (rather than all input, including LNs and other PNs)
# let's find PNs that have at least 5% of their ORN inputs from DA1 ORNs
DA1_input_ratio_thr = DA1_input_ratio[DA1_input_ratio$DA1_ORN_ratio > 0.05, ]


# we found all the cell types that get at least 5% of their ORn input from DA1 ORNs (and this 5% is more than 10 synapses)
# we will plot the absolute and relative input connectivity of these 13 cell types, and their morphology
DA1_input_ratio_thr$PN_type

# ----------------------------------------------------------------------------------------------
# tidy data frame for plotting

# calculate input synapse numbers / neuron for cell types
DA1_input_ratio_to_plot = as.data.frame(cbind(c(DA1_input_ratio_thr$PN_type, DA1_input_ratio_thr$PN_type, DA1_input_ratio_thr$PN_type),
                                              c(DA1_input_ratio_thr$ORN_weight, DA1_input_ratio_thr$other_ORN_w, DA1_input_ratio_thr$nonORN_inputs)))

DA1_input_ratio_to_plot$source = c(rep("Or67d ORNs", nrow(DA1_input_ratio_to_plot)/3),
                                   rep("other ORNs", nrow(DA1_input_ratio_to_plot)/3),
                                   rep("non ORN input", nrow(DA1_input_ratio_to_plot)/3))

colnames(DA1_input_ratio_to_plot) = c("cell_type", "input", "source")
DA1_input_ratio_to_plot$input = as.numeric(as.character(DA1_input_ratio_to_plot$input))


# count number of cells per cell type and normalise connections by that
ds_list = as.list(as.character(unique(DA1_input_ratio_to_plot$cell_type)))
ds_nl = lapply(ds_list, neuprint_read_neurons) # loading all PN morphologies might take some time
n_cells = unlist(lapply(ds_nl[1:13], length))
# for cell type OA-VMUa5 there is only one member of the cell type connected, but there are two (bilateral) neurons
n_cells[6] = 2
DA1_input_ratio_to_plot$input_perneuron = DA1_input_ratio_to_plot$input / rep(n_cells, 3)
DA1_input_ratio_to_plot$cell_type = factor(DA1_input_ratio_to_plot$cell_type, levels = as.character(unique(DA1_input_ratio_to_plot$cell_type)))
DA1_input_ratio_to_plot$source = factor(DA1_input_ratio_to_plot$source, levels = c("non ORN input", "other ORNs", "Or67d ORNs"))

all_input = DA1_input_ratio_to_plot %>%
        group_by(cell_type) %>%
        summarise(sum(input))

DA1_input_ratio_to_plot$rel_input =  DA1_input_ratio_to_plot$input / rep(all_input$`sum(input)`, 3)


# ----------------------------------------------------------------------------------------------
# input connectivity plots
# all cell types from the hemibrain dataset that:
# - don't contain the term 'LN'
# - have more than 10 inputs from DA1 ORNs
# - more than 5% of their ORN inputs is from DA1 ORNs

barwidth = 0.8

# Figure S1F
PN_DA1_ratio_abs = ggplot() + 
  geom_bar(data =  DA1_input_ratio_to_plot, 
           mapping = aes(x = cell_type, y = input_perneuron, fill = source), 
           stat="identity", 
           position='stack', 
           width = barwidth) +
  theme_classic() +
  scale_fill_manual(labels = c("non ORN input", "other ORNs", "Or67d ORNs"), values = c("grey", "#ffe404", "#865B30")) +
  guides(fill=guide_legend(title=NULL)) +
  labs(y = "# of inputs / neuron") +
  labs(x = "") +
  theme(legend.text=element_text(size=15),
        legend.position = c("none")) +
  theme(axis.title = element_text(size = 15)) + 
  theme(axis.text = element_text(size = 15)) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) # + facet_zoom(ylim = c(0, 200), show.area = T)
PN_DA1_ratio_abs


# Figure S1H
PN_DA1_ratio_rel = ggplot() + 
  geom_bar(data =  DA1_input_ratio_to_plot, 
           mapping = aes(x = cell_type, y = rel_input*100, fill = source), 
           stat="identity", 
           position='stack', 
           width = barwidth) +
  theme_classic() +
  scale_fill_manual(labels = c("non ORN input", "other ORNs", "Or67d ORNs"), values = c("grey", "#ffe404", "#865B30")) +
  guides(fill=guide_legend(title=NULL)) +
  labs(y = "% inputs") +
  labs(x = "") +
  theme(legend.text=element_text(size=15),
        legend.position = c("none")) +
  theme(axis.title = element_text(size = 15)) + 
  theme(axis.text = element_text(size = 15)) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) # + facet_zoom(ylim = c(0, 200), show.area = T)
PN_DA1_ratio_rel


pdf("Figure_panels/Figure S1F.pdf", height = 6, width = 6)
print(PN_DA1_ratio_abs)
dev.off()
pdf("Figure_panels/Figure S1H.pdf", height = 6, width = 6)
print(PN_DA1_ratio_rel)
dev.off()




# ----------------------------------------------------------------------------------------------
# morphology 3D plots

# omit DA1 lPN and DA1 lvPN as they were shown in Figure 1B, C and Figure S1A, B
ds_list = as.list(as.character(unique(DA1_input_ratio_to_plot$cell_type)))
ds_list_sub = as.character(ds_list)[c(-1, -4)]
# ds_nl = lapply(ds_list, neuprint_read_neurons)
ds_nl_plot = ds_nl[c(-1, -4)]

# set plot orientation of the hemibrain volume
zoom.hb<-0.6446092
userMatrix.hb<-structure(c(0.998458445072174, -0.0149619737640023, 0.0534479767084122, 
                           0, -0.0538569055497646, -0.02841392531991, 0.998144328594208, 
                           0, -0.0134155601263046, -0.999484062194824, -0.0291759837418795, 
                           0, 0, 0, 0, 1), .Dim = c(4L, 4L))
windowRect.hb<-c(116L, 1125L, 1121L, 2131L)
file_out = "Figure_panels/Figure S1I/"


# cols = brewer.pal(length(ds_nl_plot), "Spectral")
cols = rep("#422616", length(ds_nl_plot))

linewidth = 5
# plot_order = c(2, 1, 4, 3, 5, 6, 7, 8, 9, 10, 11)
plot_order = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)

for (i in 1:length(ds_nl_plot)) {
  nopen3d(zoom = zoom.hb, userMatrix = userMatrix.hb, windowRect=windowRect.hb)
  plot3d(hemibrainr::hemibrain.surf, alpha = .05)
  plot3d(ds_nl_plot[[i]], col = cols[i], lwd = linewidth)
  snapshot3d(paste(file_out, ds_list_sub[i], ".png", sep = ""))
}
# close all rgl windows
while (rgl.cur() > 0) { rgl.close() }
