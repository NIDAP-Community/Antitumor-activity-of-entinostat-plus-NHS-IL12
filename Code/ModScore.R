# Load relevant libraries, version numbers are in comments

library(quantmod) # 0.4.18
library(grid) # 4.0.4
library(data.table) # 1.14.0
library(Seurat) # 3.1.5
library(tidyverse) # 1.3.0
library(gridExtra) # 2.3
library(ggpmisc) # 0.3.8-1
library(Rfast) # 2.0.1
library(dplyr) # 1.0.5
library(Matrix) # 1.3-2
library(cowplot) # 1.1.1
library(RColorBrewer) # 1.1-2
library(ggrepel)

###---------------------------------------------------- ModScore and Cell Classification----------------------------------------------------------------------

# Load Inputs
# Seurat Object from Quality Control
SO <- SO_merge

SO_length <- length(SO@meta.data)

# Transcriptional markers used to identify cells (Note: This code replaces 
# gene expression information of irrelevant genes (e.g. Vamp4 and Vash2) with information from negative transcriptional markers (e.g. Cd4_neg, Sell_low)

# The csv used in this code was formatted to account for this)
Markers <- read.csv("") # load Gene_List (formatted)

# Returns the unique samples in the seurat object
samples = unique(SO@meta.data$sample_name) 

# Formatting
colnames(SO@meta.data) <- gsub("orig_ident","orig.ident",colnames(SO@meta.data))

if("active.ident" %in% slotNames(SO)){
  sample_name = as.factor(SO@meta.data$orig.ident)
  names(sample_name)=names(SO@active.ident)
  SO@active.ident <- as.factor(vector())
  SO@active.ident <- sample_name
  SO.sub = subset(SO, ident = samples)
} else {
  sample_name = as.factor(SO@meta.data$orig.ident)
  names(sample_name)=names(SO@active.ident)
  SO@active.ident <- as.factor(vector())
  SO@active.ident <- sample_name
  SO.sub = subset(SO, ident = samples)
} 

# Original markers
orig_markers <- c("Cd4","Foxp3","Cd8a","Ly6g","Nr4a1","Sell","Ly6c1","Adgre1","Mrc1","Cd38","Cd3e","Prf1")

# Random list of (irrelevant) genes to be overwritten by Cd4_neg, Foxp3_neg, Cd8a_neg, etc adjusted ModuleScore values
neg_markers_names <- c("Vamp4","Vash2","Vav2","Vapb","Vav3","Vamp3","Vamp5","Vamp8","Vax2","Vamp1","Vasp","Myb")

# Append neg_markers_names to rownames of SO.sub
neg_markers_list <- list() # Create a list for storage and retrieval

# Calculate adjusted gene expression for negative markers
for (i in seq_along(orig_markers)){
  
  SO.sub$SCT@scale.data[neg_markers_names[i],] <- max(SO.sub$SCT@scale.data[orig_markers[1],]) - SO.sub$SCT@scale.data[orig_markers[1],]  
  
}

marker.list <- as.list(Markers)

# Thresholds for bimodal distribution of Module Score density plot
thres_vec <- c(0.1,0.1,0.09,0.41,0.43,0.5,0.625,0.25,0.09,0.05,0.25,0.215)
names(thres_vec) <- colnames(Markers)

# Generate colored tsne plots
i = 0
for (i in seq_along(marker.list)) {
  print(names(marker.list[i]))
  present=lapply(marker.list[[i]], function(x) x %in% rownames(SO.sub@assays$SCT@data)) # apply function(x) x %in% rownames(SO.sub) to each element of marker.list
  absentgenes = unlist(marker.list[[i]])[present==FALSE]
  presentgenes = unlist(marker.list[[i]])[present==TRUE]
  print(paste0("Genes not present: ",paste0(absentgenes,collapse=",")))
  print(paste0("Genes present: ",paste0(presentgenes,collapse=",")))
  
  # Calculate Module Scores and append values to metadata
  if(length(presentgenes) > 0){
    SO.sub=AddModuleScore(SO.sub,marker.list[i],name = names(marker.list[i]))
    
    m = 0
    m = paste0(names(marker.list[i]),"1")
    SO.sub@meta.data[[m]] <- scales::rescale(SO.sub@meta.data[[m]], to=c(0,1))
        clusid = SO.sub@meta.data[[m]]
        d <- density(clusid)
    
        midpt.thres = thres_vec[i] 
    
        p1 <- DimPlot(SO.sub, reduction = "umap", group.by = "ident")
    
        clusmat=data.frame(umap1=p1$data$UMAP_1,umap2=p1$data$UMAP_2, clusid=as.numeric(SO.sub@meta.data[[m]]))
    
        clusmat %>% group_by(clusid) %>% summarise(umap1.mean=mean(umap1), umap2.mean=mean(umap2)) -> umap.pos
        title=as.character(m)
        clusmat %>% dplyr::arrange(clusid) -> clusmat
    
    # Dimension reduction plot
        g <- ggplot(clusmat, aes(x=umap1, y=umap2)) +
          theme_bw() +
          theme(legend.title=element_blank()) +
          geom_point(aes(colour=clusid),alpha=0.5,shape = 20,size=1) +
          scale_color_gradientn(colours = c("blue4","lightgrey", "red"), values = scales::rescale(c(0,midpt.thres/2,midpt.thres,(midpt.thres+1)/2,1), limits = c(0, 1))) + guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
          xlab("umap-1") + ylab("umap-2")
    
        clusid.df <- data.frame(id=SO.sub@meta.data$orig.ident,ModuleScore=SO.sub@meta.data[[m]])
        g1 = RidgePlot(SO.sub,features=m,group.by="orig.ident") + theme(legend.position="top", title = element_blank()) + geom_vline(xintercept = midpt.thres, linetype = "dashed", color = "red3") + scale_x_continuous(breaks = seq(0,1,0.1))
    
    # Violin Plot
        g2 = ggplot(clusid.df,aes(x=id,y=ModuleScore)) + geom_violin(aes(fill=id)) +  theme_classic() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title = element_blank(), panel.background = element_blank(), axis.text.x=element_blank(),legend.text=element_text(size=rel(0.8)),legend.position="top" ) + guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) + geom_hline(yintercept = midpt.thres, linetype = "dashed", color = "red3")
    
    # Color gradient density plot
        g3 = ggplot(data.frame(x = d$x, y = d$y), aes(x, y)) + xlab("ModuleScore") + ylab("Density") + geom_line() +
          geom_segment(aes(xend = d$x, yend = 0, colour = x)) + scale_y_log10() +
          scale_color_gradientn(colours = c("blue4","lightgrey", "red"), values = scales::rescale(c(0,midpt.thres/2,midpt.thres,(midpt.thres+1)/2,1), limits = c(0, 1))) + geom_vline(xintercept = midpt.thres, linetype = "dashed", color = "red3") + theme(legend.title = element_blank())
    
        last.figure = grid.arrange(g,g1,g2,g3, ncol=2, top=names(marker.list[i]))
  }
}

## Cell Type Classification
# Formatting metadata column names
names_repl <- substr(colnames(SO.sub@meta.data[(length(SO@meta.data) +1): length(SO.sub@meta.data)]), 1, nchar(colnames(SO.sub@meta.data[(length(SO@meta.data) +1): length(SO.sub@meta.data)]))-1)
names_orig <- colnames(SO.sub@meta.data[(length(SO@meta.data) +1): length(SO.sub@meta.data)])
setnames(SO.sub@meta.data, names_orig, names_repl) # remove the "1" from the end of the Module Score columns

### Classification of cells from Module Scores
## First classify cells in general class:
General_class <- c("CD8_T","CD4_T","Tregs","Macrophages","PMN_neutrophils","Monocytes","DCs","pDC","B_cells","NKs")

# Shorten metadata to contain only Module Score values of General Classes of Cells
trunc_vec <- c((length(SO@meta.data) +1): length(SO.sub@meta.data))
SO_Trunc_Metadata_General <- SO.sub@meta.data[trunc_vec][,General_class]

General_thres_vec <- c(0.1,0.1,0.09,0.41,0.43,0.5,0.625,0.25,0.09,0.05)
# Predict Cell Types based on Comparison of Module Scores and whether Module Scores are above threshold values 
Predict_Cell_from_ModScore <- function(ModScore_Metadata,thres_vec,output_name){
  max_col_vector <- max.col(ModScore_Metadata) # vector containing the column number with the maximum value
  
  # Returns a vector containing the highest ModScores of ModScore_Metadata
  hi_values <- c()
  i = 0
  SO_rows <- 1:nrow(ModScore_Metadata)
  for (i in seq_along(SO_rows)){
    hi_values[i] <- ModScore_Metadata[i,max_col_vector[i]]
  }
  
  # Returns a logical vector of whether the highest ModScores in ModScore_Metadata exceeds manual threshold 
  thres_logical_vector <- as.numeric(hi_values > thres_vec[max.col(ModScore_Metadata)])
  
  # Filter vectors
  filter_vector <- max_col_vector * thres_logical_vector
  filter_vector <- filter_vector + 1 # add 1 into filter_matrix to account for "unknown" cell type in appended_names
  
  # Original names appended to "unknown" classification for cells with ModScores below threshold
  appended_names <- c("unknown", names(ModScore_Metadata))
  
  # Added the names into a Likely_CellType Column
  dupl_data <- ModScore_Metadata
  dupl_data[,"Likely_CellType"] <- appended_names[filter_vector]
  assign(output_name,dupl_data,envir = globalenv())
}

Predict_Cell_from_ModScore(SO_Trunc_Metadata_General,General_thres_vec,"General_output")

### Classify cells in subclass (M1 and M2)
Sub_class <- c("Macrophages","M1","M2") # Must include Macrophages: if M1 and M2 scores don't exceed that of Macrophages, cells will remain Macrophages

# Subset out cells predicted to be Macrophages. This data will be used to further classify M1 and M2 Macrophages from General Macrophage population
Macrophage_cells <- rownames(General_output[General_output$Likely_CellType == "Macrophages",])
SO_Trunc_Metadata_M1M2 <- SO.sub@meta.data[trunc_vec][Macrophage_cells,][,Sub_class]

# Create general annotation with barcoded cell and corresponding identity. M1 and M2 identities will be appended to this table after second round of classification
SO_Trunc_Metadata_no_M1M2 <- General_output[!General_output$Likely_CellType == "Macrophages",]
non_Macrophage_cells <- rownames(SO_Trunc_Metadata_no_M1M2)
gen_annot_table <- data.frame(cells = non_Macrophage_cells, identity = SO_Trunc_Metadata_no_M1M2$Likely_CellType)

# Repeat Module Score Comparison and Cell Prediction with Macrophage Subset:
Macrophage_thres_vec <- c(0.41,0.25,0.215) # If M1 or M2 ModScore is higher than that of General Macrophages, than the cell will be classified as M1 or M2. Otherwise cell will remained with General Macrophage label
Predict_Cell_from_ModScore(SO_Trunc_Metadata_M1M2,Macrophage_thres_vec,"M1M2_output")

# Append updated Macrophage, M1, M2 metadata to gen_annot_table
Mac_M1M2_anot_table <- data.frame(cells = rownames(M1M2_output), identity = M1M2_output$Likely_CellType)
Final_annot_table <- rbind(gen_annot_table,Mac_M1M2_anot_table)

metadata <- SO.sub@meta.data
metadata$Likely_CellType <- Final_annot_table$identity[match(rownames(SO.sub@meta.data),Final_annot_table$cells)]

# Append final annotation table to Seurat Object
SO.sub@meta.data <- metadata

rm(list=setdiff(ls(), c("SO","SO.sub")))

#Generate Figure 5b

# Inputs.
so <- SO.sub
reduction = "umap"
Keep <- TRUE

## Replace dots in metadata column names with underscores.
colnames(so@meta.data) = gsub("\\.", "_", colnames(so@meta.data))

## Original color-picking code.
n <- 2e3
seed = 10
set.seed(seed)
ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
ourColorSpace <- as(ourColorSpace, "LAB")
distinctColorPalette <-function(k=1,seed) {
  currentColorSpace <- ourColorSpace@coords
  # Set iter.max to 20 to avoid convergence warnings.
  set.seed(seed)
  km <- kmeans(currentColorSpace, k, iter.max=20)
  colors <- unname(hex(LAB(km$centers)))
  return(colors)
}

Barcode <- rownames(so@meta.data)

## User-selected metadata column(s) is/are used to set idents.
metacols = c("Likely_CellType")
if(length(metacols) > 1) {
  Filter.orig = paste0(so@meta.data[[metacols[1]]],"_",so@meta.data[[metacols[2]]])
  colname <- paste0(metacols[1], "_", metacols[2])
  so <- AddMetaData(so, Filter.orig, col.name=colname)
  barcodes <- names(Idents(so))
  so@meta.data %>% dplyr::arrange(factor(Barcode), levels = barcodes) %>% pull(metacols[1]) -> metadat_sort
  Idents(so) <- factor(metadat_sort)
} else {
  Filter.orig = so@meta.data[[metacols[1]]]
  colname <- metacols[1]
  barcodes <- names(Idents(so))
  so@meta.data %>% dplyr::arrange(factor(Barcode), levels = barcodes) %>% pull(metacols[1]) -> metadat_sort
  Idents(so) <- factor(metadat_sort)
}

## Get colors from user parameter and add more if the default list is too short.
cols1 <- c("salmon1","steelblue","goldenrod1","lightskyblue3","firebrick2","dodgerblue3","darkseagreen3","darkorange1","mediumpurple4","lemonchiffon","orchid3","chocolate3","lightgrey")

# Make colored umap
p1 = DimPlot(so, reduction=reduction, 
             group.by=colname,
             pt.size=0.1) + 
  theme_classic() + 
  scale_color_manual(values=cols1) + 
  theme(legend.position="right") +
  guides(colour=guide_legend(ncol=1,override.aes = list(size = 2))) +
  ggtitle(colname)

# Generate Figure 5b:
Figure_5b <- print(p1)

###------------------------------------------------------------- DotPlot and Contingency Table ------------------------------------------------------------------------------------------------

so <- SO.sub 
Cluster_Identities <- read.csv("") # Load Cluster Identities.csv

so@meta.data %>% mutate(SCT_snn_res.2.4_mod = case_when(as.numeric(SCT_snn_res.2.8) == 32 ~ 21,                                                          as.numeric(SCT_snn_res.2.8) == 35 ~ 39,
                                                        as.numeric(SCT_snn_res.2.8) == 41 ~ 36,
                                                        as.numeric(SCT_snn_res.2.8) == 42 ~ 40,
                                                        TRUE ~ as.numeric(SCT_snn_res.2.4))) -> so@meta.data
Idents(so) <- so@meta.data$SCT_snn_res.2.4_mod

markers <- rev(c("Itgam","Ly6g","Adgre1","Cd68","Siglec1","Csf1r","Apoe","Mafb","Nr4a1",
             "F13a1","Lgals3","Cxcl9","Cxcl10","Sell","Arg1","Mmp13","Mmp12","Cd38","Mrc1",
             "Nos2","Cd40","Cd86","Cx3cr1","Cd163","Cd3e","Cd8a","Cd3d","Cd3g","Cd4",
             "Foxp3","Gzmk","Prf1","Il2ra","Ncr1","Ifng","Klrb1c","Ly6c1","Vcan",
             "Fn1","Ccr2","S100a9","S100a8","Il1b","G0s2","Csf3r","Cst3","Atox1","Nccrp1",
             "Ccr9","Siglech","Klk1","Cox6a2","Cd79a","Fcmr"))

dp <- DotPlot(so, assay="SCT", features=markers,dot.scale=4,cols = c("lightgrey", "blue"))

clustertab <- Cluster_Identities
dp$data$id <- as.character(dp$data$id)
clusname = clustertab$Latest_Clusnames
names(clusname) = as.character(clustertab$Clusnum)
dp$data$name <- clusname[dp$data$id]
dp$data$name <- factor(dp$data$name,levels=rev(c("NK","CD8","CD4","Tregs","M1_Mac_1","M1_Mac_2","M2_Mac","Mac_1","Mac_2","Mac_3","Mac_4","Mac_5","Mon_1","Mon_2","Mon_3","Mon_4","Mon_5","Mon_6","Mon_7","Mon_8","Mon_9","N_1","N_2","N_3","N_4","N_5","N_7","N_8","DC_1","DC_2","pDC","Unk_1","Unk_2","Unk_3","Unk_4","Unk_5","Unk_6","Unk_7","Unk_8","Unk_9")))

scale.func <- switch(EXPR = "radius", size = scale_size, radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))

plot <- ggplot(data = dp$data, mapping = aes_string(x = "features.plot", y = "name")) +
  geom_point(mapping = aes_string(size = "pct.exp", color = "avg.exp.scaled")) +
  scale.func(range = c(0, 4)) +
  theme_cowplot() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90))

# Generate Dotplot (Figure S3d)
Figure_S3d <- plot
Figure_S3d

# Generate Contingency Table (Table S4)
so@meta.data$cellclusters <- clusname[as.character(so@meta.data$SCT_snn_res.2.4_mod)]
Table_S4 <- table(so@meta.data$cellclusters,so@meta.data$Likely_CellType)
Table_S4

#write.table(Supplementary_Table4,"Table_S4.csv", sep = "\t", quote = FALSE, col.names = NA)

SO.sub <- so

rm(list=setdiff(ls(), c("SO.sub","SO","Figure_S3d","Table_S4","Figure_5b")))


# Generate Figure 3a:
SO_plot = SO.sub
doCiteSeq <- FALSE

summarizeCutOff <-min(5,20) 

samples = eval(parse(text=gsub('\\[\\]','c()','c("CD8dep","Combo","ENT","NHSIL12","PBS")')))
if (length(samples) == 0) {
  print("No samples specified. Using all samples...")
  samples = unique(SO_plot@meta.data$sample_name)
}

if("active.ident" %in% slotNames(SO_plot)){
  sample_name = as.factor(SO_plot@meta.data$orig.ident)
  names(sample_name)=names(SO_plot@active.ident)
  SO_plot@active.ident <- as.factor(vector())
  SO_plot@active.ident <- sample_name
  SO.sub_plot = subset(SO_plot, ident = samples)
} else {
  sample_name = as.factor(SO_plot@meta.data$orig.ident)
  names(sample_name)=names(SO_plot@active.ident)
  SO_plot@active.ident <- as.factor(vector())
  SO_plot@active.ident <- sample_name
  SO.sub_plot = subset(SO_plot, ident = samples)
}

print("selected SO:")
print(SO)

meta.df <- SO.sub_plot@meta.data
colnames(SO.sub_plot@meta.data) = gsub("\\.","_",colnames(SO.sub_plot@meta.data))

m = eval(parse(text=gsub('\\[\\]','c()','c("cellclusters")')))
m = m[!grepl("Barcode",m)]
if (length(m) == 0) {
  print("No metadata columns specified. Plotting sample names and RNA clusters...")
  x = colnames(SO.sub_plot@meta.data)
  x = x[grepl("RNA",x)]
  m = c("sample_name",x)
}

#ERROR CATCHING
#collect valid names of valid columns
validColumns <- character()
for (i in colnames(meta.df)) {
  if (!any(is.na(meta.df[[i]]))) {
    validColumns <-c(validColumns,i)
  }
}

colsToSummarize <- eval(parse(text=gsub('\\[\\]','c()','c()')))
m = unique(c(m,colsToSummarize))

if (length(colsToSummarize)>0) {
  #Optional Summarization of Metadata
  for (i in colsToSummarize) {
    col <- meta.df[[i]]
    valCount <- length(unique(col))
    
    if ((valCount >=summarizeCutOff) & (i != 'Barcode') & (!is.element(class(meta.df[[i]][1]),c("numeric","integer")))) {
      freqVals <- as.data.frame(-sort(-table(col)))$col[1:summarizeCutOff]
      print(freqVals)
      summarized_col = list()
      count <- 0
      for (j in col) {
        
        print(j)
        print(paste("count is",count))
        
        if (is.na(j) || is.null(j) || (j =="None")) {
          count <- count + 1
          summarized_col[count] <- "NULLorNA"
          print("NULLorNA")
        } else if (j %in% freqVals){
          count <- count + 1
          summarized_col[count] <- j
          print("valid")
        } else {
          count <- count + 1
          summarized_col[count] <- "Other"
          print("Other")
        }
      }
      meta.df[[i]] <- summarized_col
    }
  }
  #assign new metadata
  colnames(meta.df) = gsub("\\.","_",colnames(meta.df))
  SO.sub_plot@meta.data <- meta.df
  colnames(SO.sub_plot@meta.data) = gsub("\\.","_",colnames(SO.sub_plot@meta.data))
}

drawMetadata <- function(m){
  #check if there are NaNs in metadata, if there are, catch 
  if (any(is.na(meta.df[[m]]))) {
    print("ERROR: Metadata column appears to contain NA values. This is not recommended for clustering plots.")
    print("Please review your selected metadata column")
    print(head(meta.df[[m]]))
    print("Below are valid metadata to select for this plot:")
    print(validColumns)
    stop("End of error message.")
  }
  
  
  reduction = "umap"
  
  if (!(doCiteSeq)) {
    if(reduction=="tsne"){
      p1 <- DimPlot(SO.sub_plot, reduction = "tsne", group.by = "ident")
    } else if(reduction=="umap"){
      p1 <- DimPlot(SO.sub_plot, reduction = "umap", group.by = "ident")
    } else { 
      p1 <- DimPlot(SO.sub_plot, reduction = "pca", group.by = "ident")
    }
  } else {
    if(reduction=="tsne"){
      p1 <- DimPlot(SO.sub_plot, reduction = "protein_tsne", group.by = "ident")
    } else if(reduction=="umap"){
      p1 <- DimPlot(SO.sub_plot, reduction = "protein_umap", group.by = "ident")
    } else { 
      p1 <- DimPlot(SO.sub_plot, reduction = "protein_pca", group.by = "ident")
    }
  }
  
  #Categorical/Qualitative Variables
  if (!is.element(class(meta.df[[m]][1]),c("numeric","integer"))){
    if (!(doCiteSeq)) {
      #plot RNA clusters
      if(reduction=="tsne"){
        clusmat=data.frame(umap1=p1$data$tSNE_1,umap2=p1$data$tSNE_2, clusid=as.character(SO.sub_plot@meta.data[[m]]))
      } else if(reduction=="umap"){
        clusmat=data.frame(umap1=p1$data$UMAP_1,umap2=p1$data$UMAP_2, clusid=as.character(SO.sub_plot@meta.data[[m]]))
      } else { 
        clusmat=data.frame(umap1=p1$data$PC_1,umap2=p1$data$PC_2, clusid=as.character(SO.sub_plot@meta.data[[m]]))
      }
      
    } else {
      #else plot Antibody clusters
      if(reduction=="tsne"){
        clusmat=data.frame(umap1=p1$data$protein_tsne_1,umap2=p1$data$protein_tsne_2, clusid=as.character(SO.sub_plot@meta.data[[m]]))
      } else if(reduction=="umap"){
        clusmat=data.frame(umap1=p1$data$protein_umap_1,umap2=p1$data$protein_umap_2, clusid=as.character(SO.sub_plot@meta.data[[m]]))
      } else { 
        clusmat=data.frame(umap1=p1$data$protein_pca_1,umap2=p1$data$protein_pca_2, clusid=as.character(SO.sub_plot@meta.data[[m]]))
      }
    }
    
    clusmat %>% group_by(clusid) %>% summarise(umap1.mean=mean(umap1), umap2.mean=mean(umap2)) -> umap.pos
    title=as.character(m)
    cols=list()
    
    Lab.palette <- colorRampPalette(brewer.pal(12,"Paired"))
    n=length(unique((SO.sub_plot@meta.data[[m]])))
    cols[[1]]=brewer.pal(8, "Set3")[-2]  #Alternative
    cols[[2]]=brewer.pal(8, "Set1")
    cols[[3]]=c(cols[[1]],brewer.pal(8,"Set2")[3:6])
    cols[[4]]=c("#F8766D","#FF9912","#a100ff","#00BA38","#619CFF","#FF1493","#010407")
    cols[[5]]=c("blue","red","grey")
    cols[[6]]=Lab.palette(n)
    cols[[7]]=c("red","green","blue","orange","cyan","purple")
    cols[[8]]=c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4","#42d4f4","#f032e6","#bfef45","#fabebe","#469990","#e6beff","#9A6324","#800000","#aaffc3","#808000","#000075","#a9a9a9","#808080","#A9A9A9","#8B7355")
    colnum = 8
    
    n = length(unique(clusmat$clusid))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    qual_col_pals = qual_col_pals[c(7,6,2,1,8,3,4,5),]
    col = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    #Select to add labels to plot or not:
    if(TRUE){
      g <- ggplot(clusmat) +
        theme_bw() +
        theme(legend.title=element_blank()) +
        geom_point(aes(x=umap1, y=umap2,colour=clusid),size=0.01) +
        scale_color_manual(values=col) +
        xlab(paste(reduction,"-1")) + ylab(paste(reduction,"-2")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="right",
              panel.background = element_blank(), legend.text=element_text(size=rel(1))) +
        guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
        ggtitle(title) +
        geom_label_repel(data=umap.pos,aes(x=umap1.mean,y=umap2.mean,label=umap.pos$clusid),size=4)
    }
    else{
      g <- ggplot(clusmat, aes(x=umap1, y=umap2)) +
        theme_bw() +
        theme(legend.title=element_blank()) +
        geom_point(aes(colour=clusid),size=0.01) +
        scale_color_manual(values=col) +
        xlab(paste(reduction,"-1")) + ylab(paste(reduction,"-2")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="right",
              panel.background = element_blank(),legend.text=element_text(size=rel(1))) +
        guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
        ggtitle(title)    
    } 
    
  } else {
    # Plot quantitative data
    m = as.character(m)
    clusid = SO.sub_plot@meta.data[[m]]
    clusid = scales::rescale(SO.sub_plot@meta.data[[m]], to=c(0,1))
    clus.quant=quantile(clusid[clusid>0],probs=c(.1,.5,.9))
    midpt = clus.quant[2]
    midpt2 = clus.quant[1]
    midpt3 = clus.quant[3]
    
    if (!(doCiteSeq)) {
      #plot RNA clusters
      if(reduction=="tsne"){
        clusmat=data.frame(umap1=p1$data$tSNE_1,umap2=p1$data$tSNE_2, clusid=as.numeric(SO.sub_plot@meta.data[[m]]))
      } else if(reduction=="umap"){
        clusmat=data.frame(umap1=p1$data$UMAP_1,umap2=p1$data$UMAP_2, clusid=as.numeric(SO.sub_plot@meta.data[[m]]))
      } else { 
        clusmat=data.frame(umap1=p1$data$PC_1,umap2=p1$data$PC_2, clusid=as.numeric(SO.sub_plot@meta.data[[m]]))
      }
      
    } else {
      #else plot Antibody clusters
      if(reduction=="tsne"){
        clusmat=data.frame(umap1=p1$data$protein_tsne_1,umap2=p1$data$protein_tsne_2, clusid=as.numeric(SO.sub_plot@meta.data[[m]]))
      } else if(reduction=="umap"){
        clusmat=data.frame(umap1=p1$data$protein_umap_1,umap2=p1$data$protein_umap_2, clusid=as.numeric(SO.sub_plot@meta.data[[m]]))
      } else { 
        clusmat=data.frame(umap1=p1$data$protein_pca_1,umap2=p1$data$protein_pca_2, clusid=as.numeric(SO.sub_plot@meta.data[[m]]))
      }
    }
    
    clusmat %>% group_by(clusid) %>% summarise(umap1.mean=mean(umap1), umap2.mean=mean(umap2)) -> umap.pos
    title=as.character(m)
    print(environmentName(environment(arrange))) 
    clusmat %>% dplyr::arrange(clusid) -> clusmat
    print(environmentName(environment(arrange))) 
    g <- ggplot(clusmat, aes(x=umap1, y=umap2)) +
      theme_bw() +
      theme(legend.title=element_blank()) +
      geom_point(aes(colour=clusid),size=1) +
      scale_color_gradientn(colours = c("blue4", "lightgrey", "red"), values = scales::rescale(c(0, midpt2,midpt,midpt3, 1), limits = c(0, 1))) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
      ggtitle(title) +
      xlab("umap-1") + ylab("umap-2")
  }
  
  return(g)
}

grobs <- lapply(m, function(x) drawMetadata(x))

Figure_3a <- grid.arrange(grobs = grobs,ncol = 1,newpage=F)

rm(list=setdiff(ls(), c("SO.sub","SO","Figure_3a","Figure_S3d","Table_S4","Figure_5b")))
