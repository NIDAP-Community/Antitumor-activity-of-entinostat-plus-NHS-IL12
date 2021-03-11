# Load Libaries
library(SingleCellExperiment) # requires 1.4.1
library(sctransform) # requires 0.2.1
library(Seurat) # requires 3.1.5
library(ggplot2) # 3.3.3
library(RColorBrewer) # 1.1-2
library(scales) # 1.1.1
library(tidyverse) # 1.3.0
library(ggrepel) # 0.9.1
library(gdata) # 2.18.0
library(reshape2) # 1.4.4
library(tools) # 4.0.4
library(grid) # 4.0.4
library(gridBase) # 0.4-7
library(gridExtra) # 2.3
library(parallel) # 4.0.4
library(l2p) # requires version 0.0-1
library(magrittr) # requires 1.5
library(dplyr) # 1.0.5
library(tidyverse) # 1.3.0
library(RCurl) # 1.98-1.2
library(tidyr) # 1.1.3
library(stringr) # 1.4.0
library(MAST) # requires 1.8.2

# Requires SO generated from ModScore_Dotplot_ContingencyTable Code:
rm(list=setdiff(ls(), c("SO.sub","SO")))
Ortholog_Map_for_RNA_Seq <- read.csv("")
contrast <- c("Combo-PBS")
contrast_col <- "avg_logFC_Combo_vs_PBS"


# Generate Figures and Tables
DEG <- function(SO_input,CellType,contrasts = c("Combo-PBS")){

#collect parameters here:
contrasts <- contrasts
useSpark <- FALSE

samples = eval(parse(text=gsub('\\[\\]','c()','c("CD8dep","Combo","ENT","NHSIL12","PBS")')))

if (length(samples) == 0) {
  samples = unique(SO@meta.data$sample_name)
}

cell_vec <- rownames(SO_input@meta.data)[SO.sub@meta.data$Likely_CellType %in% CellType]
SO <- subset(SO_input, cells = cell_vec)

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

print("selected samples:")
print(SO.sub)

meta.df <- SO.sub@meta.data
  #read.csv(file = "/Users/bianjh/Documents/R files/Josh/Chariou/l2p_meta_test.csv")
colnames(SO.sub@meta.data) = gsub("\\.","_",colnames(SO.sub@meta.data))

#define contrasts
newcont <- list()
for (i in 1:length(contrasts)){
  newcont[[i]] <- c(paste(unlist(strsplit(contrasts[i],"-"))))
}
contrasts <- newcont

#ERROR CATCHING
#collect valid names of valid columns
validColumns <- character()
for (i in colnames(meta.df)) {
  if (!any(is.na(meta.df[[i]]))) {
    validColumns <-c(validColumns,i)
  }
}

param2test <- "orig_ident"

if (param2test =="") {
  mcols = colnames(SO.sub@meta.data)
  param2test <-mcols[grepl("RNA_snn",mcols)][[1]]
  print(paste("No parameter selected, defaulting to",param2test))
}

contrastTarget <- SO.sub@meta.data[[param2test]]
contrastType <- param2test
contrastCounts = as.data.frame(table(contrastTarget))
validContrasts = subset(contrastCounts, Freq>2)[[1]]

#catch malformed contrasts
for (i in contrasts) {
  if (!(i[[1]] %in% contrastTarget)) {
    print(paste(i[[1]],"is not a valid contrast for contrast type:", contrastType))
    print("Please see below for an example of valid contrasts for your selected contrast type.")
    print(validContrasts)
    stop("You have entered an invalid group to contrast against.")
  } else if (!(i[[2]] %in% contrastTarget) & (i[[2]] != "all")) {
    print(paste(i[[2]],"is not a valid contrast for contrast type:", contrastType))
    print("Please see below for an example of valid contrasts for your selected contrast type.")
    print(validContrasts)
    stop("You have entered an invalid group to contrast against.")
  } else if (length(i)>2) {
    print("Contrasts are as follows..")
    print(i)
    stop("The console says there are too many inputs in your contrasts. A contrast should only contain Group1-Group2, but the console thinks you have inputed Group1-Group2-Group3")
  } else if (!(i[[2]] %in% validContrasts) & (i[[2]] != "all")) {
    print(paste(i[[2]],"has two few values (less than 3 cells) to contrast against. Please see below for contrasts with enough cells:", validContrasts))
    stop("You have entered an invalid group to contrast against.")
  } else if (!(i[[1]] %in% validContrasts)) {
    print(paste(i[[1]],"has two few values (less than 3 cells) to contrast against. Please see below for contrasts with enough cells:", validContrasts))
    stop("You have entered an invalid group to contrast against.")
  }
}

#print out contrast cell contrastCounts
for (i in seq_along(contrasts)) {
  firstGroup <- contrasts[[i]][[1]]
  firstGroupCount <- subset(contrastCounts, contrastTarget == firstGroup)$Freq
  if  (contrasts[[i]][[2]]!= "all") {
    secondGroup <-contrasts[[i]][[2]]
    secondGroupCount <-subset(contrastCounts, contrastTarget == secondGroup)$Freq      
    print(paste("Contrast No.",i,"contrasts cluster",firstGroup,"with",firstGroupCount,"cells vs. cluster",secondGroup,"with",secondGroupCount,"cells."))
  } else {
    secondGroupCount <-ncol(SO.sub)-firstGroupCount
    print(paste("Contrast No.",i,"contrasts cluster",firstGroup,"with",firstGroupCount,"cells vs. all other clusters, totalling",secondGroupCount,"cells."))
  } 
}

#define and call function for running DEG
get_deg_table <- function(n) {
  library(Seurat)
  
  firstCluster <-n[[1]]
  secondCluster <- n[[2]]
  
  if (n[[2]]=='all') {
    secondCluster <- NULL
  }
  
  Idents(SO.sub) <- param2test
  markers = FindMarkers(SO.sub, ident.1 = firstCluster, ident.2 = secondCluster, test.use = "MAST", logfc.threshold = 0.25, verbose=FALSE, assay = "SCT")
  colnames(markers) <- chartr(old=" ",new="_",paste(colnames(markers), n[[1]],"vs",n[[2]],sep = "_"))
  return(markers)
}

if (useSpark) {
  deg_tables <- spark.lapply(contrasts, get_deg_table) 
} else {
  deg_tables <- lapply(contrasts, get_deg_table) 
}

for(i in seq_along(deg_tables)){
  degtab <- deg_tables[[i]]
  degtab %>% dplyr::filter(.[[1]] < 0.05) %>% dplyr::filter(.[[2]] > 0) %>% dim() -> pos 
  degtab %>% dplyr::filter(.[[1]] < 0.05) %>% dplyr::filter(.[[2]] < 0) %>% dim() -> neg
  print(paste0("The number of upregulated genes at p<0.05 in contrast number ", i, " is:"))
  print(pos[1])
  print(paste0("The number of downregulated genes at p<0.05 in contrast number ", i, " is:"))
  print(neg[1]) 
}

#Merge the deg tables together
out_df <- NULL
for (i in deg_tables) {
  if (is.null(out_df)) {
    out_df <- deg_tables[1]
    out_df <- as.data.frame(out_df)
  } else {
    out_df <- merge(out_df, i, by="row.names", all=TRUE)
    rownames(out_df) <- out_df$Row.names #set the rownames
    out_df$Row.names <- NULL #drop the row.names columns which we no longer need
  }
}

out_df$Gene <- rownames(out_df)
out_df$Row.names <- NULL
out_df <- out_df %>% dplyr::select(Gene, everything())
return(out_df)
}
l2p <- function(DEG_result,sort_descending,number_genes,contrast_type = "Combo-PBS"){

out_df <- DEG_result
if (contrast_type == "Combo-PBS"){
  out_df %>% select("Gene", "avg_logFC_Combo_vs_PBS") -> genesmat
  genesmat %>% arrange(desc(`avg_logFC_Combo_vs_PBS`)) -> genesmat
  genesmat %>% filter(!is.na(`avg_logFC_Combo_vs_PBS`)) -> genesmat
} else {
  out_df %>% select("Gene", "avg_logFC_CD8dep_vs_Combo") -> genesmat
  genesmat %>% arrange(desc(`avg_logFC_CD8dep_vs_Combo`)) -> genesmat
  genesmat %>% filter(!is.na(`avg_logFC_CD8dep_vs_Combo`)) -> genesmat
}
sort_descending = sort_descending

if (sort_descending) {
  genes_to_include = head(genesmat["Gene"], number_genes)
} else {
  genes_to_include = tail(genesmat["Gene"], number_genes)
}
genes_to_include <- as.vector(unique(unlist(genes_to_include)))
genes_universe = as.vector(unique(unlist(genesmat["Gene"])))
gene_set_sources_to_include = c("GO","KEGG","H","REACTOME")
categories_string <- paste(gene_set_sources_to_include, collapse=",")

organism = "Mouse"
if (organism != "Human") {
  organism_of_interest = "Mouse"
  orthology_table = filter(Ortholog_Map_for_RNA_Seq,     
                           Ortholog_Map_for_RNA_Seq$Organism==organism_of_interest)
  
  if ("from human"=="from human") {
    orthology_reference_column = "Human_gene_name"
    orthology_conversion_column = "Nonhuman_gene_name"
  } else {
    orthology_reference_column = "Nonhuman_gene_name"
    orthology_conversion_column = "Human_gene_name"
  }
  
  names(orthology_table)[names(orthology_table) == orthology_reference_column] <- "orthology_reference"
  names(orthology_table)[names(orthology_table) == orthology_conversion_column] <- "orthology_conversion"
  
  orthology_table %>% select("orthology_reference", "orthology_conversion") -> orthology_table
  orthology_table <- orthology_table
  orthology_table %>% filter(orthology_conversion %in% genes_to_include) %>% select(orthology_reference) -> genes_to_include
  genes_to_include <- as.vector(genes_to_include$orthology_reference)
  genes_to_include <- unique(genes_to_include)
  orthology_table %>% filter(orthology_conversion %in% genes_universe) %>% select(orthology_reference) -> genes_universe
  genes_universe <- as.vector(genes_universe$orthology_reference)
  genes_universe <- unique(genes_universe)
  
}else{
  genes_to_include = genes_to_include
  genes_universe = genes_universe
}

use_built_in_gene_universe = FALSE
if (use_built_in_gene_universe) {
  x <- l2pwcats(genes_to_include, categories_string)
  print("Using built-in gene universe.")
} else {
  x <- l2puwcats(genes_to_include, genes_universe, categories_string)
  print("Using all genes in differential expression analysis as gene universe.")
}
x %>%
  arrange(pval) %>% mutate(hitsPerc=(pwhitcount/(pwnohitcount+pwhitcount))*100) %>% mutate(pathtotal=pwnohitcount+pwhitcount) %>% filter(ratio >= 0) %>%
  select("pathwayname", "category", "pathwayaccessionidentifier", "pval", "fdr", "pwhitcount", "genesinpathway", "pwnohitcount","pathtotal","hitsPerc", "inputcount", "pwuniverseminuslist","ratio")  %>% dplyr::rename(diff_ratio = ratio) %>% dplyr::filter(pval < 0.05)%>% dplyr::filter(pwhitcount >= 5) -> x
#allgenes <- lapply(x$pathwayaccessionidentifier,function(x) l2pgetgenes4acc(x)) 
#x %>% mutate("allgenes" = paste(allgenes,sep="",collapse=" ")) -> x
print(paste0("Total number of pathways: ", nrow(x)))
goResults <- x
goResults %>% top_n(10, wt=-log(pval)) %>%
  arrange(-log(pval)) -> goResults
minp = min(goResults$pval) - 0.1*min(goResults$pval)
maxp = max(goResults$pval) + 0.1*max(goResults$pval)
#print(goResults$pval)
sizemax = ceiling(max(goResults$pwhitcount)/10)*10  
goResults %>% mutate(pathwayname2 = stringr::str_replace_all(pathwayname, "_", " ")) -> goResults
goResults %>% mutate(pathwayname2 = stringr::str_wrap(pathwayname2,30)) -> goResults

if (FALSE){
  goResults %>% dplyr::mutate(percorder = order(goResults$pval)) -> goResults
  goResults$pathwayname2 <- factor(goResults$pathwayname2, levels = goResults$pathwayname2[goResults$percorder])
  xmin = floor(min(goResults$pval))
  xmax = max(goResults$pval) 
  gplot <- goResults %>% 
    ggplot(aes(x=pval,
               y=pathwayname2, 
               colour=hitsPerc, 
               size=pwhitcount)) +
    geom_point() +
    theme(text = element_text(size=8)) +
    xlim(xmin,xmax) +
    expand_limits(colour = seq(minp, maxp, by = 10),
                  size = seq(0, sizemax,by=10)) +
    labs(x="p value", y="GO term", colour="Hits (%)", size="Count") 
  print(gplot)
}else{
  goResults %>% mutate(percorder = order(goResults$hitsPerc)) -> goResults
  goResults$pathwayname2 <- factor(goResults$pathwayname2, levels = goResults$pathwayname2[goResults$percorder])
  xmin = floor(min(goResults$hitsPerc)-5)
  xmax = ceiling(max(goResults$hitsPerc)+5) 
  gplot <- goResults %>% 
    ggplot(aes(x=hitsPerc,
               y=pathwayname2, 
               colour=pval, 
               size=pwhitcount)) +
    geom_point() +
    theme_classic() +
    theme(text = element_text(size=20)) +
    xlim(xmin,xmax) +
    expand_limits(colour = seq(minp, maxp, by = 10),
                  size = seq(0, sizemax,by=10)) +
    labs(x="Hits (%)", y="GO term", colour="p value", size="Count") 
  print(gplot)
}
rownames(x) <- 1:nrow(x)
return(x)
}
Bubbleplot <- function(l2p_result,top10 = TRUE,select_rows){
  
  goResults <- l2p_result
  
  if(top10 == TRUE){
    goResults %>% top_n(10, wt=-log(pval)) %>% 
      mutate(hitsPerc=hitsPerc*100) %>%
      arrange(-log(pval)) -> goResults
  } else{
    goResults %>% dplyr::slice(select_rows) %>% 
      mutate(hitsPerc=hitsPerc) -> goResults
  }
  xmin = floor(min(goResults$hitsPerc)-5)
  xmax = ceiling(max(goResults$hitsPerc)+5) 
  sizemax = ceiling(max(goResults$pwhitcount)/10)*10  
  
  goResults %>% mutate(pathwayname2 = str_replace_all(pathwayname, "_", " ")) -> goResults
  goResults$pathwayname2 <- str_to_upper(goResults$pathwayname2)
  goResults %>% mutate(pathwayname2 = stringr::str_wrap(pathwayname2,30)) -> goResults
  goResults %>% mutate(percorder = order(goResults$hitsPerc)) -> goResults
  goResults$pathwayname2 <- factor(goResults$pathwayname2, levels = goResults$pathwayname2[goResults$percorder])
  print(class(goResults$pathwayname2))
  
  minp = min(goResults$pval) - 0.1*min(goResults$pval)
  maxp = max(goResults$pval) + 0.1*max(goResults$pval)
  
  gplot <- goResults %>% 
    ggplot(aes(x=hitsPerc,
               y=pathwayname2, 
               colour=pval, 
               size=pwhitcount)) +
    geom_point() +
    theme_classic() +
    theme(text = element_text(size=20)) +
    xlim(xmin,xmax) +
    expand_limits(colour = seq(minp, maxp, by = 10),
                  size = seq(0, sizemax,by=10)) +
    theme(legend.position = "right", legend.key.height = unit(1, "cm")) +
    labs(x="Hits (%)", y="Pathway", colour="p value", size="Count") +
    theme(axis.title.x = element_text(size = rel(1.2)),axis.title.y = element_text(size = rel(1.2)))
  
  print(gplot)
  return(goResults) 
}   

DEG_4a <- DEG(SO,"CD8_T")
Table_S5 <- l2p(DEG_4a,TRUE,215)
Figure_4a <- Bubbleplot(Table_S5,top10 = FALSE,select_rows = c(2,7,11,20,31,40,42,46:47,54,69,70,74:76,78))

# Note l2p will display both Figure and Table with the exception with Figure 4a, where clinically relevant pathways are manually selected
DEG_5f <- DEG(SO,"PMN_neutrophils")
Figure_5f_Table_S6 <- l2p(DEG_5f,TRUE,156)

DEG_6b <- DEG(SO,"Macrophages")
Figure_6b_Table_S8 <- l2p(DEG_6b,TRUE,202)

DEG_6d <- DEG(SO,"M1")
Figure_6d_Table_S9 <- l2p(DEG_6d,TRUE,104)

DEG_7d <- DEG(SO,CellType = c("Macrophages","M1","M2"),contrasts = c("CD8dep-Combo")) 
Figure_7d_Table_S10 <- l2p(DEG_7d,FALSE,133,contrast_type = "CD8dep-Combo")

DEG_S7k <- DEG(SO,"PMN_neutrophils",contrasts = c("CD8dep-Combo"))
Figure_7d_Table_S11 <- l2p(DEG_S7k,FALSE,82,contrast_type = "CD8dep-Combo")


