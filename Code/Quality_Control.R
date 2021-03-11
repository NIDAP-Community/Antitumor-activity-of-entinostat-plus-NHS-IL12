# Install/Load relevant packages:

library(Seurat) # 3.1.5
library(rhdf5) # 2.34.0
library(reshape2) # 1.4.4
library(tidyverse) # 1.3.0
library(gridExtra) # 2.3
library(RColorBrewer) # 1.1-2
library(stringr) # 1.4.0
library(ggplot2) # 3.3.3
library(gtable) # 0.3.0
library(sctransform) # 0.2.1

### Step 1: Filter & Generate Quality Control Plots from h5 files

localFilePaths <- list.files(path = "../Data/", pattern = ".h5", full.names = TRUE)

subsetRegex = eval(parse(text=gsub('\\[\\]','c()','c()')))
Keep <- TRUE
if (length(subsetRegex) > 0) {
  if (Keep == TRUE){
    for (i in length(subsetRegex)) {
      localFilePaths <- localFilePaths[grepl(subsetRegex[[i]],localFilePaths)]
    }
  }
  else{
    for (i in length(subsetRegex)) {
      localFilePaths <- localFilePaths[!grepl(subsetRegex[[i]],localFilePaths)]
    }
  }
}

obj.list <- lapply(localFilePaths, function(x) { return(Read10X_h5(x, use.names=TRUE)) })
rename = TRUE
if (rename == FALSE){
  names(obj.list) <- lapply(localFilePaths, basename)
  names(obj.list) <- sapply(names(obj.list), function(x) gsub("_filtered(\\w+)?.h5","", x))
  names(obj.list) <- sapply(names(obj.list), function(x) gsub("\\.(\\w+)?","", x))
} else{
  names(obj.list) <- c("PBS","ENT","NHSIL12","Combo","CD8dep")
  obj.list <- obj.list[sort(names(obj.list))]
}

mincells = 3
mingenes = 200
organism = "Mouse"

if (organism == "Human"){
  mitoch = "^MT-"
} else{
  mitoch = "^mt-"
  cc.genes$g2m.genes= str_to_title(cc.genes$g2m.genes)
  cc.genes$s.genes = str_to_title(cc.genes$s.genes)
}

seurat_object <- function(i) {
  if (class(obj.list[[i]]) == "dgCMatrix"){
    so.nf <- CreateSeuratObject(counts = obj.list[[i]], assay = "RNA", project=names(obj.list)[[i]], min.cells = 1, min.features = 0)
  } else {
    so.nf <- CreateSeuratObject(counts = obj.list[[i]][1]$`Gene Expression`, assay = "RNA", project=names(obj.list)[[i]], min.cells = 1, min.features = 0)
  }
  so.nf <- NormalizeData(so.nf, normalization.method = "LogNormalize", scale.factor = 10000)
  so.nf[["percent.mt"]] <- PercentageFeatureSet(object = so.nf, pattern = mitoch)
  so.nf$log10GenesPerUMI <- log10(so.nf$nFeature_RNA) / log10(so.nf$nCount_RNA)
  
  if (FALSE){
    antibodies = rownames(obj.list[[i]][2]$`Antibody Capture`)[!grepl("HTO*",rownames(obj.list[[i]][2]$`Antibody Capture`))]
    so.nf[['Protein']] <- CreateAssayObject(counts = obj.list[[i]][2]$`Antibody Capture`[antibodies, colnames(x = so.nf)])
    so.nf <- NormalizeData(so.nf, assay = "Protein", normalization.method = "CLR")
  }
  if (FALSE){
    antibodies = rownames(obj.list[[i]][2]$`Antibody Capture`)[!grepl("HTO*",rownames(obj.list[[i]][2]$`Antibody Capture`))]
    HTO = rownames(obj.list[[i]][2]$`Antibody Capture`)[grepl("HTO*",rownames(obj.list[[i]][2]$`Antibody Capture`))]
    so.nf[['HTO']] <- CreateAssayObject(counts = obj.list[[i]][2]$`Antibody Capture`[HTO, colnames(x = so.nf)])
    so.nf <- NormalizeData(so.nf, assay = "HTO", normalization.method = "CLR")
  }
  
  #Filtered Seurat Object:
  if (class(obj.list[[i]]) == "dgCMatrix"){
    so <- CreateSeuratObject(counts = obj.list[[i]], assay = "RNA", project=names(obj.list)[[i]], min.cells = mincells, min.features = mingenes)
  } else {
    so <- CreateSeuratObject(counts = obj.list[[i]][1]$`Gene Expression`, assay = "RNA", project=names(obj.list)[[i]], min.cells = mincells, min.features = mingenes)
  }
  so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
  so[["percent.mt"]] <- PercentageFeatureSet(object = so, pattern = mitoch)
  so$log10GenesPerUMI <- log10(so$nFeature_RNA) / log10(so$nCount_RNA)
  
  if (FALSE){
    rownames(obj.list[[i]][2]$`Antibody Capture`) <- gsub(pattern = "_TotalSeqC", replacement = "", rownames(obj.list[[i]][2]$`Antibody Capture`))
    antibodies = rownames(obj.list[[i]][2]$`Antibody Capture`)[!grepl("HTO*",rownames(obj.list[[i]][2]$`Antibody Capture`))]
    so[['Protein']] <- CreateAssayObject(counts = obj.list[[i]][2]$`Antibody Capture`[antibodies, colnames(x = so)])
    so <- NormalizeData(so, assay = "Protein", normalization.method = "CLR")
  }
  if (FALSE){
    rownames(obj.list[[i]][2]$`Antibody Capture`) <- gsub(pattern = "_TotalSeqC", replacement = "", rownames(obj.list[[i]][2]$`Antibody Capture`))
    HTO = rownames(obj.list[[i]][2]$`Antibody Capture`)[grepl("HTO*",rownames(obj.list[[i]][2]$`Antibody Capture`))]
    so[['HTO']] <- CreateAssayObject(counts = obj.list[[i]][2]$`Antibody Capture`[HTO, colnames(x = so)])
    so <- NormalizeData(so, assay = "HTO", normalization.method = "CLR")
  }
  if (FALSE) {
    allGenes = rownames(so)
    VDJgenes = c("TRBV","TRAV","TRBD","TRAJ","TRBJ")
    print("Removing VDJ genes. Genes removed...")
    for (j in VDJgenes) {
      print(allGenes[grepl(j, allGenes)])
      allGenes = allGenes[!grepl(j, allGenes)]  
    }
    so <- SubsetData(so,features = allGenes,assay="RNA")
  }
  
  cat("\n\n")
  cat(names(obj.list)[i],":\n")
  so.origcount = dim(so)[2]
  cat(paste0("Original Cell Count=", so.origcount),"\n")
  
  #Start with filtering here:
  maxgenes = 2500
  complexity = 0.5
  MAD_gene <- TRUE
  ngenestdev <- mad(so@meta.data$nFeature_RNA)
  ngenemed <- median(so@meta.data$nFeature_RNA)
  ngenemaxlim <- ngenemed+(3*ngenestdev)
  gl = format(round(ngenemaxlim,0),nsmall=0)
  
  maxmitoch = 25
  
  MAD_mitoch <- FALSE
  mitostdev <- mad(so@meta.data$percent.mt)
  mitomed <- median(so@meta.data$percent.mt)
  mitomaxlim <- mitomed+(3*mitostdev)
  ml = format(round(mitomaxlim,2),nsmall=2)
  
  if (MAD_gene == TRUE & MAD_mitoch == TRUE)       {
    cat(paste0("Gene Count Filter = low:",mingenes," high:",gl),"\n")
    cat(paste0("Mitochondrial Percentage Filter =",ml,"\n"))
    cat(paste0("Complexity Filter =",complexity,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < ngenemaxlim & so@meta.data$percent.mt < mitomaxlim & so@meta.data$log10GenesPerUMI > complexity), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
  }
  else if (MAD_gene == FALSE & MAD_mitoch == TRUE) {
    cat(paste0("Gene Count Filter = low:", mingenes," high:", maxgenes),"\n")
    cat(paste0("Mitochondrial Percentage Filter =",ml,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < maxgenes & so@meta.data$percent.mt < mitomaxlim & so@meta.data$log10GenesPerUMI > complexity), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
  }
  else if (MAD_gene == TRUE & MAD_mitoch == FALSE){
    cat(paste0("Gene Count Filter = low:",mingenes," high:",gl),"\n")
    cat(paste0("Mitochondrial Percentage Filter =", maxmitoch,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < ngenemaxlim & so@meta.data$nFeature_RNA > mingenes & so@meta.data$percent.mt < maxmitoch & so@meta.data$log10GenesPerUMI > complexity), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
  }
  else {
    cat(paste0("Gene Count Filter = low:", mingenes," high:", maxgenes),"\n")
    cat(paste0("Mitochondrial Percentage Filter =", maxmitoch,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < maxgenes & so@meta.data$nFeature_RNA > mingenes & so@meta.data$percent.mt < maxmitoch & so@meta.data$log10GenesPerUMI > complexity), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain),"\n\n")
  }
  
  plothist <- function(count.df,name){
    g=ggplot(count.df,aes(x=value,fill=filt)) + 
      theme_bw() +
      geom_histogram(binwidth=.05, alpha = 0.7, position="identity") +
      scale_fill_manual(values = c("#00AFBB", "#FC4E07")) +
      scale_color_manual(values = c("#00AFBB", "#FC4E07")) +
      labs(x = NULL) +
      theme(plot.title = element_text(size=6),legend.position='right',legend.text=element_text(size=10),
            legend.title=element_blank()) + 
      ggtitle(paste(name,count.df$variable[1])) +
      scale_x_continuous(trans='log10') + 
      scale_linetype_manual(values=rep(c('solid', 'dashed','dotted'),6))
    return(g)
  }
  
  plotviolin <- function(count.df,name){
    axislab = unique(count.df$filt)
    col1=brewer.pal(8, "Set3")[-2] 
    col2=c(col1,brewer.pal(8,"Set2")[3:6])
    
    v = ggplot(count.df, aes(x=filt, y=value)) +
      ggtitle(paste(name,count.df$variable[1])) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),legend.text=element_text(size=rel(1.5)),
            legend.title=element_blank(), axis.text=element_text(size=10),
            axis.title.x=element_blank(),axis.title.y=element_blank(),
            axis.text.x=element_text(angle=45,hjust=1),
            plot.title = element_text(size = 12, face = "bold")) +
      geom_violin(aes(fill=as.factor(filt))) +  
      scale_fill_manual(values = c("#00AFBB", "#FC4E07")) + 
      geom_boxplot(width=.1) +
      scale_x_discrete(limits = as.vector(axislab)) 
    return(v)
  }
  
  Runplots <- function(x,name){
    df.m %>% dplyr::filter(variable == x) -> count.df
    df2.m %>% dplyr::filter(variable == x) -> count2.df
    qc.df <- array(0,dim=c(0,4))
    qc.df <- rbind(qc.df,count2.df,count.df)
    if(FALSE){
      gg <- plothist(qc.df,name)}
    else{
      gg <- plotviolin(qc.df,name)
    }
  }
  
  RunScatter <- function(x,name){
    x <- as.character(x)
    scplot.m = so@meta.data %>% dplyr::select("nCount_RNA",x) %>% dplyr::mutate(filt = "filt")
    scplot2.m = so.nf@meta.data %>% dplyr::select("nCount_RNA",x) %>% dplyr::mutate(filt = "raw") 
    sc.plot.all = rbind(scplot2.m,scplot.m)
    g=ggplot(sc.plot.all,aes_string(x="nCount_RNA",y=x,color="filt")) + 
      geom_point(size = 0.5) + 
      theme_classic() +
      ggtitle(paste(name)) 
    return(g)
  }
  
  df.m <- melt(so@meta.data)
  df.m$filt <- "filt"
  df.m$filt <- as.factor(df.m$filt)
  df2.m <- melt(so.nf@meta.data)
  df2.m$filt <- "raw"
  df2.m$filt <- as.factor(df2.m$filt)
  
  v <- unique(df.m$variable)
  grob.list <- lapply(v,function(x){Runplots(x,so@project.name)})
  grob2.list <- lapply(v,function(x){RunScatter(x, so@project.name)})
  grob.all <- arrangeGrob(grobs = grob.list, ncol = length(grob.list))
  grob2.all <- arrangeGrob(grobs = grob2.list, ncol = length(grob2.list))
  so2.list <- list(so,so.nf,grob.all,grob2.all)
  
  return(so2.list)
}


so.list <- lapply(seq_along(obj.list), seurat_object)

so.f.list <- lapply(so.list,function(x) x[[1]])
names(so.f.list) <- sapply(names(obj.list), function(x) gsub("_filtered.h5", "", x))

so.nf.list <- lapply(so.list,function(x) x[[2]])
names(so.nf.list) <- sapply(names(obj.list), function(x) gsub("_filtered.h5", "", x))

so.final.list <<- so.f.list

cat("Final filtered samples:\n")
print(so.f.list)
cat("Final filtered sample names:\n")
print(names(so.f.list))   

rm(list=setdiff(ls(), c("so.final.list")))

#-------------------------------------------------------------------------------------------------------------------------------

### Generate Post-filter Data

obj.list <- so.final.list
qc.df <- array(0,dim=c(0,3))
for (i in 1:length(obj.list)){
  so <- obj.list[[i]]
  print(so)
  df.m <- melt(so@meta.data)
  qc.df <- rbind(qc.df,df.m)
}

qfilter <- function(x){
  library(dplyr)
  qc.df %>% dplyr::filter(variable == x)
}

qc.count <- lapply(unique(qc.df$variable), function(x) {qfilter(x)})
qc.count[[1]] %>% dplyr::filter(variable=="nCount_RNA") %>% pull(value) -> RNAcounts

so@meta.data %>% rownames_to_column("Barcode") -> meta.df

rm(list=setdiff(ls(), c("obj.list")))

#--------------------------------------------------------------------------------------------------------------------------------------------------------------

### PCA & Normalization

SO = obj.list

if (class(SO) =="Seurat") {
  x =list()
  x[[1]] <- SO
  SO <- x
}
vars_to_regress <- c()
doJackStraw <- FALSE
npcs = 30
linearScale = FALSE

# Linearly scale data without regressing anything.
scale_so <- function(so){
  so <- CellCycleScoring(object = so, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
  so$CC.Difference <- so$S.Score - so$G2M.Score
  so <- FindVariableFeatures(object = so, nfeatures = 2000, mean.cutoff = c(1, 8), dispersion.cutoff = c(1, 100000), selection.method = "vst")
  all.genes <- rownames(so)
  so <- ScaleData(so,features=all.genes)
  return(so)
}

# Make PCA without regressing anything, and using only SCTransform().
pca_noregress <- function(so) {
  so <- SCTransform(so,do.correct.umi = FALSE,return.only.var.genes = FALSE)
  so <- RunPCA(object = so, features = VariableFeatures(object = so), npcs = npcs)
  return(so)
}

# Make PCA with SCTransform() and optional ScaleData, and do so with both regression (if user requests) and on all genes.
pca <- function(so) {
  if(linearScale == TRUE) {
    all.genes <- rownames(so)
    if(is.null(vars_to_regress)){     
      so <- so
    }
    else{
      so <- ScaleData(so, features=all.genes, vars.to.regress = vars_to_regress) 
    }
  }
  # Run SCTransform().
  if(is.null(vars_to_regress)){
    so <- so
  }
  else { 
    so <- SCTransform(so,do.correct.umi = TRUE, vars.to.regress = vars_to_regress, return.only.var.genes = FALSE)
  }
  so <- RunPCA(object = so, npcs = npcs)
  return(so)
}

so_scale <- lapply(SO, scale_so) 
so_list <- lapply(so_scale, pca) 

rm(list=setdiff(ls(), c("so_list")))

#-------------------------------------------------------------------------------------------------------------------

### Combine and Renormalize

SO = so_list

doMergeData <- !FALSE
dat = vector()

integratedata = FALSE

if (length(SO) > 1) {
  for(i in 2:length(SO)){dat=c(dat,SO[[i]]) }
  SO_merge <- merge(SO[[1]], y = dat, add.cell.ids = names(SO), project = "scRNAProject", merge.data = TRUE)
  allgenes <- rownames(SO_merge)
  SO_merge <- ScaleData(SO_merge, assay = "RNA", features=allgenes)
} else {
  SO_merge <- SO[[1]]
  allgenes <- rownames(SO_merge)
  SO_merge <- ScaleData(SO_merge, assay = "RNA", features=allgenes)
}

if (!("orig.ident" %in% colnames(SO_merge@meta.data))) {
  SO_merge@meta.data$orig.ident <- SO_merge@meta.data$orig_ident
}

if ("Protein" %in% names(SO_merge@assays)){
  doCiteSeq <-TRUE
}

if(FALSE){
  SO_merge <- ScaleData(SO_merge, assay = "HTO")
}

# Regress Mitochondrial Content
npcs = 21
Do_SCTransform = TRUE
vars_to_regress = c("percent.mt")

if (Do_SCTransform){
  if(is.null(vars_to_regress)){
    SO_merge <- SCTransform(SO_merge,do.correct.umi = TRUE, return.only.var.genes = FALSE)}
  else{       
    SO_merge <- SCTransform(SO_merge,do.correct.umi = TRUE,vars.to.regress=vars_to_regress,return.only.var.genes = FALSE,verbose = TRUE) 
  }
}else{
  all.genes <- rownames(SO_merge)
  if(is.null(vars_to_regress)){
    SO_merge <- SO_merge
  }
  else{
    SO_merge <- ScaleData(SO_merge, features=all.genes, assay = "RNA", vars.to.regress=vars_to_regress) 
  }
  DefaultAssay(SO_merge) <- "RNA"   
}

if (length(SO)>1) {
  all_features <- lapply(SO, row.names) %>% Reduce(intersect, .)
  if(integratedata==TRUE){
    integ_features <- SelectIntegrationFeatures(object.list = SO, nfeatures = 3000) 
    if(!is.null(SO[[1]]@assays$SCT)){
      SO <- PrepSCTIntegration(object.list = SO, anchor.features = integ_features)
      k.filter <- min(200, min(sapply(SO, ncol)))
      integ_anchors <- FindIntegrationAnchors(object.list = SO, normalization.method = "SCT", k.filter=k.filter, anchor.features = integ_features)
      SO_merge <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT",features.to.integrate = all_features)
      SO_merge <- ScaleData(SO_merge,features=all_features)
    }
    else{
      k.filter <- min(200, min(sapply(SO, ncol)))
      integ_anchors <- FindIntegrationAnchors(object.list = SO, k.filter=k.filter, anchor.features = integ_features)
      SO_merge <- IntegrateData(anchorset = integ_anchors,features.to.integrate = all_features)
      SO_merge <- ScaleData(SO_merge,features=all_features)  
    }}
}

SO_merge <- FindVariableFeatures(object = SO_merge, nfeatures = 2000, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, 100000), selection.method = "vst", verbose = FALSE)
SO_merge <- RunPCA(object = SO_merge, npcs = npcs, verbose = FALSE,seed.use = 42)
SO_merge <- RunUMAP(object = SO_merge, reduction = "pca", dims = 1:npcs, seed.use=42)
SO_merge <- RunTSNE(object = SO_merge, reduction = "pca", dim.embed = 2, dims = 1:npcs, seed.use = 1)
SO_merge <- FindNeighbors(SO_merge, dims = 1:npcs)

for (i in seq(2.4,2.8,0.2)){
  SO_merge <- FindClusters(SO_merge, resolution = i, algorithm = 2)
}

# Save final SO from QC Steps
rm(list=setdiff(ls(), c("SO_merge")))


