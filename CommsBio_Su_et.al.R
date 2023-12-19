################################################################################
################################################################################
####                                                                        ####
####                           Single-Cell Analysis                         ####
####                                                                        ####
################################################################################
################################################################################
set.seed(11051991)





#### Load library ####
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(gprofiler2)
library(stringr)
library(RColorBrewer)
library(scales)

library(ggpubr)
library(ggdendro)
library(cowplot)
library(magick)
library(purrr)
library(raster)
library(EnhancedVolcano)
library(gridExtra)
library(grid)
library(gridtext)
library(gtable)
library(officer)
library(ggrastr)
library(ggvenn)

library(monocle3)
library(NMF)
library(snow)
library(survival)
library(survminer)




#### Common input ####
object.name <- c("Ben", "Low1", "Low2", "Med", "High_2", "High_1", "COS")
file.name <- c("L49", "L07", "L28", "L63", "L31", "L44", "L43")


# Standard intigration by MNN
mnnCT7.raw.pc <- 50
mnnCT7.raw.dims <- 40
mnnCT7.raw.res <- 0.2

mnnCT7.pc <- 50
mnnCT7.dims <- 40
mnnCT7.res <- 0.2


# Standard MNN reintegration
mnnCT7.3.pc <- 50
mnnCT7.3.dims <- 40
mnnCT7.3.res <- 0.2




#### Color vector ####
CT7.cluster <- c("Chon1", "High2", "Stro", "Cos", "High1", "Chon2", "Leuk", "Prol")
names(CT7.cluster) <- 0:7
CT7.cluster.names <- c(c(0:7), CT7.cluster)
col.CT7 <- c(brewer.pal(length(CT7.cluster), "Set1"), brewer.pal(length(CT7.cluster), "Set1"))
names(col.CT7) = CT7.cluster.names
col.CT7

col.CT7.old <- col.CT7
col.CT7[c("5", "Chon2")] <- "#FEE715FF"


# Subcluster
mnnCT7.sucluster <- c("Chon1", "Chon2", "High1a", "High1b", "High2a", "High2b", 
                      "High1", "High2", "Cos", "nStro", "Prol", "Stro1", "Leuk")


# Suclusters
col.sub <- c(`High2a`="royalblue", `High2b`="skyblue", `Hype1`="purple", `Hype2`="slateblue", 
             `High1a`="orangered", `High1b`="orange", `nStro` = "green", col.CT7, `Stro1` = col.CT7[["Stro"]])


# CNV cluster
cluster.cnv <- c("Chon1", "Chon2", "High1", "High2", "Cos", "Prol", "Stro1", "nStro", "Leuk")
col.allcnv <- c(col.CT7, col.sub)[cluster.cnv]


# Patient
patient <- c("Ben", "Low_1", "Low_2", "Med", "High_2", "High_1", "COS")
col.patient <- brewer.pal(length(patient), "Set2")
names(col.patient) = patient
col.patient


# All colors
col.all <- c(col.sub, col.patient)


# Subcluster 2
CT7.2.cluster.id <- c(0:6)
col.CT7.2 <- brewer.pal(length(CT7.2.cluster.id), "Dark2")
names(col.CT7.2) = CT7.2.cluster.id
col.CT7.2


# Stromal cell type cluster
mnnCT7.2.cluster <- setNames(c("miStro", "osStro", "gCAF", "FC", "CAF", "vSMC", "EC"), 0:6)
col.CT7.2 <- c(col.CT7.2, col.CT7.2)
names(col.CT7.2)[1:length(mnnCT7.2.cluster)] <- mnnCT7.2.cluster
col.CT7.2


# Subcluster 6
CT7.6.cluster.id <- c(0:4)
col.CT7.6 <- brewer.pal(length(CT7.6.cluster.id), "Set3")
names(col.CT7.6) = CT7.6.cluster.id
col.CT7.6


# Immune cell type cluster
mnnCT7.6.cluster <- c("M1", "M2", "TC", "SP", "OC")
col.CT7.6 <- c(col.CT7.6, col.CT7.6)
names(col.CT7.6)[1:length(mnnCT7.6.cluster)] <- mnnCT7.6.cluster
col.CT7.6


# Fetal femur
FF.cluster <- setNames(c("TC", "RC", "HC", "PC", "fibroblast", "endothelial cell"), 0:5)
FF.cluster <- FF.cluster[c(2,1,3:6)]
col.ff <- scales::hue_pal()(length(FF.cluster))
col.ff <- c(col.ff, col.ff)
names(col.ff) <- c(FF.cluster, names(FF.cluster))
col.ff


# Histology
his <- c("benign", "low", "medium", "high", "dedifferentiated")
his.col <- brewer.pal(6, "Spectral")[5:1]
names(his.col) <- his
his.col


# Mutation col
mut <- c("TRUE", "FALSE", "#N/A")
mut.col <- c("red", "green", "gray")
names(mut.col) <- mut
mut.col


# Bulk group col
group.name <- c("5" = "Ben_bulk", "3" = "Med_bulk", "4" = "High_bulk", "2" = "Ded_bulk" , "1" = "Stro_bulk")
col.group <- brewer.pal(length(group.name), "Dark2")
names(col.group) <- group.name
col.group


# Heatmap col
ann_colors <- list(COL2A1.mut = mut.col, P53.mutation = mut.col, IDH.mut = mut.col, 
                   histology = his.col, group = col.group)


# Chromosome color
col.chr <- c(brewer.pal(8, "Pastel1")[-6], brewer.pal(8, "Pastel2")[-6], brewer.pal(9, "Set3")[-2], "#000000")
names(col.chr) <- c(1:22, "anno")
# col.chr[-23] <- adjustcolor(col.chr[-23], alpha.f = 0.5)
col.chr


# CNV group color
group.cnv <- c("Norm", "Low", "Med", "High")
col.cnv <-  brewer.pal(6, "Spectral")[5:2]
names(col.cnv) <- group.cnv
col.cnv


# NEW color
NEW.cluster <- c(paste0("L80-", 1:3), paste0("L81-", c(1,3,2)), paste0("L96-", 1:5), "Stro", "Leuk")
col.NEW <- c(brewer.pal(length(NEW.cluster)-2, "Set3"), col.CT7[c("Stro", "Leuk")] )
names(col.NEW) <- NEW.cluster
col.NEW




#### Functions ####
#### My Heatmap
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)




my_heatmap <- function(object, genes, bar = c("seurat_clusters", "orig.ident"), 
                       row_split = NULL,
                       column_title = NULL,
                       row_title = NULL,
                       show_rownames = T, 
                       show_annotation_name = F,
                       right_annotation = NULL,
                       cols = list("x" = setNames("x", "x")), 
                       min = -2.5, max = 2.5, fontsize = 12, res = 1){
  
  mdata <- GetAssayData(object, slot = "scale.data")[genes, ]
  # mdata <- FetchData(object, genes, slot = "scale.data") %>% t
  
  # Scale min and max values
  pmed <- function(x) pmin(pmax(x, min), max)
  mdata <- apply(mdata, 2, pmed)
  
  
  # Order cells according to input features
  if(length(bar) == 1) {
    bar_col <- object@meta.data[,c(bar, bar)]
    bar_col <- bar_col[do.call(order, bar_col), ][1]
  }else {
    bar_col <- object@meta.data[, bar]
    bar_col <- bar_col[do.call(order, bar_col), ]}
  
  mdata <- mdata[ , rownames(bar_col)]
  
  # # Add gaps_col
  # gaps_col <- bar_col %>% group_by(bar_col[, 1]) %>% tibble::rowid_to_column("rownames") %>%
  #   top_n(1, rownames) %>% pull(rownames)
  
  
  # Heatmap annotation
  anno_top <- HeatmapAnnotation(df = bar_col, col = cols,
                                show_annotation_name = show_annotation_name,
                                annotation_name_gp = gpar(fontsize = (fontsize-2)),
                                simple_anno_size = unit(3, "mm"),
                                annotation_legend_param = list(title_gp = gpar(fontsize = (fontsize-2), fontface = "bold"),
                                                               labels_gp = gpar(fontsize = (fontsize-2)),
                                                               grid_height = unit(3.5, "mm"),
                                                               grid_width = unit(3.5, "mm")))
  
  Heatmap(mdata,
          heatmap_legend_param = list(title = "Z-Score", 
                                      title_position = "leftcenter-rot",
                                      title_gp = gpar(fontsize = (fontsize-4), fontface = "bold"),
                                      labels_gp = gpar(fontsize = (fontsize-2)),
                                      grid_height = unit(fontsize-7, "mm"),
                                      grid_width = unit(1, "mm")),
          top_annotation = anno_top,
          right_annotation = right_annotation,
          col = c("violet", "black", "yellow"),
          row_title_gp = gpar(fontsize = fontsize),
          column_title_gp = gpar(fontsize = fontsize, fontface = "bold"),
          row_names_gp = gpar(fontsize = (fontsize-2)),
          column_names_gp = gpar(fontsize = (fontsize-2)),
          column_title = column_title,
          row_title = row_title,
          column_split =  bar_col[, 1],
          row_split = row_split,
          row_gap = unit(0.1, "pt"),
          cluster_rows = FALSE, cluster_columns = FALSE,
          show_column_names	= F,
          show_row_names	= F,
          use_raster = T,
          raster_device = "CairoPNG",
          raster_quality = res)
}




#### Other functions
vln <- function(object = mnnCT7, gene = gene, col = NULL, group.by = NULL){
  
  # single vln plot
  vln.single <- function(gene = NULL, object = NULL, col = NULL, group.by = NULL) {
    
    VlnPlot(object, gene, col, pt.size = 0, group.by = group.by) + NoLegend() + xlab(NULL) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      theme(axis.text.x = element_text(size = 10))
  }
  
  # Check input genes
  gene.check <- gene %in% rownames(object)
  if(0 %in% gene.check) message("Warning: ", paste0(gene[!gene.check], collapse = ","), " are not in features")
  gene <- gene[gene.check]
  if(S4Vectors::isEmpty(gene)) break
  
  # Ploting
  gene.list <- as.list(gene)
  vln.list <- lapply(gene.list, vln.single, object, col, group.by)
  grid.arrange(grobs = vln.list, ncol = 1)
}


# Permutation function
FDR.test <- function(v1, v2, nperm){
  cosi <- lsa::cosine(v1, v2)
  cosis <- sapply(1:nperm, function(x) lsa::cosine(v1, sample(v2, length(v2))))
  return(mean(cosis > as.numeric(cosi)))
}


# Box plots for gene expression
bplt <- function(gene = NULL, object = NULL,  cluster = "group", ncol = NULL, nrow = NULL, col = NULL, combine = T){
  
  # Single box plot function
  bplt.single <- function(gene = NULL, object = NULL,  cluster = "group", col = NULL) {
    
    crat_expr <- Biobase::exprs(object)[gene, ]
    crat_data <- as.data.frame(crat_expr)
    colnames(crat_data)[1] <- "org_value"
    crat_data <- cbind(crat_data, Biobase::pData(object)[cluster])
    
    if(is.null(col)) col <- scales::hue_pal()(Biobase::pData(bulk)[, cluster] %>% unique() %>% length())
    
    ggplot(data = crat_data, aes(x = !!ensym(cluster), y = org_value, fill = !!ensym(cluster))) +
      theme_classic() +
      geom_boxplot() +
      scale_fill_manual(values = col)+
      ylab("log2(intensity)") +
      xlab(cluster) +
      ggtitle(gene) + 
      theme(plot.title = element_text(hjust=0.5)) +
      NoLegend() 
  }
  
  # Check input genes
  gene.check <- gene %in% rownames(object)
  if(0 %in% gene.check) message("Warning: ", paste0(gene[!gene.check], collapse = ","), " are not in features")
  gene <- gene[gene.check]
  if(isEmpty(gene)) break
  
  # Ploting
  gene.list <- as.list(gene)
  plt.list <- lapply(gene.list, bplt.single, object = object,  cluster = cluster, col = col)
  
  if(combine == F) {return(plt.list)} else {
    if(is.null(ncol) && is.null(nrow)) nrow <- sqrt(length(gene))%/%1
    grid.arrange(grobs = plt.list, ncol = ncol, nrow = nrow)
  }
}




############################################
#### Batch effects correction using MNN ####
############################################
#### Read and preprocess 10x data #### 
data.list <- list()

for (i in 1:length(object.name)) {
  print(object.name[i])
  data.10x <- Read10X(paste("./count.v3/", file.name[i], sep = ""))
  colnames(data.10x) <- paste(colnames(data.10x), object.name[i], sep = ".")
  data.list[[i]] <- CreateSeuratObject(counts = data.10x, project = object.name[i], min.cells = 3, min.features = 200)
  data.list[[i]][["percent.mt"]] <- PercentageFeatureSet(object = data.list[[i]], pattern = "^MT-")
  data.list[[i]][["percent.rp"]] <- PercentageFeatureSet(object = data.list[[i]], pattern = "^RP")
  # data.list[[i]] <- subset(x = data.list[[i]], features = grep("^RP", rownames(data.list[[i]]), value = T, invert = T))
  data.list[[i]] <- subset(x = data.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt < 5)
  data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
  data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", 
                                         nfeatures = 5000, verbose = FALSE)
}

names(data.list) <- object.name
data.list




#### MNN Integratoin ####
# mnnCT7.raw.pc <- 50, mnnCT7.raw.dims <- 40, resolution = 0.2
mnnCT7.raw <- merge(data.list[[1]], data.list[-1])
mnnCT7.raw <- NormalizeData(mnnCT7.raw)
mnnCT7.raw <- FindVariableFeatures(mnnCT7.raw)


mnnCT7.raw.pc <- 50
mnnCT7.raw <- RunFastMNN(object.list = SplitObject(mnnCT7.raw, split.by = "orig.ident"), d = mnnCT7.raw.pc)

mnnCT7.raw.dims <- 40
mnnCT7.raw.res <- 0.2
mnnCT7.raw <- RunUMAP(mnnCT7.raw, reduction = "mnn", dims = 1:mnnCT7.raw.dims)
mnnCT7.raw <- FindNeighbors(mnnCT7.raw, reduction = "mnn", dims = 1:mnnCT7.raw.dims)
mnnCT7.raw <- FindClusters(mnnCT7.raw, resolution = mnnCT7.raw.res)
DimPlot(mnnCT7.raw, group.by = c("orig.ident", "ident"))
DimPlot(mnnCT7.raw, split.by = c("orig.ident"), group.by = "ident", label = T, ncol = 4, cols = col.CT7)
table(mnnCT7.raw$seurat_clusters)




#### Add tumor grade, clusters, and patients to meta data ####
# Add tumor grade infromation
mnnCT7.raw$histology <- str_split(mnnCT7.raw@meta.data$orig.ident, "_") %>% sapply("[[", 1)
mnnCT7.raw$cartilage <- ifelse(mnnCT7.raw@meta.data$orig.ident %in% c("Ben", "Low1", "Low2", "Med"),
                               "TURE", "FALSE")

mnnCT7.raw$patient <- factor(mnnCT7.raw$orig.ident, levels = patient)
mnnCT7.raw$cluster.id <- mnnCT7.raw$seurat_clusters
mnnCT7.raw$cluster <- factor(CT7.cluster[mnnCT7.raw$seurat_clusters], levels = CT7.cluster[c(0,5,4,1,3,7,2,6)+1])
DimPlot(mnnCT7.raw, group.by = c("patient",  "cluster"), cols = c(col.patient, col.CT7))




#### Finding differentially expressed features ####
DefaultAssay(mnnCT7.raw) <- "RNA"
mnnCT7.raw <- ScaleData(mnnCT7.raw, features = rownames(mnnCT7.raw))
mnnCT7.raw.markers <- FindAllMarkers(object = mnnCT7.raw, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mnnCT7.raw.markers <- subset(mnnCT7.raw.markers, p_val_adj < 0.05 & (pct.1 - pct.2) > 0)
table(mnnCT7.raw.markers$cluster)
# write.csv(mnnCT7.raw.markers, paste("./integration/mnnCT7.raw.markers.pc", mnnCT7.raw.pc, ".d", mnnCT7.raw.dims,
#                                     ".r", mnnCT7.raw.res, ".csv", sep = ""))

# Remove shared markers
mnnCT7.raw.unique.id <- which(table(mnnCT7.raw.markers$gene) == 1)
mnnCT7.raw.specific.markers <- subset(mnnCT7.raw.markers, gene %in% names(mnnCT7.raw.unique.id))
table(mnnCT7.raw.specific.markers$cluster)
# write.csv(mnnCT7.raw.specific.markers, paste("./integration/mnnCT7.raw.specific.markers.pc", mnnCT7.raw.pc, ".d",
#                                              mnnCT7.raw.dims, ".r", mnnCT7.raw.res, ".csv", sep = ""))


# Plot markers
mnnCT7.raw.top10 <- mnnCT7.raw.specific.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = mnnCT7.raw, features = as.character(mnnCT7.raw.top10$gene)) + NoLegend()




#### Enrichment analysis ####
mnnCT7.raw.go.list <- split.data.frame(mnnCT7.raw@misc$specific.markers, mnnCT7.raw@misc$specific.markers$cluster) %>% 
  lapply("[[", "gene") %>%
  lapply(gost)

# Extract GO results
mnnCT7.raw.result.list <- 
  lapply(mnnCT7.raw.go.list, function(go) {go$result %>% filter(term_size <1000) %>% dplyr::select(1:13)}) %>%
  lapply(function(x) {x$cluster <- (parent.frame()$i[]-1); return(x)})

mnnCT7.raw.go <- do.call(rbind, mnnCT7.raw.result.list)
table(mnnCT7.raw.go$source, mnnCT7.raw.go$cluster)




#### Infercnv ####
# prepare annotation file and count matrix
Leuk.marker <- FindMarkers(mnnCT7.raw, group.by = "cluster", ident.1 = "Leuk", min.pct = 0.1, logfc.threshold = 0.1, only.pos = T)
Leuk.gene <- Leuk.marker %>% tibble::rownames_to_column("gene") %>% filter(p_val_adj < 0.05)  %>% pull(gene)
length(Leuk.gene)

write.table(mnnCT7.raw[setdiff(rownames(mnnCT7.raw), Leuk.gene),]@meta.data["cluster"], "./infercnv/mnnCT7.raw_Leuk_anno.txt", 
            quote = F, row.names = T, col.names = F, sep = "\t")
write.table(GetAssayData(mnnCT7.raw[setdiff(rownames(mnnCT7.raw), Leuk.gene),], "counts"), './infercnv/mnnCT7.raw_Leuk_counts.matrix', quote = F, sep = "\t")


# read denoised observation and reference file
mnnCT7.raw_refe <- read.table("./infercnv/mnnCT7.rawcnv/infercnv.21_denoised.references.txt")
mnnCT7.raw_obse <- read.table("./infercnv/mnnCT7.rawcnv/infercnv.21_denoised.observations.txt")


# Calculate and save total cnv level
mnnCT7.raw_cnv <- cbind(mnnCT7.raw_obse, mnnCT7.raw_refe)
mnnCT7.raw_cnv <- mnnCT7.raw_cnv[rownames(mnnCT7.raw_cnv) %in% rownames(infercnv::genes), ]
mnnCT7.raw_cnv <- colSums((mnnCT7.raw_cnv-1)^2, na.rm = T)
write.table(mnnCT7.raw_cnv, "./infercnv/mnnCT7.raw_cnv.txt", quote = F, col.names = F)


# Read cnv_level
mnnCT7.raw_cnv <- read.delim("./infercnv/mnnCT7.raw_cnv.txt", sep = " ", header = F)
mnnCT7.raw_cnv <- setNames(mnnCT7.raw_cnv$V2, mnnCT7.raw_cnv$V1)


# Add cnv to meta data
mnnCT7.raw@meta.data$cnv <- mnnCT7.raw_cnv[colnames(mnnCT7.raw)]
VlnPlot(mnnCT7.raw, "cnv", pt.size = 0, split.by = "patient", group.by = "cluster", cols = col.patient)




#### AddModuleScore ####
# Marker gene score
gene.list <- split.data.frame(mnnCT7.raw.specific.markers, mnnCT7.raw.specific.markers$cluster)  %>% 
  lapply("[[", "gene")

mnnCT7.raw <- AddModuleScore(mnnCT7.raw, gene.list, assay = "RNA")

cluster.score <- paste0("cluster", levels(mnnCT7.raw), ".score")
colnames(mnnCT7.raw@meta.data)[grep("Cluster", colnames(mnnCT7.raw@meta.data))] <- cluster.score


cluster.score.plot.list <- mapply(function(feature, ident, col){
  (VlnPlot(mnnCT7.raw, features = feature, idents = ident, group.by = "orig.ident", pt.size = 0.1, cols = col) + 
     NoLegend()) %>% list()
}, cluster.score, levels(mnnCT7.raw), list(col.patient))

CombinePlots(cluster.score.plot.list, ncol = 4)



# Mean score matrix
score.cut.off <- 0.35
mean.score.matrix <- data.frame(patient = object.name)

for (i in 0:(length(levels(mnnCT7.raw))-1)) {
  scroe <- paste0("cluster", i, ".score")
  mean.score <- mnnCT7.raw@meta.data %>%
    # tibble::rownames_to_column("cell.id") %>%
    group_by(patient) %>%
    filter(seurat_clusters == i) %>%
    summarise(mean(!!ensym(scroe)))
  colnames(mean.score)[2] <- i
  mean.score.matrix <- right_join(mean.score, mean.score.matrix)
}
mean.score.matrix


# Average cluster score 
mean.score.melt <- reshape::melt(mean.score.matrix %>% as.data.frame())
colnames(mean.score.melt)[3] <- "mean"
mean.score.melt$mean.cut.off <- mean.score.melt$mean > score.cut.off


# Distance to highest scroe
mean.score.dist.matrix <- t(rowMaxs(t(mean.score.matrix[-1]), na.rm = T) - t(mean.score.matrix[-1])) %>%
  cbind(mean.score.matrix[1])


mean.score.dist.melt <- reshape::melt(mean.score.dist.matrix %>% as.data.frame())
colnames(mean.score.dist.melt)[3] <- "dist"
mean.score.dist.melt$dist.cut.off <- mean.score.dist.melt$dist < score.cut.off

mean.score.melt <- left_join(mean.score.melt, mean.score.dist.melt)
mean.score.melt$denosed <- mean.score.melt$mean.cut.off & mean.score.melt$dist.cut.off
mean.score.melt$mark <- ifelse(mean.score.melt$denosed == T, "*", "")
mean.score.melt$patient <- factor(mean.score.melt$patient, levels = object.name)
mean.score.melt$variable <- factor(mean.score.melt$variable, levels = levels(mnnCT7.raw))

mean.score.melt$cluster <- CT7.cluster[mean.score.melt$variable]
mean.score.melt$cluster <- factor(mean.score.melt$cluster, levels = CT7.cluster[c(0,5,4,1,3,7,2,6)+1])
# mnnCT7.raw@misc$mean.score.melt <- mean.score.melt

ggplot(mean.score.melt, aes(x = cluster, y = mean, fill = patient)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept  = score.cut.off, col = "red") +
  theme_classic() +
  xlab("cluster") + ylab("average cluster score") + labs(fill = "patient")+
  geom_text(aes(label = mark), position=position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_manual(values = col.patient)



ggplot(mean.score.melt, aes(x = variable, y = dist, fill = patient)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept  = score.cut.off, col = "red") +
  theme_classic() +
  xlab("cluster") + ylab("distance to highest scroe") + labs(fill = "patient") +
  geom_text(aes(label = mark), position=position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_manual(values = col.patient)




#### Denoising ####
# Include all patient of cluster 3
mean.score.melt$denosed.2 <- mean.score.melt$denosed
mean.score.melt[which(mean.score.melt$variable == 2), "denosed.2"] <- TRUE
mean.score.denosed <- subset(mean.score.melt, denosed.2 == T) %>% mutate(test = paste(patient, variable, sep = "-"))
cell.denosed <- paste(mnnCT7.raw$orig.ident, mnnCT7.raw$seurat_clusters, sep = "-") %in% mean.score.denosed$test
table(cell.denosed)



# Sebset mnnCHS7.raw with clear cells
mnnCT7 <- mnnCT7.raw[, cell.denosed]


# corrected cluster score plot
cluster.score.plot.list <- mapply(function(feature, ident, col){
  (VlnPlot(mnnCT7, features = feature, idents = ident, group.by = "orig.ident", pt.size = 0.1, cols = col) + 
     NoLegend()) %>% list()
}, cluster.score, levels(mnnCT7), list(col.patient))

CombinePlots(cluster.score.plot.list, ncol = 4)

# Name cluster
mnnCT7$cluster.id <- mnnCT7$seurat_clusters
names(CT7.cluster) <- levels(mnnCT7)
mnnCT7$cluster <- factor(recode(mnnCT7$cluster.id, !!!CT7.cluster), levels = CT7.cluster[c(0,5,4,1,3,7,2,6)+1])
mnnCT7$subcluster <- factor(mnnCT7$subcluster, levels = c(c("High1a", "High1b", "High2a", "High2b", "nStro", "Stro1"), 
                                                          CT7.cluster[c(0,5,4,1,3,7,2,6)+1]))



#########################
#### Denoised mnnCT7 ####
#########################
#### Finding differentially expressed features ####
# Finding differentially expressed features (cluster biomarkers)
DefaultAssay(mnnCT7) <- "RNA"
mnnCT7 <- ScaleData(mnnCT7, features = rownames(mnnCT7))
mnnCT7.markers <- FindAllMarkers(object = mnnCT7, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mnnCT7.markers <- subset(mnnCT7.markers, p_val_adj < 0.05 & (pct.1 - pct.2) > 0)
table(mnnCT7.markers$cluster)


# Remove shared markers
mnnCT7.unique.id <- which(table(mnnCT7.markers$gene) == 1)
mnnCT7.specific.markers <- subset(mnnCT7.markers, gene %in% names(mnnCT7.unique.id))
table(mnnCT7.specific.markers$cluster)



# Plot markers
mnnCT7.top10 <- mnnCT7.specific.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = mnnCT7, features = as.character(mnnCT7.top10$gene)) + NoLegend()





#### UMAP and Heatmap ####
mnnCT7@meta.data$orig.ident <- factor(mnnCT7@meta.data$orig.ident, levels = patients)
DimPlot(mnnCT7, group.by = c("cluster",  "patient"), cols = c(col.CT7, col.patient))


DimPlot(mnnCT7, split.by = c("orig.ident"), group.by = "ident", label = T, cols = col.CT7, ncol = 4) +
  NoLegend()

ggplot(mnnCT7@meta.data, aes(x = patient, fill = cluster)) +
  geom_bar(stat = "count", position = "fill") +
  scale_fill_manual(values = col.CT7) +
  xlab("patient") + ylab("percentage") + labs(fill = "cluster") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))





#### Enrichment analysis ####
mnnCT7.go.list <- split.data.frame(mnnCT7@misc$specific.markers, mnnCT7@misc$specific.markers$cluster) %>% 
  lapply("[[", "gene") %>%
  lapply(gost)

# Extract GO results
mnnCT7.result.list <- 
  lapply(mnnCT7.go.list, function(go) {go$result %>% filter(term_size <1000) %>% dplyr::select(1:13)}) %>%
  lapply(function(x) {x$cluster <- (parent.frame()$i[]-1); return(x)})

mnnCT7.go <- do.call(rbind, mnnCT7.result.list)
table(mnnCT7.go$source, mnnCT7.go$cluster)




#### Cell cycle analysis ####
exp.mat <- read.table(file = "./cell_cycle/nestorawa_forcellcycle_expressionMatrix.txt", 
                      header = TRUE, as.is = TRUE, row.names = 1)
cc.genes <- readLines(con = "./cell_cycle/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

mnnCT7 <- CellCycleScoring(object = mnnCT7, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
levels(mnnCT7$Phase) <- c("G1", "S", "G2M")

DimPlot(mnnCT7, reduction = "umap", group.by = "Phase")

# Plot cell cycle and clusters
ggplot(mnnCT7@meta.data, aes(x = seurat_clusters, fill = Phase)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  xlab("cluster") + ylab("percentage") +
  # facet_grid(. ~ orig.ident) +
  theme_classic() +
  NoLegend()

ggplot(mnnCT7@meta.data, aes(x = orig.ident, fill = Phase)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  xlab("patient") + ylab("percentage") +
  # facet_grid(. ~ orig.ident) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 





#### Infercnv filtering marker genes ####
# Finding differentially expressed features (cluster biomarkers)
DefaultAssay(mnnCT7) <- "RNA"
mnnCT7 <- ScaleData(mnnCT7, features = rownames(mnnCT7))
mnnCT7.markers <- FindAllMarkers(object = mnnCT7, only.pos = T, min.pct = 0.1, logfc.threshold = 0.1)
mnnCT7.markers <- subset(mnnCT7.markers, p_val_adj < 0.05 & (pct.1 - pct.2) > 0)
table(mnnCT7.markers$cluster)

# write.csv(mnnCT7.markers, "integration/mnnCT7.markers.min.pct.1.logfc.1.csv")
mnnCT7.markers <- read.csv("integration/mnnCT7.markers.min.pct.1.logfc.1.csv", header = T, row.names = 1, stringsAsFactors = F)
mnnCT7.markers$cluster.name <- CT7.cluster[mnnCT7.markers$cluster %>% as.character()]
length(unique(mnnCT7.markers$gene))


markers <- setdiff(rownames(mnnCT7), mnnCT7@misc$specific.markers$gene)
length(markers)

# markers <- setdiff(rownames(mnnCT7), mnnCT7.markers$gene[mnnCT7.markers$cluster.name == "Leuk"])
# length(markers)

write.table(mnnCT7[markers,]@meta.data["cluster"], "./infercnv/mnnCT7_specific.marker_anno.txt", 
            quote = F, row.names = T, col.names = F, sep = "\t")
write.table(GetAssayData(mnnCT7[markers,], "counts"), './infercnv/mnnCT7_specific.marker_counts.matrix', quote = F, sep = "\t")




#### Infercnv line plot ####
# read denoised observation and reference file
mnnCT7_refe <- read.table("./infercnv/mnnCT7_Leuk_cnv/infercnv.21_denoised.references.txt")
mnnCT7_obse <- read.table("./infercnv/mnnCT7_Leuk_cnv/infercnv.21_denoised.observations.txt")

cnv <- cbind(mnnCT7_obse, mnnCT7_refe)


cnv.genes <- read.table("infercnv/gene.coor.csv", header = T, row.names = 1, sep = ",", stringsAsFactors = F)
rownames(cnv.genes) <- cnv.genes$gene_name

# Subset cnv
gene <- cnv.genes[rownames(cnv.genes) %in% rownames(cnv), ]


mnnCT7.cnv <- cnv[rownames(gene), colnames(mnnCT7)] %>% t()
mnnCT7.cnv[1:5, 1:5]
dim(mnnCT7.cnv)


# Subset cells according to patients and clusters (cluster with more than 20 cells)
cell.list <- split.data.frame(mnnCT7@meta.data, mnnCT7$patient) %>%
  lapply(function(x) split(x$cell.id, x$cnvcluster))


cell.list <- lapply(cell.list, function(x) x[sapply(x, length) > 20])


# Plot chondoid tumors
# Add gaps_col
gaps_col <- gene %>% group_by(gene[, "rank"]) %>% tibble::rowid_to_column("rownames") %>% 
  top_n(1, rownames) %>% pull(rownames)

gaps_col[c(5,10,18,20)] <- gaps_col[c(5,10,18,20)] + 1
gaps_col[1:10] <- gaps_col[1:10] - 1

gene$chr <- ""
gene$chr[gaps_col - table(gene$rank)%/%2] <- 1:22


# Add highlighted genes
gene$col <- gene$rank

# gene$col[which(gene$gene_name == "COL10A1"):(which(gene$gene_name == "COL10A1") + 10)] <- "anno"
# 
# gene$chr[which(gene$gene_name == "COL10A1")] <- "COL10A1"


# Mean CNV level of each regin within cell list
patient.cnv <- lapply(cell.list[1:4], function(patient) {
  sapply(patient, function(x) rowMeans(cnv[,x, drop = F]))
})


# chon.cnv <- patient.cnv[1:4]


# Create plot list
cnv.list <- lapply(patient.cnv, function(x) {
  # Melt the data for ploting
  join.mtx <- cbind(gene, x[gene$gene_name,])
  join.melt <- reshape::melt(join.mtx, id.vars = colnames(gene))
  join.melt$gene_name <- factor(join.melt$gene_name, levels = gene$gene_name)
  
  # Plot data
  ggplot(join.melt) +
    geom_line(aes(x = gene_name, y = value, colour = variable, group = variable)) +
    scale_color_manual(values = col.allcnv) +
    xlab("chromosome") + ylab(names(patient.cnv)[parent.frame()$i[]]) +
    theme(axis.ticks.x = element_line(colour = col.chr[join.melt$col]), legend.title = element_blank(),
          axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 7), 
          axis.title = element_text(size = 12), axis.line.y = element_line()) +
    scale_x_discrete(labels= join.melt$chr, guide = guide_axis(n.dodge=2)) +
    geom_hline(yintercept = c(0.9, 1.1), col = "gray", linetype = "dashed") +
    ylim(min(mnnCT7.cnv), max(mnnCT7.cnv)) +
    NoLegend()
})


# Add annotation box for ch6
cnv.list[2:4] <- lapply(cnv.list[2:4], function(x) x + 
                          geom_rect(aes(xmin= "SENP6", xmax= "PDCD2", 
                                        ymin= min(mnnCT7.cnv), ymax= max(mnnCT7.cnv)), 
                                    color = "gray50", alpha = 0, linetype = "longdash")
)


# Add annotation box for ch11 and Ch19 of Low_L28
gene.11 <- range(which(gene$rank == 11))
cnv.list[3] <- lapply(cnv.list[3], function(x) x + 
                        geom_rect(aes(xmin= gene.11[1], xmax= gene.11[2], 
                                      ymin= min(mnnCT7.cnv), ymax= max(mnnCT7.cnv)), 
                                  color = "blue", alpha = 0, linetype = "longdash")
)


gene.19 <- range(which(gene$rank == 19))
cnv.list[3] <- lapply(cnv.list[3], function(x) x + 
                        geom_rect(aes(xmin= gene.19[1], xmax= gene.19[2], 
                                      ymin= min(mnnCT7.cnv), ymax= max(mnnCT7.cnv)), 
                                  color = "red", alpha = 0, linetype = "longdash")
)


# Remove x axis annotation
cnv.list[1:3] <- lapply(cnv.list[1:3], function(x) {
  x + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
})


# Merge Chon plots
cnv.plot.null <- plot_grid(plotlist = cnv.list[1:4], ncol = 1, rel_heights = c(1,1,1,1.6), greedy = F)
cnv.plot.ylab <- textGrob("mean CNV acrosss cluster", rot = 90)
cnv.plot.L <- (ggplot(mnnCT7@meta.data[mnnCT7$histology %in% c("Ben", "Low", "Med"),]) +
                 geom_line(aes(x= patient, y = 1, group = subcluster, color = subcluster)) + theme_classic() +
                 scale_color_manual(values = col.all, breaks = c("Chon1", "Chon2", "Stro", "Leuk")) +
                 labs(col = "cluster")) %>% get_legend()

cnv.plot <- plot_grid(cnv.plot.ylab, cnv.plot.null, cnv.plot.L, nrow = 1, rel_widths = c(2,50,7))
# saveRDS(cnv.plot, "figure/chon.cnv.plot_thesis.rds")


# Save CNV gene at chromosome 6
# write.csv(cnv.genes[which(cnv.genes$gene_name == "SENP6"):which(cnv.genes$gene_name == "PDCD2"), ], 
#           "infercnv/cnv_gene/cnv_gene_ch6.csv")




# Plot high grade tumors
# Add gaps_col
gaps_col <- gene %>% group_by(gene[, "rank"]) %>% tibble::rowid_to_column("rownames") %>% 
  top_n(1, rownames) %>% pull(rownames)

gaps_col[c(5,10,18,20)] <- gaps_col[c(5,10,18,20)] + 1

gene$chr <- ""
gene$chr[gaps_col - table(gene$rank)%/%2] <- 1:22


# Add highlighted genes
gene$col <- gene$rank

gene$col[which(gene$gene_name == "TGFB2"):(which(gene$gene_name == "TGFB2") + 4)] <- "anno"

gene$chr[which(gene$gene_name == "TGFB2") + 1] <- "TGFB2"


# Mean CNV level of each regin within cell list
patient.cnv <- lapply(cell.list[5:7], function(patient) {
  sapply(patient, function(x) rowMeans(cnv[,x, drop = F]))
})


# Create plot list
cnv.list <- lapply(patient.cnv, function(x) {
  # Melt the data for ploting
  join.mtx <- cbind(gene, x[gene$gene_name,])
  join.melt <- reshape::melt(join.mtx, id.vars = colnames(gene))
  join.melt$gene_name <- factor(join.melt$gene_name, levels = gene$gene_name)
  
  # Plot data
  ggplot(join.melt) +
    geom_line(aes(x = gene_name, y = value, colour = variable, group = variable)) +
    scale_color_manual(values = col.allcnv) +
    xlab("chromosome") + ylab(names(patient.cnv)[parent.frame()$i[]]) +
    theme(axis.ticks.x = element_line(colour = col.chr[join.melt$col]), legend.title = element_blank(),
          axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 7), 
          axis.title = element_text(size = 10), axis.line.y = element_line()) +
    scale_x_discrete(labels= join.melt$chr, guide = guide_axis(n.dodge=2)) +
    geom_hline(yintercept = c(0.9, 1.1), col = "gray", linetype = "dashed") +
    ylim(min(mnnCT7.cnv), max(mnnCT7.cnv)) +
    NoLegend()
})


# Add annotation box for L31
cnv.list[[1]] <- cnv.list[[1]] +
  geom_rect(aes(xmin= "NBPF15", xmax= "ZNF672",
                ymin= min(mnnCT7.cnv), ymax= max(mnnCT7.cnv)),
            color = "dodgerblue", alpha = 0)


# Remove x axis annotation
cnv.list[1:2] <- lapply(cnv.list[1:2], function(x) {
  x + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
})


# Merge Chon plots
cnv.plot.null <- plot_grid(plotlist = cnv.list, ncol = 1, rel_heights = c(1,1,1.3), greedy = F)
cnv.plot.ylab <- textGrob("mean CNV acrosss cluster", rot = 90)
cnv.plot.L <- (ggplot(mnnCT7@meta.data) +
                 geom_line(aes(x= patient, y = 1, group = cnvcluster, color = cnvcluster)) + theme_classic() +
                 scale_color_manual(values = col.allcnv, breaks = c("Undi", "Hype", "Myxo", "nStro", "Stro1", "Leuk")) +
                 labs(col = "cluster")) %>% get_legend()

cnv.plot <- plot_grid(cnv.plot.ylab, cnv.plot.null, cnv.plot.L, nrow = 1, rel_widths = c(2,50,5))
# saveRDS(cnv.plot, "figure/high.cnv.plot.rds")




#### Analysis on cnv of specific chromosomes of Low_L28 ####
# read denoised observation and reference file
mnnCT7_refe <- read.table("./infercnv/mnnCT7_Leuk_cnv/infercnv.21_denoised.references.txt")
mnnCT7_obse <- read.table("./infercnv/mnnCT7_Leuk_cnv/infercnv.21_denoised.observations.txt")

cnv <- cbind(mnnCT7_obse, mnnCT7_refe)
cnv[1:5, 1:5]

cnv.genes <- read.table("infercnv/gene.coor.csv", header = T, row.names = 1, sep = ",", stringsAsFactors = F)
rownames(cnv.genes) <- cnv.genes$gene_name

# Subset cnv
gene <- cnv.genes[rownames(cnv.genes) %in% rownames(cnv), ]
head(gene)


# Chromosome 19
cells <- mnnCT7@meta.data %>% filter(patient == "Low_L28" & cluster %in% c("Chon1", "Chon2")) %>% pull(cell.id)

L28.19 <- cnv[gene$chr == "19", cells]
dim(L28.19)

saveRDS(L28.19, "infercnv/heatmap/L28.19.rds")


ha <- HeatmapAnnotation(cluster = mnnCT7@meta.data[cells, "cluster"], 
                        which = "column",
                        col = list(cluster = col.all[c("Chon1", "Chon2")]))

L28.19.heatmap <- Heatmap(L28.19, 
                          name = "CNV",
                          col = circlize::colorRamp2(c(1.2,1,0.8), c("red", "white", "blue")),
                          top_annotation = ha, 
                          cluster_rows = F, 
                          show_column_names = F, 
                          show_row_names = F, 
                          column_split = 4,
                          clustering_method_columns = "ward.D2")



L28.19.hc <- hclust(dist(t(L28.19)), method = "ward.D2")

L28.19.cluster <- cutree(L28.19.hc, k = 2)
table(L28.19.cluster)

mnnCT7$L28.19.cluster <- NA

mnnCT7@meta.data[names(L28.19.cluster), "L28.19.cluster"] <- L28.19.cluster

DimPlot(mnnCT7[, cells], group.by = "L28.19.cluster")


L28.19.cnv <- L28.19 -1
table(L28.19.cnv<0)
L28.19.cnv[L28.19.cnv<0] <- 0
table(L28.19.cnv<0)

L28.19.cnv <- abs(colSums(L28.19.cnv))
mnnCT7$L28.19.cluster <- L28.19.cnv
mnnCT7@meta.data[names(L28.19.cnv), "L28.19.cnv"] <- L28.19.cnv
FeaturePlot(mnnCT7[, cells], features = "L28.19.cnv", cols = c("whitesmoke", "red"), max.cutoff = 70, min.cutoff = 40)



# Chromosome 11
cells <- mnnCT7@meta.data %>% filter(patient == "Low_L28" & cluster %in% c("Chon1", "Chon2")) %>% pull(cell.id)

L28.11 <- cnv[gene$chr == "11", cells]
dim(L28.11)

saveRDS(L28.11, "infercnv/heatmap/L28.11.rds")


ha <- HeatmapAnnotation(cluster = mnnCT7@meta.data[cells, "cluster"], 
                        which = "column",
                        col = list(cluster = col.all[c("Chon1", "Chon2")]))

L28.11.heatmap <- Heatmap(L28.11, 
                          name = "CNV",
                          col = circlize::colorRamp2(c(1.2,1,0.8), c("red", "white", "blue")),
                          top_annotation = ha, 
                          cluster_rows = F, 
                          show_column_names = F, 
                          show_row_names = F, 
                          column_split = 2,
                          clustering_method_columns = "ward.D2")



L28.11.hc <- hclust(dist(t(L28.11)), method = "ward.D2")

L28.11.cluster <- cutree(L28.11.hc, k = 2)
table(L28.11.cluster)

mnnCT7$L28.11.cluster <- NA

mnnCT7@meta.data[names(L28.11.cluster), "L28.11.cluster"] <- L28.11.cluster

DimPlot(mnnCT7[, cells], group.by = "L28.11.cluster")


L28.11.cnv <- L28.11 -1
table(L28.11.cnv>0)
L28.11.cnv[L28.11.cnv>0] <- 0
table(L28.11.cnv>0)

L28.11.cnv <- abs(colSums(L28.11.cnv))
mnnCT7$L28.11.cluster <- L28.11.cnv
mnnCT7@meta.data[names(L28.11.cnv), "L28.11.cnv"] <- L28.11.cnv
FeaturePlot(mnnCT7[, cells], features = "L28.11.cnv", cols = c("whitesmoke", "blue"), max.cutoff = 12, min.cutoff = 4)




#### Chon2 cnv gain region at chromosome 11 of Low_2 ####
# read denoised observation and reference file
mnnCT7_refe <- read.table("./infercnv/mnnCT7_Leuk_cnv/infercnv.21_denoised.references.txt")
mnnCT7_obse <- read.table("./infercnv/mnnCT7_Leuk_cnv/infercnv.21_denoised.observations.txt")

cnv <- cbind(mnnCT7_obse, mnnCT7_refe)
cnv[1:5, 1:5]

cnv.genes <- read.table("infercnv/gene.coor.csv", header = T, row.names = 1, sep = ",", stringsAsFactors = F)
rownames(cnv.genes) <- cnv.genes$gene_name

# Subset cnv
gene <- cnv.genes[rownames(cnv.genes) %in% rownames(cnv), ]
head(gene)


# Chromosome 19
mnnCT7$cell.id <- colnames(mnnCT7)
cells <- mnnCT7@meta.data %>% filter(patient == "Low_2" & cluster %in% c("Chon1", "Chon2")) %>% pull(cell.id)

L28.19 <- cnv[gene$chr == "19", cells]
dim(L28.19)

# saveRDS(L28.19, "infercnv/heatmap/L28.19.rds")


Chon1.cells <- mnnCT7@meta.data %>% filter(patient == "Low_2" & cluster %in% c("Chon1")) %>% pull(cell.id)
Chon2.cells <- mnnCT7@meta.data %>% filter(patient == "Low_2" & cluster %in% c("Chon2")) %>% pull(cell.id)

L28.19.mean <- data.frame(Chon1 = rowMeans(L28.19[, Chon1.cells]), Chon2 = rowMeans(L28.19[, Chon2.cells]))


table(L28.19.mean$Chon2 - L28.19.mean$Chon1 > 0.1)

cnv.genes[rownames(L28.19)[L28.19.mean$Chon2 - L28.19.mean$Chon1 > 0.1],] %>% 
  write.csv("infercnv/heatmap/L28.19_CNV_gain_Chon2.csv")



#### Analysis on the Ch6.q ####
# Copy number loss gene
ch6.cnv.genes <- cnv.genes[which(cnv.genes$gene_name == "SENP6"):which(cnv.genes$gene_name == "PDCD2"), "gene_name"]


# Apoptosis gene
apo.genes <- gconvert("GO:0043065")$name # positive regulation of apoptotic process

ggvenn(list("Ch6.CNV       " = ch6.cnv.genes, "  Apoptosis" = apo.genes), show_percentage = F, set_name_size = 5)
intersect(ch6.cnv.genes, apo.genes)

# Load TF
tf <- read.table("tf/hgTables_DDIT3_1000.txt")[,c(1:5,10)]
head(tf)

colnames(tf) <- c("seqname", "source", "feature", "start", "end", "SYMBOL")
head(tf)

ggvenn(list("Ch6.CNV" = ch6.cnv.genes, "TF.CHOP" = tf$SYMBOL), show_percentage = F, set_name_size = 5)
intersect(ch6.cnv.genes, tf$SYMBOL)



ch6.apo.gene <- intersect(ch6.cnv.genes, apo.genes) 
bulk.ch6.apo.gene <- intersect(rownames(bulk), ch6.apo.gene)

F3b2.list <- list()

for (gene in bulk.ch6.apo.gene) {
  
  bulk.level <- bulk[, bulk$OS.delay != "#N/A"]
  
  level <- Biobase::exprs(bulk.level)[gene,]
  
  bulk.level$level <- ifelse(level < quantile(level, 0.25), "Low", 
                             ifelse(level > quantile(level, 0.75), "High", "Medium"))
  table(bulk.level$level)
  
  bulk.level <- bulk.level[, bulk.level$level != "Medium"]
  # bulk.level$level <- factor(bulk.level$level, levels = c("High", "low"))
  
  
  meta.data <- Biobase::pData(bulk.level)
  meta.data$OS.delay <- as.numeric(as.character(meta.data$OS.delay))
  meta.data$OS.event <- as.numeric(as.character(meta.data$OS.event))
  meta.data <- arrange(meta.data, level)
  
  
  # OS survival
  for (i in 1:nrow(meta.data)) {
    if(meta.data$OS.delay[i] >= 60) {
      meta.data$OS.delay[i] <- 60
      meta.data$OS.event[i] <- 0
    }
  }
  
  
  sd <- survdiff(Surv(OS.delay, OS.event) ~ level, meta.data)
  pvalue <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  if(pvalue >= 0.05) next
  pvalue <- format(pvalue, digits = 1, nsmall = 3)
  # 
  # col.sur <- c(col.group[-5], alpha(col.group[5], 0.3))
  # names(col.sur) <- paste0(levels(meta.data$level), " (n=", table(meta.data$level), ")")
  
  sur.plot <- ggsurvplot(
    fit = survfit(Surv(OS.delay, OS.event) ~ level, meta.data), 
    # pval = "\n         p value = 3e-10",
    xlab = "month",
    ylab = "overall survival rate",
    legend = "right",
    legend.title = "expression",
    # legend.labs = meta.data$level %>% unique(),
    # palette = col.sur,
    # linetype = c(1,1,1,1,5),
    ggtheme = theme_survminer(base_size = 10),
    font.legend = c(10, "plain", "black"),
    font.tickslab = c(10, "plain", "black"),
    font.x = c(12, "plain", "black"),
    font.y = c(12, "plain", "black"))
  
  F3b2.list[[gene]] <- sur.plot$plot + 
    # geom_vline(xintercept = 60, col = "red") +
    annotate(geom = "text", label = paste0("p value = ", pvalue), x = 32, y = 0.1, size = 4) +
    scale_y_continuous(labels = scales::percent_format(), limits=c(0,1)) +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size = 9), 
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) + 
    labs(title = gene)  
}




#### DEG analysis on Chon1 and Chon2 ####
Chon2.Chon1 <- FindMarkers(mnnCT7, group.by = "cluster", ident.1 = "Chon2", ident.2 = "Chon1",
                           min.pct = 0, logfc.threshold = 0)

# Tidy marker genes table
Chon2.Chon1 <- Chon2.Chon1 %>% 
  tibble::rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05) %>%
  mutate(cluster = ifelse(avg_logFC > 0, "Chon2", "Chon1")) %>%
  group_by(cluster) %>%
  arrange(desc(abs(avg_logFC)), .by_group = T)
rownames(Chon2.Chon1) <- Chon2.Chon1$gene
table(Chon2.Chon1$cluster)

write.csv(Chon2.Chon1, "integration/Chon2.Chon1.markers.csv")

Chon2.Chon1 <- read.csv("integration/Chon2.Chon1.markers.csv", row.names = 1)


# Vocanol plot
features <-  intersect(gconvert(c("GO:0002062", "GO:0034976"))$name, Chon2.Chon1$gene) %>% as.character()
Chon2.Chon1$label <- Chon2.Chon1$gene
Chon2.Chon1$label[!(Chon2.Chon1$gene %in% features)] <- NA
Chon2.Chon1$adj_logFC <- pmin(pmax(Chon2.Chon1$avg_logFC, -2), 2)


ggplot(Chon2.Chon1, aes(x = adj_logFC, y = -log10(p_val_adj))) +
  geom_point(aes(col = cluster)) +
  geom_vline(xintercept = c(0.25, -0.25), col = "grey50", linetype = "dashed") +
  geom_hline(yintercept = 0.05, col = "grey50", linetype = "dashed") +
  geom_text_repel(aes(label = label), size = 4,
                  box.padding   = 0.35,
                  point.padding = 0.5,
                  segment.color = 'grey50') +
  scale_x_continuous(limits = c(-2, 2)) +
  scale_color_manual(values = col.all) +
  labs(x = "log2 fold change", y = "-log10(p value)") +
  theme_classic()



# CHOP regulators
Chon2.Chon1.filtered <- read.csv("integration/Chon2.Chon1.markers.filtered.csv", row.names = 1)

Chon2.marker <- Chon2.Chon1.filtered$gene[Chon2.Chon1.filtered$cluster == "Chon2"]

# Load TF
tf <- read.table("tf/hgTables_DDIT3_1000.txt")[,c(1:5,10)]
head(tf)

colnames(tf) <- c("seqname", "source", "feature", "start", "end", "SYMBOL")
head(tf)


ggvenn(list(CHOP.TF = tf$SYMBOL, Chon1.marker = Chon2.marker))

intersect(Chon2.marker, tf$SYMBOL)



#####################
####      OA     ####
#####################
#### Load data, integration, and clustering ####
# Read data
OA.raw_counts <- read.table("OA/GSE104782_allcells_UMI_count.txt.gz", header = T, row.names = 1)
OA.meta.data <- read.csv("OA/GSE104782_Table_Cell_quality_information_and_clustering_information.csv",
                         header = T, row.names = 1)


# Create Seurat object
OA <- CreateSeuratObject(OA.raw_counts, min.cells = 3, min.features = 200, project = "OA", meta.data = OA.meta.data)


# Add meta data
OA@meta.data$histology <- str_split(colnames(OA), "[_.]+") %>% sapply("[[", 2)


# QC and selecting cells for further analysis
OA[["percent.mt"]] <- PercentageFeatureSet(object = OA, pattern = "^MT")

# Visualize QC metrics as a violin plot
VlnPlot(object = OA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(object = OA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = OA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# OA <- subset(x = OA, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt < 5)

# Normalizing the data
OA <- NormalizeData(object = OA, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection)
OA <- FindVariableFeatures(object = OA, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
OA.top10 <- head(x = VariableFeatures(object = OA), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = OA)
LabelPoints(plot = plot1, points = OA.top10, repel = TRUE)


# Scaling the data
OA.all.genes <- rownames(x = OA)
OA <- ScaleData(object = OA, features = OA.all.genes)
#OA <- ScaleData(object = OA, features = OA.all.genes, vars.to.regress = "percent.mt")

# Perform linear dimensional reduction
OA <- RunPCA(object = OA, features = VariableFeatures(object = OA))

# Examine and visualize PCA results a few different ways
print(x = OA[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object = OA, dims = 1:2, reduction = "pca")
DimPlot(object = OA, reduction = "pca")
DimHeatmap(object = OA, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the 'dimensionality' of the dataset
# OA <- JackStraw(object = OA, num.replicate = 100)
# OA <- ScoreJackStraw(object = OA, dims = 1:20)
# JackStrawPlot(object = OA, dims = 1:15)

ElbowPlot(object = OA, ndims = 30)


# Cluster the cells
OA.dims <- 5
OA.res <- 0.1
OA <- FindNeighbors(object = OA, dims = 1:OA.dims)
OA <- FindClusters(object = OA, resolution = OA.res)

# Run non-linear dimensional reduction (UMAP/tSNE)
OA <- RunUMAP(object = OA, dims = 1:OA.dims)
DimPlot(object = OA, reduction = "umap")
table(OA@meta.data$seurat_clusters)
DimPlot(object = OA, reduction = "umap", group.by = c("seurat_clusters", "Cluster", "histology"))




#### Gene markers and gene set enrichemnt ####
# Finding differentially expressed features (cluster biomarkers)
Idents(OA) <- OA$Cluster
OA.markers <- FindAllMarkers(object = OA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
table(OA.markers$cluster)
OA.specific.markers <- subset(OA.markers, p_val_adj < 0.05 & (pct.1 - pct.2) > 0)
OA.unique.id <- which(table(OA.specific.markers$gene) == 1)
OA.specific.markers <- subset(OA.specific.markers, gene %in% names(OA.unique.id))
table(OA.specific.markers$cluster)

OA@misc$specific.markers <- OA.specific.markers



# Enrichment analysis
OA.go.list <- split(OA.specific.markers$gene, OA.specific.markers$cluster)[-1] %>% 
  lapply(gost)


# Extract GO results with term_size < 1000
OA.result.list <- 
  lapply(OA.go.list, function(go) {go$result %>% filter(term_size <1000) %>% dplyr::select(1:13)}) %>%
  lapply(function(x) {x$cluster <- (parent.frame()$i[]-1); return(x)})

OA.go <- do.call(rbind, OA.result.list)
table(OA.go$source, OA.go$cluster)


# Vlnplot of selected markers
DimPlot(OA, reduction = "umap")
VlnPlot(OA, c("SLC2A1", "MT-ND3", "TOP2A", "ACTA2", "CD68"), pt.size = 0, ncol = 4)




#### Similarity analysis ####
# Calculate similarity matrix
rib.cluster <- unique(OA$Cluster)[-1]
chon.cluster <- intersect(levels(mnnCT7$cluster), mnnCT7$cluster)

nperm <- 10000
simatrix <- matrix(nrow = length(rib.cluster), ncol = length(chon.cluster), dimnames = list(rib.cluster, chon.cluster))
FDRmatrix <- simatrix

for (i in rib.cluster) {
  shared.gene <- intersect(rownames(mnnCT7), OA@misc$specific.markers$gene[OA@misc$specific.markers$cluster == i])
  rib.ctl <- OA$Cluster != i
  rib.dif.vec <- rowMeans(OA[shared.gene,!rib.ctl]@assays$RNA@scale.data) - rowMeans(OA[shared.gene,rib.ctl]@assays$RNA@scale.data)
  for (j in chon.cluster) {
    chon.ctl <- mnnCT7$cluster != j
    chon.dif.vec <- rowMeans(mnnCT7[shared.gene,!chon.ctl]@assays$RNA@scale.data) - rowMeans(mnnCT7[shared.gene,chon.ctl]@assays$RNA@scale.data)
    simatrix[i,j] <- lsa::cosine(rib.dif.vec, chon.dif.vec)
    FDRmatrix[i,j] <- FDR.test(chon.dif.vec, rib.dif.vec, nperm)
  }
}

FDRmatrix

# Radar plot
data <- simatrix %>% as.data.frame()
data[data<0] <- 0
grid.min <- min(data)
grid.max <- max(data)
radar.data <- cbind(cluster = rownames(data), data)

ggradar2::ggradar2(t(data), grid.max = 0.75, grid.min = 0, centre.y = -0.05, 
                   base.size =0, axis.label.size = 5, grid.label.size = 5,
                   group.line.width = 1, group.point.size = 2, gridline.max.colour = "gray50",
                   gridline.max.linetype = 1, fullscore = rep(0.75, 7),
                   gridline.label = c(0,0.15,"", 0.45, "",0.75), legend.text.size = 14) +
  theme(legend.position = "right")




#### Cell cycle (cluster cell by stage for further analysis) ####
exp.mat <- read.table(file = "./cell_cycle/nestorawa_forcellcycle_expressionMatrix.txt", 
                      header = TRUE, as.is = TRUE, row.names = 1)
cc.genes <- readLines(con = "./cell_cycle/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

OA <- CellCycleScoring(object = OA, s.features = s.genes, g2m.features = g2m.genes, 
                       set.ident = FALSE)

DimPlot(OA, reduction = "umap", group.by = "Phase")

# Plot cell cycle and clusters
ggplot(OA@meta.data, aes(x = seurat_clusters, fill = Phase)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  xlab("cluster") + ylab("percentage") +
  theme_classic()




#####################
##### cluster 1 #####
#####################
#### Subclustering ####
mnnCT7.1 <- subset(mnnCT7, idents = "1")
DimPlot(mnnCT7.1)
# mnnCT7.1 <- FindVariableFeatures(mnnCT7.1)
# 
# 
# # Scaling the data
# mnnCT7.1.all.genes <- rownames(x = mnnCT7.1)
# mnnCT7.1 <- ScaleData(object = mnnCT7.1, features = mnnCT7.1.all.genes)
# 
# 
# # Perform linear dimensional reduction
# mnnCT7.1 <- RunPCA(object = mnnCT7.1, features = VariableFeatures(object = mnnCT7.1))
# 
# 
# ElbowPlot(object = mnnCT7.1, ndims = 20)


# mnnCT7.1.dims <- 17
# mnnCT7.1.res <- 0.2

mnnCT7.1.dims <- 17
mnnCT7.1.res <- 0.2
mnnCT7.1 <- RunUMAP(mnnCT7.1, reduction = "mnn", dims = 1:mnnCT7.1.dims)
mnnCT7.1 <- FindNeighbors(mnnCT7.1, reduction = "mnn", dims = 1:mnnCT7.1.dims)
mnnCT7.1 <- FindClusters(mnnCT7.1, resolution = mnnCT7.1.res)
DimPlot(mnnCT7.1, group.by = c("orig.ident", "ident"))


# Add subclusters
mnnCT7.1$subcluster <- paste0(mnnCT7.1$subcluster, ifelse(mnnCT7.1$seurat_clusters == 0, 2, 1))
mnnCT7@meta.data["subcluster"][colnames(mnnCT7.1),] <- mnnCT7.1$subcluster
DimPlot(mnnCT7, group.by = "subcluster", cells = grep("L31", colnames(mnnCT7), value = T), cols = col.sub)




#### DE analysis on subclusters ####
mnnCT7.1 <- ScaleData(mnnCT7.1, features = rownames(mnnCT7.1), assay = "RNA")
mnnCT7.1.markers <- FindAllMarkers(object = mnnCT7.1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                                   assay = "RNA")
mnnCT7.1.markers <- subset(mnnCT7.1.markers, p_val_adj < 0.05 & (pct.1 - pct.2) > 0)
table(mnnCT7.1.markers$cluster)


# Remove shared markers
mnnCT7.1.unique.id <- which(table(mnnCT7.1.markers$gene) == 1)
mnnCT7.1.specific.markers <- subset(mnnCT7.1.markers, gene %in% names(mnnCT7.1.unique.id))
table(mnnCT7.1.specific.markers$cluster)

# Add cluster
mnnCT7.1.specific.markers$cluster.name <- mnnCT7.1.specific.markers$cluster
mnnCT7.1.cluster <- sort(unique(mnnCT7.1$subcluster))
mnnCT7.1.specific.markers$cluster.name <- ifelse(mnnCT7.1.specific.markers$cluster == 0, 
                                                 mnnCT7.1.cluster[2], mnnCT7.1.cluster[1])


# Plot markers
mnnCT7.1.top10 <- mnnCT7.1.specific.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
my_heatmap(mnnCT7.1, genes = mnnCT7.1.top10$gene)




#### Enrichment assay ####
mnnCT7.1.go.list <- split.data.frame(mnnCT7.1.specific.markers, mnnCT7.1.specific.markers$cluster.name) %>% 
  lapply("[[", "gene") %>%
  lapply(gost)

# Extract GO results
mnnCT7.1.result.list <- 
  lapply(mnnCT7.1.go.list, function(go) {go$result %>% filter(term_size <1000) %>% dplyr::select(1:13)}) %>%
  lapply(function(x) {x$cluster.name <- (mnnCT7.1.cluster[parent.frame()$i[]]); return(x)})

mnnCT7.1.go <- do.call(rbind, mnnCT7.1.result.list)
table(mnnCT7.1.go$source, mnnCT7.1.go$cluster)




#####################
##### cluster 3 #####
#####################
#### L43.inte ####
VlnPlot(mnnCT7[, grep("L43", colnames(mnnCT7))], "cnv", group.by = "subcluster") +
  geom_hline(yintercept = 30, col = "red")

L43.inte <- mnnCT7[, grep("L43", colnames(mnnCT7))]

cnv_cluster <- as.character(L43.inte$cluster)
L43.inte$subcluster <- ifelse(cnv_cluster=="Stro" & L43.inte$cnv>30, "nStro", "Strol")
L43.inte$subcluster <- factor(L43.inte$subcluster, levels = mnnCT7.sucluster)

Idents(L43.inte) <- as.factor(L43.inte@meta.data["subcluster"])
VlnPlot(L43.inte, "cnv", group.by = "subcluster") + geom_hline(yintercept = 30, col = "red")
DimPlot(L43.inte)


# Add subcluster
mnnCT7@meta.data["subcluster"][colnames(L43.inte),] <- L43.inte$subcluster
DimPlot(mnnCT7, group.by = "subcluster", cells = grep("L43", colnames(mnnCT7), value = T), cols = col.sub)




#####################
##### cluster 4 #####
#####################
#### Subclustering ####
mnnCT7.4 <- subset(mnnCT7, idents = "4")
DimPlot(mnnCT7.4)


mnnCT7.4.dims <- 12
mnnCT7.4.res <- 0.1
mnnCT7.4 <- RunUMAP(mnnCT7.4, reduction = "mnn", dims = 1:mnnCT7.4.dims)
mnnCT7.4 <- FindNeighbors(mnnCT7.4, reduction = "mnn", dims = 1:mnnCT7.4.dims)
mnnCT7.4 <- FindClusters(mnnCT7.4, resolution = mnnCT7.4.res)
DimPlot(mnnCT7.4, group.by = c("orig.ident", "ident"))
# DimPlot(mnnCT7.4, split.by = c("patients"), group.by = "ident", label = T)

# Add subclusters
mnnCT7.4$subcluster <- paste0(mnnCT7.4$subcluster, ifelse(mnnCT7.4$seurat_clusters == 0, 1, 2))
mnnCT7@meta.data["subcluster"][colnames(mnnCT7.4),] <- mnnCT7.4$subcluster
DimPlot(mnnCT7, group.by = "subcluster", cells = grep("L44", colnames(mnnCT7), value = T), cols = col.sub)




#### DE analysis on subclusters ####
mnnCT7.4 <- ScaleData(mnnCT7.4, features = rownames(mnnCT7.4), assay = "RNA")
mnnCT7.4.markers <- FindAllMarkers(object = mnnCT7.4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                                   assay = "RNA")
mnnCT7.4.markers <- subset(mnnCT7.4.markers, p_val_adj < 0.05 & (pct.1 - pct.2) > 0)
table(mnnCT7.4.markers$cluster)


# Remove shared markers
mnnCT7.4.unique.id <- which(table(mnnCT7.4.markers$gene) == 1)
mnnCT7.4.specific.markers <- subset(mnnCT7.4.markers, gene %in% names(mnnCT7.4.unique.id))
table(mnnCT7.4.specific.markers$cluster)


# Add cluster
mnnCT7.4.specific.markers$cluster.name <- mnnCT7.4.specific.markers$cluster
mnnCT7.4.cluster <- sort(unique(mnnCT7.4$subcluster))
mnnCT7.4.specific.markers$cluster.name <- ifelse(mnnCT7.4.specific.markers$cluster == 0, 
                                                 mnnCT7.4.cluster[1], mnnCT7.4.cluster[2])



# Plot markers
mnnCT7.4.top10 <- mnnCT7.4.specific.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
my_heatmap(mnnCT7.4, genes = mnnCT7.4.top10$gene)




#### Enrichment assay ####
mnnCT7.4.go.list <- split.data.frame(mnnCT7.4.specific.markers, mnnCT7.4.specific.markers$cluster.name) %>% 
  lapply("[[", "gene") %>%
  lapply(gost)

# Extract GO results
mnnCT7.4.result.list <- 
  lapply(mnnCT7.4.go.list, function(go) {go$result %>% filter(term_size <1000) %>% dplyr::select(1:13)}) %>%
  lapply(function(x) {x$cluster.name <- (mnnCT7.4.cluster[parent.frame()$i[]]); return(x)})

mnnCT7.4.go <- do.call(rbind, mnnCT7.4.result.list)
table(mnnCT7.4.go$source, mnnCT7.4.go$cluster)




####################
##### mnnCT7.2 #####
####################
#### Subclustering ####
# Subset data
mnnCT7.2.raw <- subset(mnnCT7, idents = c(2))
mnnCT7.2.raw <- FindVariableFeatures(mnnCT7.2.raw)


# mnnCT7.2.raw.pc <- 50, mnnCT7.2.raw.dims <- 40, mnnCT7.2.raw.res <- 0.2/0.3
mnnCT7.2.raw.pc <- 50
mnnCT7.2.raw <- RunFastMNN(object.list = SplitObject(mnnCT7.2.raw, split.by = "orig.ident"), d = mnnCT7.2.raw.pc)


mnnCT7.2.raw.dims <- 40
mnnCT7.2.raw.res <- 0.3
mnnCT7.2.raw <- RunUMAP(mnnCT7.2.raw, reduction = "mnn", dims = 1:mnnCT7.2.raw.dims)
mnnCT7.2.raw <- FindNeighbors(mnnCT7.2.raw, reduction = "mnn", dims = 1:mnnCT7.2.raw.dims)
mnnCT7.2.raw <- FindClusters(mnnCT7.2.raw, resolution = mnnCT7.2.raw.res)
DimPlot(mnnCT7.2.raw, group.by = c("orig.ident", "ident"))
# DimPlot(mnnCT7.2.raw, split.by = c("orig.ident"), group.by = "ident", label = T)

# Add subclusters
names(mnnCT7.2.cluster) <- levels(mnnCT7.2.raw)

mnnCT7.2.raw$subcluster <- recode(mnnCT7.2.raw$seurat_clusters, !!!mnnCT7.2.cluster)
levels(mnnCT7.2.raw$subcluster) <- mnnCT7.2.cluster[c(1,2,4,3,5,6,7)]
DimPlot(mnnCT7.2.raw, group.by = "subcluster", cols = col.CT7.2)




#### DE analysis on subclusters ####
mnnCT7.2.raw <- ScaleData(mnnCT7.2.raw, features = rownames(mnnCT7.2.raw), assay = "RNA")
mnnCT7.2.raw.markers <- FindAllMarkers(object = mnnCT7.2.raw, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                                       assay = "RNA")
mnnCT7.2.raw.markers <- subset(mnnCT7.2.raw.markers, p_val_adj < 0.05 & (pct.1 - pct.2) > 0)
table(mnnCT7.2.raw.markers$cluster)
# write.csv(mnnCT7.2.raw.markers, paste("./integration/mnnCT7.2.raw.markers.pc",
#                                  mnnCT7.2.raw.pc, ".d", mnnCT7.2.raw.dims, ".r", mnnCT7.2.raw.res, ".csv", sep = ""))

# Remove shared markers
mnnCT7.2.raw.unique.id <- which(table(mnnCT7.2.raw.markers$gene) == 1)
mnnCT7.2.raw.specific.markers <- subset(mnnCT7.2.raw.markers, gene %in% names(mnnCT7.2.raw.unique.id))
table(mnnCT7.2.raw.specific.markers$cluster)
# write.csv(mnnCT7.2.raw.specific.markers, paste("./integration/mnnCT7.2.raw.specific.markers.pc",
#                                           mnnCT7.2.raw.pc, ".d", mnnCT7.2.raw.dims, ".r", mnnCT7.2.raw.res, ".csv", sep = ""))


# Plot markers
mnnCT7.2.raw.top10 <- mnnCT7.2.raw.specific.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
my_heatmap(mnnCT7.2.raw, genes = mnnCT7.2.raw.top10$gene)




#### Enrichment assay ####
mnnCT7.2.raw.go.list <- split.data.frame(mnnCT7.2.raw.specific.markers, mnnCT7.2.raw.specific.markers$cluster) %>% 
  lapply("[[", "gene") %>%
  lapply(gost)

# Extract GO results
mnnCT7.2.raw.result.list <- 
  lapply(mnnCT7.2.raw.go.list, function(go) {go$result %>% filter(term_size <1000) %>% dplyr::select(1:13)}) %>%
  lapply(function(x) {x$cluster <- (parent.frame()$i[]-1); return(x)})

mnnCT7.2.raw.go <- do.call(rbind, mnnCT7.2.raw.result.list)
table(mnnCT7.2.raw.go$source, mnnCT7.2.raw.go$cluster)

# write.csv(mnnCT7.2.raw.go, paste("./integration/mnnCT7.2.raw.go.pc",mnnCT7.2.raw.pc, ".d", mnnCT7.2.raw.dims, ".r", mnnCT7.2.raw.res, ".csv",
#                            sep = ""))






#### AddModuleScore ####
# Marker gene score
gene.list <- split.data.frame(mnnCT7.2.raw.specific.markers, mnnCT7.2.raw.specific.markers$cluster)  %>% 
  lapply("[[", "gene")

mnnCT7.2.raw <- AddModuleScore(mnnCT7.2.raw, gene.list, assay = "RNA")

cluster.score <- paste0("subcluster3.", names(gene.list), ".score")
colnames(mnnCT7.2.raw@meta.data)[grep("Cluster", colnames(mnnCT7.2.raw@meta.data))] <- cluster.score


cluster.score.plot.list <- mapply(function(feature, ident, col){
  (VlnPlot(mnnCT7.2.raw, features = feature, idents = ident, group.by = "orig.ident", pt.size = 0.1, cols = col) + 
     NoLegend()) %>% list()
}, cluster.score, levels(mnnCT7.2.raw), list(col.patient))

CombinePlots(cluster.score.plot.list, ncol = 4)





# Mean score matrix
score.cut.off <- 0.25
mean.score.matrix <- data.frame(patient = object.name)

for (i in 0:(length(levels(mnnCT7.2.raw))-1)) {
  scroe <- paste0("subcluster3.", i, ".score")
  mean.score <- mnnCT7.2.raw@meta.data %>%
    # tibble::rownames_to_column("cell.id") %>%
    group_by(patient) %>%
    filter(seurat_clusters == i) %>%
    summarise(mean(!!ensym(scroe)))
  colnames(mean.score)[2] <- i
  mean.score.matrix <- right_join(mean.score, mean.score.matrix)
}
mean.score.matrix


# Average cluster score 
mean.score.melt <- reshape::melt(mean.score.matrix %>% as.data.frame())
colnames(mean.score.melt)[3] <- "mean"
mean.score.melt$mean.cut.off <- mean.score.melt$mean > score.cut.off


# Distance to highest scroe
mean.score.dist.matrix <- t(rowMaxs(t(mean.score.matrix[-1]), na.rm = T) - t(mean.score.matrix[-1])) %>%
  cbind(mean.score.matrix[1])


mean.score.dist.melt <- reshape::melt(mean.score.dist.matrix %>% as.data.frame())
colnames(mean.score.dist.melt)[3] <- "dist"
mean.score.dist.melt$dist.cut.off <- mean.score.dist.melt$dist < score.cut.off

mean.score.melt <- left_join(mean.score.melt, mean.score.dist.melt)
mean.score.melt$denosed <- mean.score.melt$mean.cut.off & mean.score.melt$dist.cut.off
mean.score.melt$mark <- ifelse(mean.score.melt$denosed == T, "*", "")
mean.score.melt$patients <- factor(mean.score.melt$patient, levels = object.name)
mean.score.melt$variable <- factor(mean.score.melt$variable, levels = levels(mnnCT7.2.raw))



ggplot(mean.score.melt, aes(x = variable, y = mean, fill = patient)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept  = score.cut.off, col = "red") +
  theme_classic() +
  xlab("cluster") + ylab("average cluster score") + labs(fill = "patient")+
  geom_text(aes(label = mark), position=position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_manual(values = col.patient)



ggplot(mean.score.melt, aes(x = variable, y = dist, fill = patient)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept  = score.cut.off, col = "red") +
  theme_classic() +
  xlab("cluster") + ylab("distance to highest scroe") + labs(fill = "patient") +
  geom_text(aes(label = mark), position=position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_manual(values = col.patient)




#### Denoising ####
mean.score.denosed <- subset(mean.score.melt, denosed == T) %>% mutate(test = paste(patients, variable, sep = "-"))
cell.denosed <- paste(mnnCT7.2.raw$patient, mnnCT7.2.raw$subcluster, sep = "-") %in% mean.score.denosed$test
table(cell.denosed)



# Sebset mnnCHS7.raw with clear cells
mnnCT7.2 <- mnnCT7.2.raw[, cell.denosed]


# corrected cluster score plot
cluster.score.plot.list <- mapply(function(feature, ident, col){
  (VlnPlot(mnnCT7.2, features = feature, idents = ident, group.by = "orig.ident", pt.size = 0.1, cols = col) + 
     NoLegend()) %>% list()
}, cluster.score, levels(mnnCT7.2), list(col.patient))

CombinePlots(cluster.score.plot.list, ncol = 4)


# Denoising by CNV level
VlnPlot(mnnCT7.2, group.by = "seurat_clusters", "cnv", split.by = "patient", cols = col.patient) +
  geom_hline(yintercept = 40, col = "red")
VlnPlot(mnnCT7[, grep("L43", colnames(mnnCT7))], "cnv", group.by = "subcluster")


mnnCT7.2.cnv <- mnnCT7.2[, mnnCT7.2$cnv>40]
mnnCT7.2 <- mnnCT7.2[, mnnCT7.2$cnv<40]


###########################
#### Denoised mnnCT7.2 ####
#### DE analysis on subclusters ####
mnnCT7.2 <- ScaleData(mnnCT7.2, features = rownames(mnnCT7.2), assay = "RNA")
mnnCT7.2.markers <- FindAllMarkers(object = mnnCT7.2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                                   assay = "RNA")
mnnCT7.2.markers <- subset(mnnCT7.2.markers, p_val_adj < 0.05 & (pct.1 - pct.2) > 0)
table(mnnCT7.2.markers$cluster)
# write.csv(mnnCT7.2.markers, paste("./integration/mnnCT7.2.markers.pc",
#                                  mnnCT7.2.pc, ".d", mnnCT7.2.dims, ".r", mnnCT7.2.res, ".csv", sep = ""))

# Remove shared markers
mnnCT7.2.unique.id <- which(table(mnnCT7.2.markers$gene) == 1)
mnnCT7.2.specific.markers <- subset(mnnCT7.2.markers, gene %in% names(mnnCT7.2.unique.id))
table(mnnCT7.2.specific.markers$cluster)
# write.csv(mnnCT7.2.specific.markers, paste("./integration/mnnCT7.2.specific.markers.pc",
#                                           mnnCT7.2.pc, ".d", mnnCT7.2.dims, ".r", mnnCT7.2.res, ".csv", sep = ""))


# Plot markers
mnnCT7.2.top10 <- mnnCT7.2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
my_heatmap(mnnCT7.2, genes = mnnCT7.2.top10$gene)




#### Enrichment assay ####
mnnCT7.2.go.list <- split.data.frame(mnnCT7.2.specific.markers, mnnCT7.2.specific.markers$cluster) %>% 
  lapply("[[", "gene") %>%
  lapply(gost)

# Extract GO results
mnnCT7.2.result.list <- 
  lapply(mnnCT7.2.go.list, function(go) {go$result %>% filter(term_size <1000) %>% dplyr::select(1:13)}) %>%
  lapply(function(x) {x$cluster <- (parent.frame()$i[]-1); return(x)})

mnnCT7.2.go <- do.call(rbind, mnnCT7.2.result.list)
table(mnnCT7.2.go$source, mnnCT7.2.go$cluster)

# write.csv(mnnCT7.2.go, paste("./integration/mnnCT7.2.go.pc",mnnCT7.2.pc, ".d", 
#                                       mnnCT7.2.dims, ".r", mnnCT7.2.res, ".csv",sep = ""))




####################
##### mnnCT7.6 #####
####################
#### Subclustering ####
mnnCT7.6.raw <- subset(mnnCT7, idents = "6")
DimPlot(mnnCT7.6.raw)


# mnnCT7.6.raw.pc <- 50, mnnCT7.6.raw.dims <- 40, mnnCT7.6.raw.res <- 0.2
mnnCT7.6.raw.pc <- 50
mnnCT7.6.raw <- RunFastMNN(object.list = SplitObject(mnnCT7.6.raw, split.by = "orig.ident"), d = mnnCT7.6.raw.pc)
mnnCT7.6.raw$patient <- factor(mnnCT7.6.raw$patient, levels = patient)


mnnCT7.6.raw.dims <- 40
mnnCT7.6.raw.res <- 0.4
mnnCT7.6.raw <- RunUMAP(mnnCT7.6.raw, reduction = "mnn", dims = 1:mnnCT7.6.raw.dims)
mnnCT7.6.raw <- FindNeighbors(mnnCT7.6.raw, reduction = "mnn", dims = 1:mnnCT7.6.raw.dims)
mnnCT7.6.raw <- FindClusters(mnnCT7.6.raw, resolution = mnnCT7.6.raw.res)
DimPlot(mnnCT7.6.raw, group.by = c("orig.ident", "ident"))
DimPlot(mnnCT7.6.raw, split.by = c("patient"), group.by = "ident", label = T)


# Set cluster name
names(mnnCT7.6.cluster) <- levels(mnnCT7.6.raw)

mnnCT7.6.raw$subcluster <- recode(mnnCT7.6.raw$seurat_clusters, !!!mnnCT7.6.cluster)
mnnCT7.6.raw$subcluster <- factor(mnnCT7.6.raw$subcluster, levels = mnnCT7.6.cluster[c(1,2,5,3,4)])
DimPlot(mnnCT7.6.raw, group.by = "subcluster", cols = col.CT7.6)



#### DE analysis on subclusters ####
mnnCT7.6.raw <- ScaleData(mnnCT7.6.raw, features = rownames(mnnCT7.6.raw), assay = "RNA")
mnnCT7.6.raw.markers <- FindAllMarkers(object = mnnCT7.6.raw, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                                       assay = "RNA")
mnnCT7.6.raw.markers <- subset(mnnCT7.6.raw.markers, p_val_adj < 0.05 & (pct.1 - pct.2) > 0)
table(mnnCT7.6.raw.markers$cluster)
# write.csv(mnnCT7.6.raw.markers, paste("./integration/mnnCT7.6.raw.markers.pc",
#                                  mnnCT7.6.raw.pc, ".d", mnnCT7.6.raw.dims, ".r", mnnCT7.6.raw.res, ".csv", sep = ""))

# Remove shared markers
mnnCT7.6.raw.unique.id <- which(table(mnnCT7.6.raw.markers$gene) == 1)
mnnCT7.6.raw.specific.markers <- subset(mnnCT7.6.raw.markers, gene %in% names(mnnCT7.6.raw.unique.id))
table(mnnCT7.6.raw.specific.markers$cluster)
# write.csv(mnnCT7.6.raw.specific.markers, paste("./integration/mnnCT7.6.raw.specific.markers.pc",
#                                           mnnCT7.6.raw.pc, ".d", mnnCT7.6.raw.dims, ".r", mnnCT7.6.raw.res, ".csv", sep = ""))


# Plot markers
mnnCT7.6.raw.top10 <- mnnCT7.6.raw.specific.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
my_heatmap(mnnCT7.6.raw, genes = mnnCT7.6.raw.top10$gene)




#### Enrichment assay ####
mnnCT7.6.raw.go.list <- split.data.frame(mnnCT7.6.raw.specific.markers, mnnCT7.6.raw.specific.markers$cluster) %>% 
  lapply("[[", "gene") %>%
  lapply(gost)

# Extract GO results
mnnCT7.6.raw.result.list <- 
  lapply(mnnCT7.6.raw.go.list, function(go) {go$result %>% filter(term_size <1000) %>% dplyr::select(1:13)}) %>%
  lapply(function(x) {x$cluster <- (parent.frame()$i[]-1); return(x)})

mnnCT7.6.raw.go <- do.call(rbind, mnnCT7.6.raw.result.list)
table(mnnCT7.6.raw.go$source, mnnCT7.6.raw.go$cluster)

# write.csv(mnnCT7.6.raw.go, paste("./integration/mnnCT7.6.raw.go.pc",mnnCT7.6.raw.pc, ".d", mnnCT7.6.raw.dims, ".r", mnnCT7.6.raw.res, ".csv",
#                            sep = ""))




#### AddModuleScore ####
# Marker gene score
gene.list <- split.data.frame(mnnCT7.6.raw.specific.markers, mnnCT7.6.raw.specific.markers$cluster)  %>% 
  lapply("[[", "gene")

mnnCT7.6.raw <- AddModuleScore(mnnCT7.6.raw, gene.list, assay = "RNA")

cluster.score <- paste0("subcluster6.", names(gene.list), ".score")
colnames(mnnCT7.6.raw@meta.data)[grep("Cluster", colnames(mnnCT7.6.raw@meta.data))] <- cluster.score


cluster.score.plot.list <- mapply(function(feature, ident, col){
  (VlnPlot(mnnCT7.6.raw, features = feature, idents = ident, group.by = "orig.ident", pt.size = 0.1, cols = col) + 
     NoLegend()) %>% list()
}, cluster.score, levels(mnnCT7.6.raw), list(col.patient))

CombinePlots(cluster.score.plot.list, ncol = 4)





# Mean score matrix
score.cut.off <- 0.25
mean.score.matrix <- data.frame(patient = object.name)

for (i in 0:(length(levels(mnnCT7.6.raw))-1)) {
  scroe <- paste0("subcluster6.", i, ".score")
  mean.score <- mnnCT7.6.raw@meta.data %>%
    # tibble::rownames_to_column("cell.id") %>%
    group_by(patient) %>%
    filter(seurat_clusters == i) %>%
    summarise(mean(!!ensym(scroe)))
  colnames(mean.score)[2] <- i
  mean.score.matrix <- right_join(mean.score, mean.score.matrix)
}
mean.score.matrix


# Average cluster score 
mean.score.melt <- reshape::melt(mean.score.matrix %>% as.data.frame())
colnames(mean.score.melt)[3] <- "mean"
mean.score.melt$mean.cut.off <- mean.score.melt$mean > score.cut.off


# Distance to highest scroe
mean.score.dist.matrix <- t(rowMaxs(t(mean.score.matrix[-1]), na.rm = T) - t(mean.score.matrix[-1])) %>%
  cbind(mean.score.matrix[1])


mean.score.dist.melt <- reshape::melt(mean.score.dist.matrix %>% as.data.frame())
colnames(mean.score.dist.melt)[3] <- "dist"
mean.score.dist.melt$dist.cut.off <- mean.score.dist.melt$dist < score.cut.off

mean.score.melt <- left_join(mean.score.melt, mean.score.dist.melt)
mean.score.melt$denosed <- mean.score.melt$mean.cut.off & mean.score.melt$dist.cut.off
mean.score.melt$mark <- ifelse(mean.score.melt$denosed == T, "*", "")
mean.score.melt$patients <- factor(mean.score.melt$patient, levels = object.name)
mean.score.melt$variable <- factor(mean.score.melt$variable, levels = levels(mnnCT7.6.raw))



ggplot(mean.score.melt, aes(x = variable, y = mean, fill = patient)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept  = score.cut.off, col = "red") +
  theme_classic() +
  xlab("cluster") + ylab("average cluster score") + labs(fill = "patient")+
  geom_text(aes(label = mark), position=position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_manual(values = col.patient)



ggplot(mean.score.melt, aes(x = variable, y = dist, fill = patient)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept  = score.cut.off, col = "red") +
  theme_classic() +
  xlab("cluster") + ylab("distance to highest scroe") + labs(fill = "patient") +
  geom_text(aes(label = mark), position=position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_manual(values = col.patient)




#### Denoising ####
mean.score.denosed <- subset(mean.score.melt, denosed == T) %>% mutate(test = paste(patients, variable, sep = "-"))
cell.denosed <- paste(mnnCT7.6.raw$patient, mnnCT7.6.raw$seurat_clusters, sep = "-") %in% mean.score.denosed$test
table(cell.denosed)



# Sebset mnnCHS7.raw with clear cells
mnnCT7.6 <- mnnCT7.6.raw[, cell.denosed]


# corrected cluster score plot
cluster.score.plot.list <- mapply(function(feature, ident, col){
  (VlnPlot(mnnCT7.6, features = feature, idents = ident, group.by = "orig.ident", pt.size = 0.1, cols = col) + 
     NoLegend()) %>% list()
}, cluster.score, levels(mnnCT7.6), list(col.patient))

CombinePlots(cluster.score.plot.list, ncol = 4)


# Denoising by CNV level
VlnPlot(mnnCT7.6, group.by = "seurat_clusters", "cnv", split.by = "patient", cols = col.patient) +
  geom_hline(yintercept = 40, col = "red")
VlnPlot(mnnCT7[, grep("L43", colnames(mnnCT7))], "cnv", group.by = "subcluster")


mnnCT7.6.cnv <- mnnCT7.6[, mnnCT7.6$cnv>40]
mnnCT7.6 <- mnnCT7.6[, mnnCT7.6$cnv<40]


#### ImmCluster ####
# Save txt file
ImmCluster <- as.matrix(mnnCT7.6@assays[["RNA"]]@counts)
colnames(ImmCluster) <- str_replace(colnames(ImmCluster), "\\.", "_")
write.table(ImmCluster, "ImmCluster/mnnCT7.6.txt", sep = "\t", quote = F)


# Read ImmCluster annotation
Imm.cluster <- read.table("ImmCluster/sctype_cellPredict.txt", header = T)
rownames(Imm.cluster) <- str_split(Imm.cluster$cellname, "_", n = 2) %>% sapply(paste, collapse = ".")


mnnCT7.6$ImmCluster <- Imm.cluster[colnames(mnnCT7.6), "celltype"]
DimPlot(mnnCT7.6, group.by = "ImmCluster")





###########################
#### Denoised mnnCT7.6 ####
#### DE analysis on subclusters ####
mnnCT7.6 <- ScaleData(mnnCT7.6, features = rownames(mnnCT7.6), assay = "RNA")
mnnCT7.6.markers <- FindAllMarkers(object = mnnCT7.6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                                   assay = "RNA")
mnnCT7.6.markers <- subset(mnnCT7.6.markers, p_val_adj < 0.05 & (pct.1 - pct.2) > 0)
table(mnnCT7.6.markers$cluster)
# write.csv(mnnCT7.6.markers, paste("./integration/mnnCT7.6.markers.pc",
#                                  mnnCT7.6.pc, ".d", mnnCT7.6.dims, ".r", mnnCT7.6.res, ".csv", sep = ""))

# Remove shared markers
mnnCT7.6.unique.id <- which(table(mnnCT7.6.markers$gene) == 1)
mnnCT7.6.specific.markers <- subset(mnnCT7.6.markers, gene %in% names(mnnCT7.6.unique.id))
table(mnnCT7.6.specific.markers$cluster)
# write.csv(mnnCT7.6.specific.markers, paste("./integration/mnnCT7.6.specific.markers.pc",
#                                           mnnCT7.6.pc, ".d", mnnCT7.6.dims, ".r", mnnCT7.6.res, ".csv", sep = ""))


# Plot markers
mnnCT7.6.top10 <- mnnCT7.6.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
my_heatmap(mnnCT7.6, genes = mnnCT7.6.top10$gene)




#### Enrichment assay ####
mnnCT7.6.go.list <- split.data.frame(mnnCT7.6.specific.markers, mnnCT7.6.specific.markers$cluster) %>% 
  lapply("[[", "gene") %>%
  lapply(gost)

# Extract GO results
mnnCT7.6.result.list <- 
  lapply(mnnCT7.6.go.list, function(go) {go$result %>% filter(term_size <1000) %>% dplyr::select(1:13)}) %>%
  lapply(function(x) {x$cluster <- (parent.frame()$i[]-1); return(x)})

mnnCT7.6.go <- do.call(rbind, mnnCT7.6.result.list)
table(mnnCT7.6.go$source, mnnCT7.6.go$cluster)

# write.csv(mnnCT7.6.go, paste("./integration/mnnCT7.6.go.pc",mnnCT7.6.pc, ".d", 
#                                       mnnCT7.6.dims, ".r", mnnCT7.6.res, ".csv",sep = ""))




#####################
#### Fetal Femur ####
#####################
#### QC and clustering ####
# Load the Matrix file from Cell Ranger
FF_L75.data <- Read10X(data.dir = "./count.v3/L75")

# Create Seurat Object
FF_L75 <- CreateSeuratObject(counts = FF_L75.data, project = "FF_L75", min.cells = 3, min.features = 200)
FF_L75

# QC and selecting cells for further analysis
FF_L75[["percent.mt"]] <- PercentageFeatureSet(object = FF_L75, pattern = "^MT-")
FF_L75[["percent.rp"]] <- PercentageFeatureSet(object = FF_L75, pattern = "^RP")

# Visualize QC metrics as a violin plot
VlnPlot(object = FF_L75, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(object = FF_L75, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = FF_L75, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

FF_L75 <- subset(x = FF_L75, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt < 5)

# Normalizing the data
FF_L75 <- NormalizeData(object = FF_L75, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection)
FF_L75 <- FindVariableFeatures(object = FF_L75, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
FF_L75.top10 <- head(x = VariableFeatures(object = FF_L75), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = FF_L75)
LabelPoints(plot = plot1, points = FF_L75.top10, repel = TRUE)


# Scaling the data
FF_L75.all.genes <- rownames(x = FF_L75)
FF_L75 <- ScaleData(object = FF_L75, features = FF_L75.all.genes)

# Perform linear dimensional reduction
FF_L75 <- RunPCA(object = FF_L75, features = VariableFeatures(object = FF_L75))


ElbowPlot(object = FF_L75, ndims = 20)


# Cluster the cells
FF_L75.dims <- 20
FF_L75.res <-  0.1 #0.4 
FF_L75 <- FindNeighbors(object = FF_L75, dims = 1:FF_L75.dims)
FF_L75 <- FindClusters(object = FF_L75, resolution = FF_L75.res)

# Run non-linear dimensional reduction (UMAP/tSNE)
FF_L75 <- RunUMAP(object = FF_L75, dims = 1:FF_L75.dims)
DimPlot(object = FF_L75, reduction = "umap", label = T)
table(FF_L75@meta.data$seurat_clusters)


FF.cluster <- c("undifferentiated C", "resting C", "hypertrophic C", "proliferating C", "fibroblast", "endothelial cell")
names(FF.cluster) <- 0:5
FF_L75$cluster <- factor(FF.cluster[FF_L75$seurat_clusters], levels = FF.cluster[c(2,1,3:6)])
DimPlot(object = FF_L75, reduction = "umap", group.by = "cluster")


#### Finding differentially expressed features ####
FF_L75.markers <- FindAllMarkers(object = FF_L75, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
table(FF_L75.markers$cluster)
FF_L75@misc$markers <- FF_L75.markers


FF_L75.specific.markers <- subset(FF_L75.markers, p_val_adj < 0.05 & (pct.1 - pct.2) > 0)
FF_L75.unique.id <- which(table(FF_L75.specific.markers$gene) == 1)
FF_L75.specific.markers <- subset(FF_L75.specific.markers, gene %in% names(FF_L75.unique.id))
table(FF_L75.specific.markers$cluster)
FF_L75@misc$specific.markers <- FF_L75.specific.markers




#### Enrichment analysis ####
FF_L75.go.list <- split.data.frame(FF_L75.specific.markers, FF_L75.specific.markers$cluster) %>%
  lapply("[[", "gene") %>%
  lapply(gost)

# Extract GO results
FF_L75.result.list <- 
  lapply(FF_L75.go.list, function(go) {go$result %>% filter(term_size <1000) %>% dplyr::select(1:13)}) %>%
  lapply(function(x) {if(nrow(x) != 0) x$cluster <- (parent.frame()$i[]-1) ; return(x)})

FF_L75.go <- do.call(rbind, FF_L75.result.list)
table(FF_L75.go$source, FF_L75.go$cluster)




#### CCA (cosine similarity)####
# Calculate similarity matrix
rib.cluster <- levels(FF_L75$cluster)
chon.cluster <- intersect(levels(mnnCT7$cluster), mnnCT7$cluster)
nperm <- 10000

simatrix <- matrix(nrow = length(rib.cluster), ncol = length(chon.cluster), dimnames = list(rib.cluster, chon.cluster))
FDRmatrix <- simatrix

for (i in rib.cluster) {
  shared.gene <- intersect(rownames(mnnCT7), FF_L75@misc$specific.markers$gene[FF_L75@misc$specific.markers$cluster.name == i])
  
  rib.ctl <- FF_L75$cluster != i
  rib.dif.vec <- rowMeans(FF_L75[shared.gene,!rib.ctl]@assays$RNA@scale.data) - rowMeans(FF_L75[shared.gene,rib.ctl]@assays$RNA@scale.data)
  
  for (j in chon.cluster) {
    chon.ctl <- mnnCT7$cluster != j
    chon.dif.vec <- rowMeans(mnnCT7[shared.gene,!chon.ctl]@assays$RNA@scale.data) - rowMeans(mnnCT7[shared.gene,chon.ctl]@assays$RNA@scale.data)
    
    simatrix[i,j] <- lsa::cosine(rib.dif.vec, chon.dif.vec)
    FDRmatrix[i,j] <- FDR.test(chon.dif.vec, rib.dif.vec, nperm)
  }
}

FDRmatrix[1:4,1:5]


# Radar plot
data <- simatrix %>% as.data.frame()
data[data<0] <- 0
grid.min <- min(data)
grid.max <- max(data)

radar.data <- cbind(cluster = rownames(data), data)


ggradar2::ggradar2(t(data), grid.max = 1, grid.min = 0, centre.y = -0.05, 
                   base.size =0, axis.label.size = 5, grid.label.size = 5,
                   group.line.width = 1, group.point.size = 2, gridline.max.colour = "gray50",
                   gridline.max.linetype = 1, group.colours = col.CT7, fullscore = rep(0.75, 6),
                   webtype = "lux", gridline.label = c(0,0.15,"", 0.45, "",0.75), legend.text.size = 14) +
  theme(legend.position = "right")




################################################################################
################################################################################
####                                                                        ####
####                              Bulk RNA-seq Analysis                     ####
####                                                                        ####
################################################################################
################################################################################
set.seed(11051991)





#### Load libraries ####
#General Bioconductor packages
library(Biobase)
library(oligoClasses)

#Annotation and data import packages
library(ArrayExpress)
library(pd.hugene.2.0.st)
library(hugene20sttranscriptcluster.db)

#Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)

#Analysis and statistics packages
library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)
library(WGCNA)

#Plotting and color options packages
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(grid)
library(ggdendro)
library(scales)
library(ggdendro)

#Formatting/documentation packages
#library(rmarkdown)
#library(BiocStyle)
library(dplyr)
library(tidyr)
library(dendextend)

#Helpers:
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)
#library(devtools)

# Single cell packages
library(Seurat)
library(gprofiler2)


# Deconvolution packages
library(MuSiC)
library(xbioc)
library(bseqsc)
library(scBio)



#### Color vector ####
# Integrated object
CT7.cluster.id <- c(0:7)
CT7.cols <- brewer.pal(length(CT7.cluster.id), "Set1")
names(CT7.cols) = CT7.cluster.id
CT7.cols


# Patient
patients <- c("Ben_L49", "Low_L07", "Low_L28", "High_L63", "Med_L31", "High_L44", "Cos_L43")
patient.cols <- brewer.pal(length(patients), "Set2")
names(patient.cols) = patients
patient.cols


# Subcluster 3
CT7.3.cluster.id <- c(0:4)
CT7.3.cols <- brewer.pal(length(CT7.3.cluster.id), "Set3")
names(CT7.3.cols) = CT7.3.cluster.id
CT7.3.cols


# Histology
his <- c("benign", "G1", "G2", "G3", "dedifferentiated")
his.col <- brewer.pal(6, "Spectral")[5:1]
names(his.col) <- his
his.col


# Mutation col
mut <- c("TRUE", "FALSE", "#N/A")
mut.col <- c("red", "green", "gray")
names(mut.col) <- mut
mut.col


# Bulk group col
group <- 1:5
group.col <- brewer.pal(length(group), "Dark2")[c(1, 2, 5, 4, 3)]
names(group.col) <- group
group.col


# Heatmap col
ann_colors <- list(COL2A1.mut = mut.col, P53.mut = mut.col, IDH.mut = mut.col, 
                   histology = his.col, group = group.col)



#### Functions ####
getprobeid <- function(x) AnnotationDbi::select(hugene20sttranscriptcluster.db,
                                                keys = x,
                                                columns = "PROBEID",
                                                keytype = "SYMBOL") %>% pull(PROBEID)



getgene <- function(x) AnnotationDbi::select(hugene20sttranscriptcluster.db,
                                             keys = x,
                                             columns = "SYMBOL",
                                             keytype = "PROBEID") %>% pull(SYMBOL)




# Box plots for gene expression
bplt <- function(gene = NULL, object = NULL,  cluster = "group", ncol = NULL, nrow = NULL, col = NULL, combine = T){
  
  # Single box plot function
  bplt.single <- function(gene = NULL, object = NULL,  cluster = "group", col = NULL) {
    
    crat_expr <- Biobase::exprs(object)[gene, ]
    crat_data <- as.data.frame(crat_expr)
    colnames(crat_data)[1] <- "org_value"
    crat_data <- cbind(crat_data, Biobase::pData(object)[cluster])
    
    if(is.null(col)) col <- scales::hue_pal()(Biobase::pData(object)[, cluster] %>% unique() %>% length())
    
    ggplot(data = crat_data, aes(x = !!ensym(cluster), y = org_value, fill = !!ensym(cluster))) +
      theme_classic() +
      geom_boxplot() +
      scale_fill_manual(values = col)+
      ylab("log2(intensity)") +
      xlab(cluster) +
      ggtitle(gene) + 
      theme(plot.title = element_text(hjust=0.5)) +
      NoLegend() 
  }
  
  # Check input genes
  gene.check <- gene %in% rownames(object)
  if(0 %in% gene.check) message("Warning: ", paste0(gene[!gene.check], collapse = ","), " are not in features")
  gene <- gene[gene.check]
  if(isEmpty(gene)) break
  
  # Ploting
  gene.list <- as.list(gene)
  plt.list <- lapply(gene.list, bplt.single, object = object,  cluster = cluster, col = col)
  
  if(combine == F) {return(plt.list)} else {
    if(is.null(ncol) && is.null(nrow)) nrow <- sqrt(length(gene))%/%1
    grid.arrange(grobs = plt.list, ncol = ncol, nrow = nrow)
  }
}




#### Load data ####
bulk <- readRDS("bulk/bulk.rds")
# mnnCT7 <- readRDS("mnnCT7.rds")


# bulk_eset_norm <- readRDS("bulk_eset_norm.rds")
# colnames(bulk_eset_norm) <- str_sub(colnames(bulk_eset_norm), 5, 9)



################################################
#### Bulk mRNA microarray QC and annotation ####
################################################
#### Import of microarray expression data as "ExpressionSet" ####
# # Download data
# dir.create("E-MTAB-7264/report")
raw_data_dir <- "./bulk/E-MTAB-7264"


# anno_AE <- getAE("E-MTAB-7264", path = raw_data_dir, type = "raw")


# Create pData
bulk.meta.data <- read.csv("bulk/meta.data.csv", header = T)

sdrf_location <- file.path(raw_data_dir, "E-MTAB-7264.sdrf.txt")
SDRF <- read.delim(sdrf_location, stringsAsFactors = F)
rownames(SDRF) <- SDRF$Source.Name

SDRF <- subset(SDRF, Source.Name %in% bulk.meta.data$?..sampleID)

file.name <- SDRF[bulk.meta.data$?..sampleID, "Array.Data.File"]
rownames(bulk.meta.data) <- bulk.meta.data$Array.Data.File <-  file.name
bulk.meta.data <- subset(bulk.meta.data, !(Histology == "#N/A"))
bulk.meta.data$Histology <- factor(bulk.meta.data$Histology, levels = his)

SDRF <- AnnotatedDataFrame(bulk.meta.data)



# Read CEL files
raw_data <- oligo::read.celfiles(filenames = file.path(raw_data_dir, 
                                                       SDRF$Array.Data.File),
                                 verbose = FALSE, phenoData = SDRF)
stopifnot(validObject(raw_data))



#### Quality control of the raw data ####
# Boxplot before QC
oligo::boxplot(raw_data, target = "core", 
               main = "Boxplot of log2-intensitites for the raw data")



#### RMA calibration of the data ####
bulk_eset_norm <- oligo::rma(raw_data, target = "core")


# Boxplot on calibrated data
row_medians_assayData <- 
  Biobase::rowMedians(as.matrix(Biobase::exprs(bulk_eset_norm)))

RLE_data <- sweep(Biobase::exprs(bulk_eset_norm), 1, row_medians_assayData)

RLE_data <- as.data.frame(RLE_data)
RLE_data_gathered <- 
  tidyr::gather(RLE_data, patient_array, log2_expression_deviation)

ggplot2::ggplot(RLE_data_gathered, aes(patient_array,
                                       log2_expression_deviation)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-2, 2)) + 
  theme(axis.text.x = element_text(colour = "aquamarine4", 
                                   angle = 60, size = 6.5, hjust = 1 ,
                                   face = "bold"))


# PCA on calibrated data
exp_bulk <- Biobase::exprs(bulk_eset_norm)
PCA <- prcomp(t(exp_bulk), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Histology = 
                       Biobase::pData(bulk_eset_norm)$Histology)


ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Histology)) +
  ggtitle("PCA plot of the calibrated, summarized data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) +
  scale_color_manual(values = brewer.pal(7, "Set1"))




#### Outlier Removal ####
tree <- hclust(dist(t(exprs(bulk_eset_norm))), method="average")
plot(tree)


# Outlier plot
normadj <- (0.5 + 0.5*bicor(exprs(bulk_eset_norm)))^2
netsummary <- fundamentalNetworkConcepts(normadj)
C <- netsummary$Connectivity
Z.C <- (C-mean(C))/sqrt(var(C))

datLabel <- Biobase::pData(bulk_eset_norm)$?..sampleID
plot(1:length(Z.C),Z.C,main="Outlier Plot",xlab = "Samples",ylab="Connectivity Z Score")
text(1:length(Z.C),Z.C,label=datLabel,pos=3,cex=0.6)
abline(h= -1.5, col="red")

to_keep <- abs(Z.C) < 1.5
table(to_keep)
colnames(exprs(bulk_eset_norm))[!to_keep]


# Hclust plot shows cell to keep
dend <- as.dendrogram(tree)
colors_to_use <- as.numeric(to_keep) +2
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend) <- colors_to_use
dend.list<-as.character(Biobase::pData(bulk_eset_norm)$?..sampleID)
labels(dend)<- dend.list[order.dendrogram(dend)]
plot(dend)


# Filter outliner
bulk_eset_norm <- bulk_eset_norm[, to_keep]

# saveRDS(bulk_eset_norm, "bulk_eset_norm.rds")



#### Filtering based on intensity ####
# Histograph
bulk_medians <- rowMedians(Biobase::exprs(bulk_eset_norm))

hist_res <- hist(bulk_medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")


man_threshold <- 5

hist_res <- hist(bulk_medians, 100, col = "cornsilk", freq = FALSE, 
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

abline(v = man_threshold, col = "coral4", lwd = 2)


# Individual histogram
i=1 
plot(density((exprs(bulk_eset_norm)[,i]),na.rm=T),col = as.numeric(Biobase::pData(bulk_eset_norm)$Histology)[i],
     main = "Histogram", xlab="log2 exp",xlim=c(2,15),ylim=c(0,0.5))
for(i in 2:dim(exprs(bulk_eset_norm))[2]){
  lines(density((exprs(bulk_eset_norm)[,i]),na.rm=T), col = as.numeric(Biobase::pData(bulk_eset_norm)$Histology)[i],)
}
legend("topright",legend = levels(Biobase::pData(bulk_eset_norm)$Histology),
       fill = as.numeric(as.factor(levels(Biobase::pData(bulk_eset_norm)$Histology))))
abline(v = man_threshold, col = "coral4", lwd = 2)


# Filtering
samples_cutoff <- ncol(bulk_eset_norm)/4

idx_man_threshold <- apply(Biobase::exprs(bulk_eset_norm), 1,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)

bulk_manfiltered <- subset(bulk_eset_norm, idx_man_threshold)


# Individual histogram after filtering based on intensity
i=1 
plot(density((exprs(bulk_manfiltered)[,i]),na.rm=T),col = as.numeric(Biobase::pData(bulk_manfiltered)$Histology)[i],
     main = "Histogram", xlab="log2 exp",xlim=c(2,15),ylim=c(0,0.5))
for(i in 2:dim(exprs(bulk_manfiltered))[2]){
  lines(density((exprs(bulk_manfiltered)[,i]),na.rm=T), col = as.numeric(Biobase::pData(bulk_manfiltered)$Histology)[i],)
}
legend("topright",legend = levels(Biobase::pData(bulk_manfiltered)$Histology),
       fill = as.numeric(as.factor(levels(Biobase::pData(bulk_manfiltered)$Histology))))
abline(v = man_threshold, col = "coral4", lwd = 2)



#### Heatmap clustering analysis after QC ####
Biobase::pData(bulk_manfiltered)[Biobase::pData(bulk_manfiltered) == "#N/A"] <- NA

bulk_manfiltered$idh.mut <- bulk_manfiltered$IDH1.AAmut != "wt" | bulk_manfiltered$IDH2.AAmut != "wt"
bulk_manfiltered$col2a1.mut <- bulk_manfiltered$COL2A1  != "wt"
bulk_manfiltered$p53.mut <- (bulk_manfiltered$TP53  != "wt")

Biobase::pData(bulk_manfiltered)[is.na(Biobase::pData(bulk_manfiltered))] <- "#N/A"

annotation_for_heatmap <- 
  data.frame(COL2A1.mut=as.character(bulk_manfiltered$col2a1.mut),
             P53.mut=as.character(bulk_manfiltered$p53.mut),
             IDH.mut = as.character(bulk_manfiltered$idh.mut),
             histology = bulk_manfiltered$Histology)

row.names(annotation_for_heatmap) <- bulk_manfiltered$?..sampleID

dists <- as.matrix(dist(t(exprs(bulk_manfiltered)), method = "manhattan"))
rownames(dists) <- bulk_manfiltered$?..sampleID 

ann_colors <- list(COL2A1.mut = mut.col, P53.mut = mut.col, IDH.mut = mut.col, histology = his.col)
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))(255))
colnames(dists) <- NULL
diag(dists) <- NA

pheatmap(dists, 
         col = hmcol, 
         cutree_cols = 4,
         cutree_rows = 4,
         annotation_row = annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE, 
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = TRUE), 
                           max(dists, na.rm = TRUE)), 
         legend_labels = (c("small", "large")),
         main = "Clustering heatmap for the calibrated samples")



#### Annotation of the transcript clusters ####
anno_bulk <- AnnotationDbi::select(hugene20sttranscriptcluster.db,
                                   keys = (featureNames(bulk_manfiltered)),
                                   columns = c("SYMBOL", "GENENAME"),
                                   keytype = "PROBEID")


# Removing multiple mappings and probes without SYMBOL mapping
anno_bulk <- subset(anno_bulk, !is.na(SYMBOL))

anno_grouped <- group_by(anno_bulk, PROBEID)
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))

head(anno_summarized)


anno_filtered <- filter(anno_summarized, no_of_matches == 1)
head(anno_filtered)
nrow(anno_filtered)


# Subset bulk with annotated probes
bulk <- bulk_manfiltered[anno_filtered$PROBEID,]
bulk


# Annotate bulk_final
head(anno_bulk)

fData(bulk)$PROBEID <- rownames(fData(bulk))

fData(bulk) <- left_join(fData(bulk), anno_bulk)


rownames(fData(bulk)) <- fData(bulk)$PROBEID 
validObject(bulk)
head(fData(bulk))


# Collapse Rows: probe maps multiple gene SYMBOL
CR <- collapseRows(exprs(bulk), rowGroup = fData(bulk)$SYMBOL, rowID = fData(bulk)$PROBEID)
idx <- match(CR$group2row[,"selectedRowID"], fData(bulk)$PROBEID)
bulk <- bulk[idx,]

rownames(bulk) <- fData(bulk)$SYMBOL
colnames(bulk) <- Biobase::pData(bulk)$?..sampleID


# Save RDS of Experimentset after QC and annotation
saveRDS(bulk, "bulk.rds")



#######################
#### Deconvolution ####
#######################
#### Cibersortx ####
gene.common <- intersect(Biobase::fData(bulk)$SYMBOL, rownames(mnnCT7))

# Take signature matrix
mtx.refe <- GetAssayData(mnnCT7[gene.common, ], slot = "counts")
mtx.refe[1:5, 1:5]


subcluster <- mnnCT7@meta.data["subcluster"]
subcluster$subcluster <- subcluster$subcluster %>% as.character()
subcluster[colnames(mnnCT7.2),] <- mnnCT7.2$subcluster %>% as.character()
subcluster[colnames(mnnCT7.6),] <- mnnCT7.6$subcluster %>% as.character()
table(subcluster$subcluster)


cell_to_remove <- (subcluster$subcluster %in% c("Stro1", "Leuk")) | mnnCT7$patient == "Cos_L43"
table(cell_to_remove)
table(subcluster$subcluster[!cell_to_remove])


mtx.refe <- rbind(t(as.matrix(subcluster$subcluster[!cell_to_remove])), as.matrix(mtx.refe)[, !cell_to_remove])
mtx.refe[1:5, 1:5]

# mtx.refe[1,] <- ifelse(mtx.refe[1,] == 5, 0, mtx.refe[1,])
# table(mtx.refe[1,])

rownames(mtx.refe)[1] <- "GeneSymbol"
mtx.refe[1:5, 1:5]
table(mtx.refe[1,])

write.table(mtx.refe, file = "./cibersortx/mtx.refe.txt", quote = F, row.names = T, col.names = F, sep="\t")


# Prepare mixtur file
mtx.mix <- exprs(bulk)

mtx.mix <- rbind(colnames(mtx.mix), mtx.mix)
rownames(mtx.mix)[1] <- "Gene"

mtx.mix[1:5, 1:5]

write.table(mtx.mix, file = "./cibersortx/mtx.mix.txt", quote = F, row.names = T, col.names = F, sep="\t")


# Plot result
cibersortx <- read.csv("./cibersortx/CIBERSORTx_Job13_Results.csv", row.names = 1, header = T)[, 1:7]
colnames(cibersortx) <- str_split(colnames(cibersortx), "X") %>% sapply("[[", 2)
colnames(cibersortx) <- CT7.cluster[colnames(cibersortx)]
head(cibersortx)


data <- cbind(cibersortx[heatmap.data$tree_row$labels, ], bulk$group) %>% 
  tibble::rownames_to_column("sample.id") %>% 
  reshape::melt()


data$variable <- factor(data$variable, levels = rev(CT7.cluster))

level <- heatmap.data$tree_row$labels[heatmap.data$tree_row$order]
data$sample.id <- factor(data$sample.id, levels = level)

level <- bulk$group[heatmap.data$tree_row$order] %>% unique()
data$`bulk$group` <- factor(data$`bulk$group`, levels = level)

ggplot(data = data) +
  geom_bar(aes(x = sample.id, y = value, fill = variable), stat = "identity") +
  theme_classic() +
  facet_grid(~`bulk$group`, scales="free", space="free") +
  scale_fill_manual(values = col.CT7) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(labels = percent_format())



# Plot dist heatmap and ciborsortx
cibersortx.data <- ggplot(data = data) +
  geom_bar(aes(x = sample.id, y = value, fill = variable), stat = "identity", width = 0.85) +
  theme_void() +
  scale_fill_manual(values = col.CT7) +
  labs(fill = "cluster") +
  theme(axis.text.y = element_text(size = 8), axis.ticks.y = element_line()) +
  scale_y_continuous(labels = percent_format(), position = "right") +
  theme(legend.key.size = unit(0.3, "cm"))


heatmap.data <- pheatmap(dists, 
                         clustering_method = "ward.D2",
                         col = hmcol,
                         show_rownames = F,
                         show_colnames = F,
                         border_color = NA,
                         annotation_col = cbind(annotation_for_heatmap, group = bulk$group),
                         annotation_colors = ann_colors,
                         legend = TRUE, 
                         fontsize = 8,
                         treeheight_row = 0,
                         treeheight_col = 0,
                         legend_breaks = c(min(dists, na.rm = TRUE), 
                                           max(dists, na.rm = TRUE)), 
                         legend_labels = (c("small", "large")))


grob.title <- textGrob("Deconvolution of bulk mRNA profiles", hjust = 0.5, vjust = 0.5, gp = gpar(fontsize = 15))

cibersortx.data <- cibersortx.data + theme(plot.margin=unit(c(0,2.1,0,0.1),"cm"))

heatmap.hclust <- ggdendrogram(heatmap.data$tree_row) + theme_void() + 
  theme(plot.margin=unit(c(0,3.5,0,-0.65),"cm"))

grid.arrange(grobs = list(heatmap.hclust, cibersortx.data, heatmap.data$gtable), 
             layout_matrix = matrix(c(1,2,2,2,3,3,3,3,3,3,3,3,3,3)),
             top = grob.title)




# Average by group
data <- cbind(cibersortx[heatmap.data$tree_row$labels, ], bulk$group) %>% 
  group_by(`bulk$group`) %>%
  summarise_all(mean) %>% as.data.frame() %>%
  reshape::melt()

data$variable <- factor(data$variable, levels = 7:0)

level <- bulk$group[heatmap.data$tree_row$order] %>% unique()
data$`bulk$group` <- factor(data$`bulk$group`, levels = level)

ggplot(data = data) +
  geom_bar(aes(x = `bulk$group`, y = value, fill = variable), stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = CT7.cols) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(labels = percent_format()) +
  xlab("patient group") +
  ylab("percentage") +
  labs(fill = "cluster") +
  ggtitle("Group mean of cibersortx") 




##############################
#### NMF metagene analysis####
##############################
#### NMF analysis ####
rank <- 4
cluster  <-  5

genes <- mnnCT7@misc$specific.markers %>% filter(cluster %in% c(2,6,7)) %>% pull(gene) # 2, 6, 7, cluster = 3

res.nmf <- nmf(bulk[rownames(bulk) %in% genes,], rank = rank, nrun = 50, seed = 19911105)

consensusmap(res.nmf)

condata <- consensusmap(res.nmf)


cut <- dendextend::cutree(condata$Colv, k=cluster)
table(cut)

conheatmap <- pheatmap(res.nmf@consensus[condata$colInd, ], 
                       cluster_cols = as.hclust(condata$Colv), cluster_rows = F,
                       annotation_col = cbind(annotation_for_heatmap, rank = as.character(cut)),
                       annotation_colors = ann_colors,
                       show_rownames = F,
                       show_colnames = F)




# ciborsortx
cibersortx <- read.csv("./cibersortx/CIBERSORTx_Job13_Results.csv", row.names = 1, header = T)[, 1:7]
colnames(cibersortx) <- str_split(colnames(cibersortx), "X") %>% sapply("[[", 2)
colnames(cibersortx) <- CT7.cluster[colnames(cibersortx)]
head(cibersortx)


data <- cbind(cibersortx[colnames(bulk), ], group = factor(cut)) %>% 
  tibble::rownames_to_column("sample.id") %>% 
  reshape::melt()


data$variable <- factor(data$variable, levels = levels(mnnCT7$cluster))

level <- colnames(bulk)[condata[["rowInd"]]]
data$sample.id <- factor(data$sample.id, levels = level)

# level <- bulk$group[heatmap.data$tree_row$order] %>% unique()
# data$`bulk$group` <- factor(data$`bulk$group`, levels = level)


cibersortx.data <- ggplot(data = data) +
  geom_bar(aes(x = sample.id, y = value, fill = variable), stat = "identity", width = 0.7) +
  theme_void() +
  scale_fill_manual(values = col.all) +
  labs(fill = "cluster") +
  theme(axis.text.y = element_text(size = 8), axis.ticks.y = element_line()) +
  scale_y_continuous(labels = scales::percent_format(), position = "right") +
  theme(legend.key.size = unit(0.4, "cm"))


# Combine plots
cibersortx.data <- cibersortx.data + theme(plot.margin=unit(c(0,1.74,-0.3,0.11),"cm"))

plot_grid(cibersortx.data, conheatmap[[4]], ncol = 1, rel_heights = c(2.5,10))


# # Name group and change histology name
# group.name <- c("5" = "Low_bulk", "3" = "Med_bulk", "4" = "High_bulk", "2" = "ded_bulk" , "1" = "Stro_bulk")
# bulk$group.name <- factor(group.name[bulk$group], levels = group.name)
# 
# histology <- c("benign" = "benign", "G1" = "low", "G2" = "medium", "G3" = "high" , "dedifferentiated" = "dedifferentiated")
# bulk$histology <- factor(histology[bulk$Histology], levels = histology)




#### DE analysis ####
bulk$group <- cut %>% as.character()

pheatmap(res.nmf@consensus[condata$rowInd, condata$rowInd], 
         cluster_cols = F, cluster_rows = F,
         annotation_col = cbind(annotation_for_heatmap, rank = bulk$group),
         annotation_colors = ann_colors,
         show_rownames = F,
         show_colnames = F)



group <- as.factor(cut)
length.cut <- levels(group) %>% length()

design_bulk <- model.matrix(~ 0 + group)

head(design_bulk)


# Fitting linear model
data.lmfit <- lmFit(bulk, design_bulk)

# All pairwise comparision
contrast.matrix <- matrix(data = -1/(length.cut-1), nrow = length.cut, ncol = length.cut)
colnames(contrast.matrix) <- rownames(contrast.matrix) <- colnames(design_bulk)
diag(contrast.matrix) <- 1
contrast.matrix

data.fit.con <- contrasts.fit(data.lmfit, contrast.matrix)
data.fit.eb <- eBayes(data.fit.con)


# Extracting results 0.75 with duplicated genes
bulk.logFC <- 0.75
bulk.marker.list <- lapply(as.list(1:length.cut), function(x){topTable(data.fit.eb, number = Inf, coef = x) %>% 
    filter(logFC>bulk.logFC & adj.P.Val<0.05)})


bulk.marker.gene.list <- lapply(bulk.marker.list, "[[", 2)
lapply(bulk.marker.gene.list, length)


# Save marker gene
bulk.marker <- lapply(bulk.marker.list, function(x){
  x$group <- colnames(contrast.matrix)[parent.frame()$i[]]
  return(x)})

bulk.marker <- do.call(rbind, bulk.marker)
table(bulk.marker$group)


bulk.marker$group.name <- group.name[str_sub(bulk.marker$group, 6,6)]
table(bulk.marker$group.name)
# write.csv(bulk.marker, paste0("./bulk/bulk.specific.markers.csv"))


# Plot all markers
pheatmap(bulk[unlist(bulk.marker.gene.list[as.numeric(unique(bulk$group[condata$rowInd]))]),
              condata$rowInd],
         fontsize = 7,
         color = colorRampPalette(c("violet", "black", "yellow"))(100),
         show_rownames = F,
         annotation_colors = ann_colors,
         cluster_rows = F, cluster_cols = F, scale = "row",
         annotation_col = cbind(annotation_for_heatmap, group = bulk$group),
         breaks = seq(-2.5, 2.5, by = 5/100))


#### Enrichment assay ####
bulk.go.list <- lapply(bulk.marker.gene.list, gost)


# Extract GO results
bulk.go.result.list <- 
  lapply(bulk.go.list, function(go) {go$result %>% filter(term_size <1000) %>% dplyr::select(1:13)}) %>%
  lapply(function(x) {x$group <- (group.name[as.character(parent.frame()$i[])]); return(x)})

bulk.go <- do.call(rbind, bulk.go.result.list)
table(bulk.go$source, bulk.go$group)




#### Fisher's exact test ####
# Fisher's exact test

p.val.mtx <- matrix(nrow = length(levels(bulk$group.name)), ncol = length(levels(bulk$histology)),
                    dimnames = list(levels(bulk$group.name), levels(bulk$histology)))

for (i in rownames(p.val.mtx)) {
  group <- ifelse(bulk$group.name == i, i, "others")
  for (j in colnames(p.val.mtx)) {
    histology <- ifelse(bulk$histology == j, j, "others")
    result <- fisher.test(group, histology)
    p.val.mtx[i,j] <- result$p.value
  }
}

p.val.mtx



group <- ifelse(bulk$group == 3, "yes", "others")
histology <- ifelse(bulk$Histology %in% c("benign", "G1", "G2"), "yes", "others")
fisher.test(group, histology)



