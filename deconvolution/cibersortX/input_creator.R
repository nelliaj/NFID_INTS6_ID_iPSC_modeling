
library(Seurat)

rawCounts <- read.table("deconvolution/countmatrix_NFID.txt", header=T)

metadata_bulk <- read.table("deconvolution/metadata.txt", header=T)
metadata_bulk <- metadata_bulk[metadata_bulk$tp=="d28",]

dim(metadata_bulk)
# [1] 27 15

rownames(rawCounts) <- rawCounts$Geneid
rawCounts$Geneid <- NULL

rawCounts <- rawCounts[,colnames(rawCounts) %in% metadata_bulk$seqsample]
rawCounts <- rawCounts[rowSums(rawCounts)!=0,]

dim(rawCounts)
# [1] 70242    27

library(org.Hs.eg.db)

annot_df <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(rawCounts), column="SYMBOL", keytype="ENSEMBL")
dup_ids <- unique(annot_df$ENSEMBL[which(duplicated(annot_df$ENSEMBL))])


## Read Seurat object
obj <- readRDS("deconvolution/merged_seurat_coculture_subset_v3.rds")
meta_obj <- read.table("deconvolution/metadata_scanpy_pau_040426.txt")

obj$leiden <- meta_obj[match(rownames(obj@meta.data), rownames(meta_obj)),]$leiden

assign_df <- read.table("deconvolution/assignments_snapseed_pau_040426.txt", sep="\t", header=T)

metrics_lv1 <- read.table("deconvolution/metrics_level1_snapseed_pau_040426.txt", sep="\t", header=T)
metrics_lv2 <- read.table("deconvolution/metrics_level2_snapseed_pau_040426.txt", sep="\t", header=T)
metrics_lv3 <- read.table("deconvolution/metrics_level3_snapseed_pau_040426.txt", sep="\t", header=T)
metrics_lv4 <- read.table("deconvolution/metrics_level4_snapseed_pau_040426.txt", sep="\t", header=T)
metrics_lv5 <- read.table("deconvolution/metrics_level5_snapseed_pau_040426.txt", sep="\t", header=T)

merging <- assign_df[,colnames(assign_df)[grepl("level_", colnames(assign_df))]]
merging_mask <- t(apply(merging, 1, nchar))>0

assign_df$merged <- sapply(1:nrow(merging), function(y){
  paste(merging[y,merging_mask[y,]], collapse=",")
}, simplify=T)

obj <- obj[,obj$orig.ident=="d28_ngn2glia"]


dim(obj@assays$RNA["counts"])
# [1] 29972 28734



obj$merged_levels <- assign_df[match(obj$leiden, assign_df$leiden),]$merged


abbrev_vector <- c("neural_progenitor_cell,diencephalic_npc"="NPC_D",
                   "neural_progenitor_cell,oRG"="NPC_oRG",
                   "neural_progenitor_cell,rhombencephalic_npc,cerebellar_npc"="NPC_R_cereb",
                   "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,diencephalic_excitatory_neuron"="EN_NT_D",
                   "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,diencephalic_excitatory_neuron,hypothalamic_excitatory_neuron"="EN_NT_D_hyp",
                   "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,mesencephalic_excitatory_neuron"="EN_NT_M",
                   "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,rhombencephalic_excitatory_neuron,cerebellar_excitatory_neuron"="EN_NT_R_cereb",
                   "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,rhombencephalic_excitatory_neuron,pons_excitatory_neuron"="EN_NT_R_PONS",
                   "pericytes"="Pericytes",
                   "pns_neurons"="PNS")

obj$abbreviation <- unname(abbrev_vector[obj$merged_levels])

abbrev_vector2 <- c("neural_progenitor_cell,diencephalic_npc"="NPC",
                   "neural_progenitor_cell,oRG"="NPC",
                   "neural_progenitor_cell,rhombencephalic_npc,cerebellar_npc"="NPC",
                   "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,diencephalic_excitatory_neuron"="EN",
                   "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,diencephalic_excitatory_neuron,hypothalamic_excitatory_neuron"="EN",
                   "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,mesencephalic_excitatory_neuron"="EN",
                   "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,rhombencephalic_excitatory_neuron,cerebellar_excitatory_neuron"="EN",
                   "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,rhombencephalic_excitatory_neuron,pons_excitatory_neuron"="EN",
                   "pericytes"="Pericytes",
                   "pns_neurons"="PNS")

obj$abbreviation2 <- unname(abbrev_vector2[obj$merged_levels])


##################################
## bunch of code for collapsing ##
##################################


unique(paste0(obj$merged_levels,"---->", obj$abbreviation))

# [1] "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,rhombencephalic_excitatory_neuron,pons_excitatory_neuron---->EN_NT_R_PONS"       
# [2] "pns_neurons---->PNS"                                                                                                                          
# [3] "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,diencephalic_excitatory_neuron,hypothalamic_excitatory_neuron---->EN_NT_D_hyp"   
# [4] "neural_progenitor_cell,oRG---->NPC_oRG"                                                                                                       
# [5] "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,mesencephalic_excitatory_neuron---->EN_NT_M"                                     
# [6] "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,rhombencephalic_excitatory_neuron,cerebellar_excitatory_neuron---->EN_NT_R_cereb"
# [7] "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,diencephalic_excitatory_neuron---->EN_NT_D"                                      
# [8] "neural_progenitor_cell,rhombencephalic_npc,cerebellar_npc---->NPC_R_cereb"                                                                    
# [9] "pericytes---->Pericytes"                                                                                                                      
# [10] "neural_progenitor_cell,diencephalic_npc---->NPC_D" 

##################################
##################################
##################################


unique(paste0(obj$merged_levels,"---->", obj$abbreviation2))
# [1] "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,rhombencephalic_excitatory_neuron,pons_excitatory_neuron---->EN"      
# [2] "pns_neurons---->PNS"                                                                                                               
# [3] "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,diencephalic_excitatory_neuron,hypothalamic_excitatory_neuron---->EN" 
# [4] "neural_progenitor_cell,oRG---->NPC"                                                                                                
# [5] "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,mesencephalic_excitatory_neuron---->EN"                               
# [6] "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,rhombencephalic_excitatory_neuron,cerebellar_excitatory_neuron---->EN"
# [7] "neuron,excitatory_neuron,non_telencephalic_excitatory_neuron,diencephalic_excitatory_neuron---->EN"                                
# [8] "neural_progenitor_cell,rhombencephalic_npc,cerebellar_npc---->NPC"                                                                 
# [9] "pericytes---->Pericytes"                                                                                                           
# [10] "neural_progenitor_cell,diencephalic_npc---->NPC"                                      "neural_progenitor_cell_diencephalic_npc---->NPC_D"      

stopifnot(all(colnames(obj)==rownames(obj@meta.data)))

table(obj$abbreviation2)
# EN       NPC Pericytes       PNS 
# 19986       842       233      7673 

table(obj$abbreviation)
# EN_NT_D   EN_NT_D_hyp       EN_NT_M EN_NT_R_cereb  EN_NT_R_PONS         NPC_D       NPC_oRG   NPC_R_cereb     Pericytes           PNS 
# 3281          2535          4361          1201          8608            11           349           482           233          7673 

annot_df <- annot_df[!is.na(annot_df$SYMBOL),]

dup_ids <- unique(annot_df$SYMBOL[which(duplicated(annot_df$SYMBOL))])


for (i in dup_ids){
  # if there are duplications with symbol id, just keep the ensemble id with highest expression
  tmp <- annot_df[annot_df$SYMBOL==i,]
  id_to_keep <- names(which.max(rowSums(rawCounts)[tmp$ENSEMBL]))
  id_to_remove <- tmp[tmp$ENSEMBL!=id_to_keep,]
  annot_df <- annot_df[!annot_df$ENSEMBL %in% id_to_remove$ENSEMBL,]
}


dup_ids_ensembl <- unique(annot_df$ENSEMBL[which(duplicated(annot_df$ENSEMBL))])

## remove all duplicates from ensembles ids
annot_df <- annot_df[is.na(match(annot_df$ENSEMBL, dup_ids_ensembl)),]
stopifnot(!any(duplicated(annot_df$SYMBOL)))
stopifnot(!any(duplicated(annot_df$ENSEMBL)))

rawCounts <- rawCounts[rownames(rawCounts) %in% annot_df$ENSEMBL,]
rownames(rawCounts) <- annot_df[match(rownames(rawCounts), annot_df$ENSEMBL),]$SYMBOL

write.table(rawCounts, file="deconvolution/cibersortx/input/d28_rawCounts_bulkRNAseq.txt",
            quote=F, row.names=T, col.names=T, sep="\t")


table(rownames(obj) %in% rownames(rawCounts))
# FALSE  TRUE 
# 11544 18428 


obj <- obj[rownames(obj) %in% rownames(rawCounts),]
matrixCounts_to_export <- obj@assays$RNA["counts"]

#abbreviation

colnames(matrixCounts_to_export) <- unname(obj$abbreviation[match(colnames(matrixCounts_to_export), rownames(obj@meta.data))])

matrixCounts_to_export <- as.matrix(matrixCounts_to_export)

write.table(
  matrixCounts_to_export,
  file = "deconvolution/cibersortx/input/d28_snapseed_scRNAseq.txt",
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = TRUE
)

#abbreviation2

matrixCounts_to_export <- obj@assays$RNA["counts"]

colnames(matrixCounts_to_export) <- unname(obj$abbreviation2[match(colnames(matrixCounts_to_export), rownames(obj@meta.data))])

matrixCounts_to_export <- as.matrix(matrixCounts_to_export)

write.table(
  matrixCounts_to_export,
  file = "deconvolution/cibersortx/input/d28_snapseed_scRNAseq_collapsed.txt",
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = TRUE
)


