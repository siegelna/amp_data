suppressPackageStartupMessages({
    library(Seurat)
    library(celldex)
    library(ShinyCell)
    library(dplyr)
    library(SingleR)
    library(dplyr)
})

setwd("/data/user/projects/amp_data/phase2")

load(file=file.path("objects", "02.rda"))

# Columns of interest
obj@meta.data$Stim <- obj@meta.data$Type
obj@meta.data$Tissue_Type <- obj@meta.data$broad.group
obj@meta.data$Fine_Cell_Type <- obj@meta.data$fine.type
obj@meta.data$Broad_Cell_Type <- obj@meta.data$broad.type

columns_to_keep <-  c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'Stim', 'Tissue_Type', 'Fine_Cell_Type', 'Broad_Cell_Type')
obj@meta.data <- obj@meta.data[, columns_to_keep, drop = FALSE]

# Subset cell types
bcell <- subset(obj, subset = Broad_Cell_Type == "B Cell")
dn <- subset(obj, subset = Broad_Cell_Type == "DN")
glom <- subset(obj, subset = Broad_Cell_Type == "GLOM")
intl <- subset(obj, subset = Broad_Cell_Type == "INTL")
loh <- subset(obj, subset = Broad_Cell_Type == "LOH")
myeloid <- subset(obj, subset = Broad_Cell_Type == "Myeloid Cell")
nk <- subset(obj, subset = Broad_Cell_Type == "NK Cell")
plasma <- subset(obj, subset = Broad_Cell_Type == "Plasma Cell")
pt <- subset(obj, subset = Broad_Cell_Type == "PT")
tcell <- subset(obj, subset = Broad_Cell_Type == "T Cell")
und <- subset(obj, subset = Broad_Cell_Type == "UND")

rm(obj)

# Clustering

bcell <- RunPCA(bcell, npcs = 30, verbose = TRUE)
bcell <- RunUMAP(bcell, reduction = "pca", dims = 1:20)
bcell <- FindNeighbors(bcell, reduction = "pca", dims = 1:20)
bcell <- FindClusters(bcell, resolution = 0.5)
save(glom, file=file.path("objects", "bcell.rda"))

dn <- RunPCA(dn, npcs = 30, verbose = TRUE)
dn <- RunUMAP(dn, reduction = "pca", dims = 1:20)
dn <- FindNeighbors(dn, reduction = "pca", dims = 1:20)
dn <- FindClusters(dn, resolution = 0.5)
save(glom, file=file.path("objects", "dn.rda"))

glom <- RunPCA(glom, npcs = 30, verbose = TRUE)
glom <- RunUMAP(glom, reduction = "pca", dims = 1:20)
glom <- FindNeighbors(glom, reduction = "pca", dims = 1:20)
glom <- FindClusters(glom, resolution = 0.5)
save(glom, file=file.path("objects", "glom.rda"))

intl <- RunPCA(intl, npcs = 30, verbose = TRUE)
intl <- RunUMAP(intl, reduction = "pca", dims = 1:20)
intl <- FindNeighbors(intl, reduction = "pca", dims = 1:20)
intl <- FindClusters(intl, resolution = 0.5)
save(intl, file=file.path("objects", "intl.rda"))

loh <- RunPCA(loh, npcs = 30, verbose = TRUE)
loh <- RunUMAP(loh, reduction = "pca", dims = 1:20)
loh <- FindNeighbors(loh, reduction = "pca", dims = 1:20)
loh <- FindClusters(loh, resolution = 0.5)
save(loh, file=file.path("objects", "loh.rda"))

myeloid <- RunPCA(myeloid, npcs = 30, verbose = TRUE)
myeloid <- RunUMAP(myeloid, reduction = "pca", dims = 1:20)
myeloid <- FindNeighbors(myeloid, reduction = "pca", dims = 1:20)
myeloid <- FindClusters(myeloid, resolution = 0.5)
save(myeloid, file=file.path("objects", "myeloid.rda"))

nk <- RunPCA(nk, npcs = 30, verbose = TRUE)
nk <- RunUMAP(nk, reduction = "pca", dims = 1:20)
nk <- FindNeighbors(nk, reduction = "pca", dims = 1:20)
nk <- FindClusters(nk, resolution = 0.5)
save(nk, file=file.path("objects", "nk.rda"))

plasma <- RunPCA(plasma, npcs = 30, verbose = TRUE)
plasma <- RunUMAP(plasma, reduction = "pca", dims = 1:20)
plasma <- FindNeighbors(plasma, reduction = "pca", dims = 1:20)
plasma <- FindClusters(plasma, resolution = 0.5)
save(plasma, file=file.path("objects", "plasma.rda"))

pt <- RunPCA(pt, npcs = 30, verbose = TRUE)
pt <- RunUMAP(pt, reduction = "pca", dims = 1:20)
pt <- FindNeighbors(pt, reduction = "pca", dims = 1:20)
pt <- FindClusters(pt, resolution = 0.5)
save(pt, file=file.path("objects", "pt.rda"))

tcell <- RunPCA(tcell, npcs = 30, verbose = TRUE)
tcell <- RunUMAP(tcell, reduction = "pca", dims = 1:20)
tcell <- FindNeighbors(tcell, reduction = "pca", dims = 1:20)
tcell <- FindClusters(tcell, resolution = 0.5)
save(tcell, file=file.path("objects", "tcell.rda"))

und <- RunPCA(und, npcs = 30, verbose = TRUE)
und <- RunUMAP(und, reduction = "pca", dims = 1:20)
und <- FindNeighbors(und, reduction = "pca", dims = 1:20)
und <- FindClusters(und, resolution = 0.5)
save(und, file=file.path("objects", "und.rda"))

# Load annotation
monaco.ref <- celldex::MonacoImmuneData()
blueprint.ref <- celldex::BlueprintEncodeData()

# bcell
## Blueprint annnotation
load(file=file.path("objects", "bcell.rda"))
sce <- LayerData(bcell)
blueprint.main <- SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.main)
bcell@meta.data$blueprint.main <- blueprint.main$pruned.labels
save(bcell, file=file.path("objects", "bcell04.rda"))

## Monaco annotation
sce <- LayerData(bcell)
monaco.fine <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.fine)
bcell@meta.data$monaco.fine <- monaco.fine$pruned.labels
save(bcell, file=file.path('objects', 'bcell05.rda'))

# dn
## Blueprint annnotation
load(file=file.path("objects", "dn.rda"))
sce <- LayerData(dn)
blueprint.main <- SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.main)
dn@meta.data$blueprint.main <- blueprint.main$pruned.labels
save(dn, file=file.path("objects", "dn04.rda"))

# glom
## Blueprint annnotation
load(file=file.path("objects", "glom.rda"))
sce <- LayerData(glom)
blueprint.main <- SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.main)
glom@meta.data$blueprint.main <- blueprint.main$pruned.labels
save(glom, file=file.path("objects", "glom04.rda"))


# intl
## Blueprint annnotation
load(file=file.path("objects", "intl.rda"))
sce <- LayerData(intl)
blueprint.main <- SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.main)
intl@meta.data$blueprint.main <- blueprint.main$pruned.labels
save(intl, file=file.path("objects", "intl04.rda"))


# loh
## Blueprint annnotation
load(file=file.path("objects", "loh.rda"))
sce <- LayerData(loh)
blueprint.main <- SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.main)
loh@meta.data$blueprint.main <- blueprint.main$pruned.labels
save(loh, file=file.path("objects", "loh04.rda"))


# myeloid
## Blueprint annnotation
load(file=file.path("objects", "myeloid.rda"))
sce <- LayerData(myeloid)
blueprint.main <- SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.main)
myeloid@meta.data$blueprint.main <- blueprint.main$pruned.labels
save(myeloid, file=file.path("objects", "myeloid04.rda"))

## Monaco annotation
sce <- LayerData(myeloid)
monaco.fine <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.fine)
myeloid@meta.data$monaco.fine <- monaco.fine$pruned.labels
save(myeloid, file=file.path('objects', 'myeloid05.rda'))


# nk
## Blueprint annnotation
load(file=file.path("objects", "nk.rda"))
sce <- LayerData(nk)
blueprint.main <- SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.main)
nk@meta.data$blueprint.main <- blueprint.main$pruned.labels
save(nk, file=file.path("objects", "nk04.rda"))

## Monaco annotation
sce <- LayerData(nk)
monaco.fine <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.fine)
nk@meta.data$monaco.fine <- monaco.fine$pruned.labels
save(nk, file=file.path('objects', 'nk05.rda'))


# plasma
## Blueprint annnotation
load(file=file.path("objects", "plasma.rda"))
sce <- LayerData(plasma)
blueprint.main <- SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.main)
plasma@meta.data$blueprint.main <- blueprint.main$pruned.labels
save(plasma, file=file.path("objects", "plasma04.rda"))

## Monaco annotation
sce <- LayerData(plasma)
monaco.fine <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.fine)
plasma@meta.data$monaco.fine <- monaco.fine$pruned.labels
save(plasma, file=file.path('objects', 'plasma05.rda'))


# pt
## Blueprint annnotation
load(file=file.path("objects", "pt.rda"))
sce <- LayerData(pt)
blueprint.main <- SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.main)
pt@meta.data$blueprint.main <- blueprint.main$pruned.labels
save(pt, file=file.path("objects", "pt04.rda"))


# tcell
## Blueprint annnotation
load(file=file.path("objects", "tcell.rda"))
sce <- LayerData(tcell)
blueprint.main <- SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.main)
tcell@meta.data$blueprint.main <- blueprint.main$pruned.labels
save(tcell, file=file.path("objects", "tcell04.rda"))

## Monaco annotation
sce <- LayerData(tcell)
monaco.fine <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.fine)
tcell@meta.data$monaco.fine <- monaco.fine$pruned.labels
save(tcell, file=file.path('objects', 'tcell05.rda'))


# und
## Blueprint annnotation
load(file=file.path("objects", "und.rda"))
sce <- LayerData(und)
blueprint.main <- SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.main)
und@meta.data$blueprint.main <- blueprint.main$pruned.labels
save(und, file=file.path("objects", "und04.rda"))


### Write shiny files

#bcell
seu =  bcell
scConf1 = createConfig(seu)
makeShinyFiles(seu, scConf1,
             gene.mapping = TRUE,
             shiny.prefix = "sc1",
             shiny.dir = "kidney_amp2_sle_scRNAseq/",
             gex.assay = "SCT",
             default.multigene = present_genes)

#dn
seu =  dn
scConf2 = createConfig(seu)
makeShinyFiles(seu, scConf2,
             gene.mapping = TRUE,
             shiny.prefix = "sc2",
             shiny.dir = "kidney_amp2_sle_scRNAseq/",
             gex.assay = "SCT",
             default.multigene = present_genes)

#glom
seu =  glom
scConf3 = createConfig(seu)
makeShinyFiles(seu, scConf3,
             gene.mapping = TRUE,
             shiny.prefix = "sc3",
             shiny.dir = "kidney_amp2_sle_scRNAseq/",
             gex.assay = "SCT",
             default.multigene = present_genes)

#intl
seu =  intl
scConf4 = createConfig(seu)
makeShinyFiles(seu, scConf4,
             gene.mapping = TRUE,
             shiny.prefix = "sc4",
             shiny.dir = "kidney_amp2_sle_scRNAseq/",
             gex.assay = "SCT",
             default.multigene = present_genes)

#loh
seu =  loh
scConf5 = createConfig(seu)
makeShinyFiles(seu, scConf5,
             gene.mapping = TRUE,
             shiny.prefix = "sc5",
             shiny.dir = "kidney_amp2_sle_scRNAseq/",
             gex.assay = "SCT",
             default.multigene = present_genes)

#myeloid
seu =  myeloid
scConf6 = createConfig(seu)
makeShinyFiles(seu, scConf6,
             gene.mapping = TRUE,
             shiny.prefix = "sc6",
             shiny.dir = "kidney_amp2_sle_scRNAseq/",
             gex.assay = "SCT",
             default.multigene = present_genes)

#nk
seu =  nk
scConf7 = createConfig(seu)
makeShinyFiles(seu, scConf7,
             gene.mapping = TRUE,
             shiny.prefix = "sc7",
             shiny.dir = "kidney_amp2_sle_scRNAseq/",
             gex.assay = "SCT",
             default.multigene = present_genes)

#plasma
seu =  plasma
scConf8 = createConfig(seu)
makeShinyFiles(seu, scConf8,
             gene.mapping = TRUE,
             shiny.prefix = "sc8",
             shiny.dir = "kidney_amp2_sle_scRNAseq/",
             gex.assay = "SCT",
             default.multigene = present_genes)

#pt
seu =  pt
scConf9 = createConfig(seu)
makeShinyFiles(seu, scConf9,
             gene.mapping = TRUE,
             shiny.prefix = "sc9",
             shiny.dir = "kidney_amp2_sle_scRNAseq/",
             gex.assay = "SCT",
             default.multigene = present_genes)

#tcell
seu =  tcell
scConf10 = createConfig(seu)
makeShinyFiles(seu, scConf10,
             gene.mapping = TRUE,
             shiny.prefix = "sc10",
             shiny.dir = "kidney_amp2_sle_scRNAseq/",
             gex.assay = "SCT",
             default.multigene = present_genes)

#und   
seu =  und
scConf11 = createConfig(seu)
makeShinyFiles(seu, scConf11,
             gene.mapping = TRUE,
             shiny.prefix = "sc11",
             shiny.dir = "kidney_amp2_sle_scRNAseq/",
             gex.assay = "SCT",
             default.multigene = present_genes)


makeShinyCodesMulti(
  shiny.title = "AMP Phase II SLE", shiny.footnotes = NULL,
  shiny.prefix = c("sc1", "sc2", "sc3", "sc4", "sc5", "sc6", "sc7", "sc8", "sc9", "sc10", "sc11"),
  shiny.headers = c("B Cell", "DN", "GLOM", "INTL", "LOH", "Myeloid", "NK", "Plasma Cell", "PT", "T Cell", "UND"), 
  shiny.dir = "kidney_amp2_sle_scRNAseq/") 