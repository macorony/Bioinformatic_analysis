library(ELMER)
library(ELMER.data)
library(MultiAssayExperiment)
# Data input
distal.probes <- get.feature.probe(genome = "hg19", 
                                   met.platform = "450K", 
                                   rm.chr = paste0("chr", c(2:22, "X", "Y")))
data(LUSC_RNA_refined, package = "ELMER.data")
data(LUSC_meth_refined, package = "ELMER.data")
GeneExp
Meth

# Create MAE object
mae <- createMAE(exp = GeneExp, met = Meth, save = T, linearize.exp = TRUE, save.filename = "mae.rda", 
                 met.platform = "450K", genome = "hg19", TCGA = T)
as.data.frame(colData(mae)[1:5,])
as.data.frame(sampleMap(mae)[1:5,])

# Non-TCGA data
met <- matrix(rep(0,15), ncol = 5)
colnames(met) <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")
rownames(met) <- c("cg26928153", "cg16269199", "cg13869341")

exp <- matrix(rep(0,15), ncol = 5)
colnames(exp) <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")
rownames(exp) <- c("ENSG00000073282", "ENSG00000078900", "ENSG00000141510")

assay <- c(
  rep("DNA methylation", ncol(met)), 
  rep("Gene expression", ncol(exp))
)
primary <- c(colnames(met), colnames(exp))
colname <- c(colnames(met), colnames(exp))
sampleMap <- data.frame(assay, primary, colname)

distal.probes <- get.feature.probe(
  genome = "hg19", met.platform = "EPIC"
)

colData <- DataFrame(primary = colnames(met),
                     sample = colnames(met))

mae <- createMAE(exp = exp, met = met, save = T, filter.probes = distal.probes, colData = colData,
                 sampleMap = sampleMap, linearize.exp = T, save.filename = "mae.rda", met.platform = "EPIC", 
                 genome = "hg19", TCGA = F)


# Identifying differentially methylated probes

sig.diff <- get.diff.meth(data = mae, group.col = "definition",
                          group1 = "Primary solid Tumor",
                          group2 = "Solid Tissue Normal", 
                          minSubgroupFrac = 0.2, 
                          sig.dif = 0.3, 
                          diff.dir = "hypo", 
                          cores = 1, 
                          dir.out = "result", 
                          pvalue = 0.01
                          )
head(sig.diff)

