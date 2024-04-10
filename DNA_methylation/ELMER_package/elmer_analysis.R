library(ELMER)
library(ELMER.data)
library(MultiAssayExperiment)
library(GenomicInteractions)
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

mae <- get(load("mae.rda"))

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

# Identifying putative probe-gene pairs

nearGenes <- GetNearGenes(data = mae, probes = sig.diff$probe, numFlankingGenes = 20)

hypo.pair <- get.pair(data = mae, group.col = "definition",
                      group1 = "Primary solid Tumor", group2 = "Solid Tissue Normal", 
                      nearGenes = nearGenes, mode = "unsupervised", 
                      permu.dir = "result/permu", permu.size = 100, raw.pvalue = 0.05, 
                      Pe = 0.01, filter.probes = T, filter.percentage = 0.05, 
                      filter.portion = 0.3, dir.out = "result", cores = 1, label = "hypo"
                      )
hypo.pair

# Motif enrichment analysis on the selected probes
enriched.motif <- get.enriched.motif(data = mae, probes = hypo.pair$Probe, 
                                     dir.out = "result", labe = "hypo", min.incidence = 10,
                                     lower.OR = 1.1)

# Identifying regulatory TFS
TF <- get.TFs(data = mae, group.col = "definition", 
              group1 = "Primary solid Tumor", group2 = "Solid Tissue Normal", 
              mode = "unsupervised", enriched.motif = enriched.motif, 
              dir.out = "result", cores = 1, label = "hypo", save.plots = T)

# Plot

scatter.plot(data = mae, 
             byProbe = list(probe = c("cg19403323"), numFlankingGenes = 20),
             category = "definition", 
             lm = T, save = F)

scatter.plot(data = mae, 
             byPair = list(probe = c("cg19403323"), gene = c("ENSG00000143469")), 
             category = "definition", save = T, lm_line = T)

enriched.motif[[names(enriched.motif)[1]]]

scatter.plot(data = mae, 
             byTF = list(
               TF = c("TP53", "SOX2"), 
               probe = enriched.motif[[names(enriched.motif)[1]]]
               ),
             category = "definition", 
             save = T, 
             lm_line = T)

schematic.plot(pair = hypo.pair, data = mae, group.col = "definition", 
               byProbe = hypo.pair$Probe[1], save = F)

schematic.plot(pair = hypo.pair, data = mae, group.col = "definition", 
               byGene = hypo.pair$GeneID[1], save = F)


heatmapPairs(data = mae, 
             group.col = "definition", 
             group1 = "Primary solid Tumor", 
             annotation.col = c("years_smoked", "gender"),
             group2 = "Solid Tissue Normal", 
             pairs = hypo.pair,
             filename = NULL)
