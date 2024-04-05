# DNA methylation analysis
BiocManager::install("methyKit")
BiocManager::install("genomation")
library(methylKit)
library(genomation)
# 1. Reading the methylation call files and store them as flat file database. 
file.list = list(system.file("extdata", "test1.myCpG.txt", package = "methylKit"), 
                 system.file("extdata", "test2.myCpG.txt", package = "methylKit"), 
                 system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                 system.file("extdata", "control2.myCpG.txt", package = "methylKit")
                 )
myobj = methRead(file.list, 
                 sample.id = list("test1", "test2", "ctrl1", "ctrl2"),
                 assembly = "hg18",
                 treatment = c(1,1,0,0),
                 context = "CpG",
                 mincov = 10
                 )


myobjDB = methRead(file.list, 
                 sample.id = list("test1", "test2", "ctrl1", "ctrl2"),
                 assembly = "hg18", 
                 treatment = c(1,1,0,0), 
                 context = "CpG", 
                 dbtype = "tabix",
                 dbdir = "methylDB")

print(myobjDB[[1]]@dbpath)

# 2. Descriptive statistics on samples. 
# 2.1 Percent methylation histogram should have two peaks on both ends, any given base are either methylated or not. 
getMethylationStats(myobj[[1]], plot = T, both.strands = FALSE)
getMethylationStats(myobj[[2]], plot = T, both.strands = FALSE)
getMethylationStats(myobj[[3]], plot = T, both.strands = FALSE)
getMethylationStats(myobj[[4]], plot = T, both.strands = FALSE)
# 2.2 Histogram of coverage
getCoverageStats(myobj[[1]], plot = T, both.strands = FALSE)
# 2.3 Filtering samples based on read coverage
filtered.myobj = filterByCoverage(myobj, lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)

# 3. Comparative analysis
# 3.1 Merge samples 
meth = unite(myobj, destrand = FALSE)
head(meth)
# 3.2 Sample correlation
getCorrelation(meth, plot = T)
# 3.3 Clustering samples
clusterSamples(meth, dist = "correlation", method = "ward", plot = T)
PCASamples(meth)
# 3.4 Batch effects
sampleAnnotation = data.frame(bath_id = c("a", "a", "b", "b"), 
                              age = c(19, 34, 23, 40))
as = assocComp(mBase = meth, sampleAnnotation = sampleAnnotation)

newObj = removeComp(meth, comp = 1)
getCorrelation(newObj, plot = T)
# 3.5 Tiling windows analysis
# tiles the genome with windows of 1000bp length and 1000bp step-size and summarizes the methylation information on those tiles
myobj_lowCov = methRead(file.list, 
                        sample.id = list("test1", "test2", "ctrl1", "ctrl2"), 
                        assembly = "hg18", 
                        treatment = c(1,1,0,0),
                        context = "CpG", 
                        mincov = 3)

tiles = tileMethylCounts(myobj_lowCov, win.size = 1000, step.size = 1000, cov.bases = 10)

# 3.6 Finding differentially methylated bases or regions
myDiff = calculateDiffMeth(meth)

myDiff25p.hyper = getMethylDiff(myDiff, difference = 25, qvalue = 0.01, type = "hyper")

myDiff25p.hypo = getMethylDiff(myDiff, difference = 25, qvalue=0.01, type = "hypo")

myDiff25p = getMethylDiff(myDiff, difference = 25, qvalue = 0.01)

diffMethPerChr(myDiff, plot = T, qvalue.cutoff = 0.01, meth.cutoff = 25)

# 3.7 Correcting for overdispersion
sim.methylBase1 = dataSim(replicates = 6, sites = 1000, 
                          treatment = c(rep(1, 3), rep(0, 3)),
                          sample.ids = c(paste0("test", 1:3), paste0("ctrl", 1:3))
                          )
my.diffMeth = calculateDiffMeth(sim.methylBase1, overdispersion = "MN", test = "Chisq")
# 3.8 Accounting for covariates
covariates = data.frame(age = c(30, 80, 34, 30, 80, 40))
sim.methylBase <- dataSim(replicates = 6, sites = 1000, 
                          treatment = c(rep(1,3), rep(0,3)), 
                          covariates = covariates, 
                          sample.ids = c(paste0("test", 1:3), paste0("ctrl", 1:3)))
my.diffMeth3 <- calculateDiffMeth(sim.methylBase, covariates=covariates, 
                                  overdispersion = "MN", test = "Chisq")

# Annotation
library(genomation)
gene.obj = readTranscriptFeatures(system.file("extdata", "refseq.hg18.bed.txt", 
                                              package = "methylKit"))
annotateWithGeneParts(as(myDiff25p, "GRanges"), gene.obj)

cpg.obj = readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt", 
                                       package = "methylKit"), 
                           feature.flank.name = c("CpGi", "shores"))
diffCpGann = annotateWithFeatureFlank(as(myDiff25p, "GRanges"), 
                                      cpg.obj$CpGi, cpg.obj$shores, 
                                      feature.name = "CpGi", flank.name = "shores")
diffCpGann

