
### GBM 450k dataset

The 450k methylation array dataset is from [Strum et al., 2012](http://dx.doi.org/10.1016/j.ccr.2012.08.024). 
The dataset is available from GEO database with id [GSE36278](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36278).

Following files are needed for the analysis:

- [450K_annotation.txt](https://jokergoo.github.io/cola_examples/GBM_450K/450K_annotation.txt). Phenotype data is from [Sturm et al., 2012](http://dx.doi.org/10.1016/j.ccr.2012.08.024), supplementary table S1.

HTML reports for cola analysis:

- [GBM_450K_cgi_all_subgroup_cola_report](https://jokergoo.github.io/cola_examples/GBM_450K/GBM_450K_cgi_all_subgroup_cola_report/cola_report.html)
- [GBM_450K_cgi_island_subgroup_cola_report](https://jokergoo.github.io/cola_examples/GBM_450K/GBM_450K_cgi_island_subgroup_cola_report/cola_report.html)
- [GBM_450K_cgi_shore_subgroup_cola_report](https://jokergoo.github.io/cola_examples/GBM_450K/GBM_450K_cgi_shore_subgroup_cola_report/cola_report.html)
- [GBM_450K_cgi_sea_subgroup_cola_report](https://jokergoo.github.io/cola_examples/GBM_450K/GBM_450K_cgi_sea_subgroup_cola_report/cola_report.html)

Following code performs the analysis.

```r
library(ComplexHeatmap)
library(matrixStats)
library(circlize)
library(RColorBrewer)
library(GenomicRanges)

library(cola)
```

The methylation profiles have been measured by Illumina HumanMethylation450 BeadChip arrays.
First,  load probe data via [the **IlluminaHumanMethylation450kanno.ilmn12.hg19** package](https://bioconductor.org/packages/IlluminaHumanMethylation450kanno.ilmn12.hg19/).

```r
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data("IlluminaHumanMethylation450kanno.ilmn12.hg19",
    package = "IlluminaHumanMethylation450kanno.ilmn12.hg19")
probe = IlluminaHumanMethylation450kanno.ilmn12.hg19 # change to a short name
```

Methylation profiles can be download from [GEO database](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36278).
The [**GEOquery** package](http://bioconductor.org/packages/release/bioc/html/GEOquery.html) is used to retrieve expression data from GEO.


```r
library(GEOquery)
if(file.exists("GSE36278_450K.RData")) {
    load("GSE36278_450K.RData")
} else {
    gset = getGEO("GSE36278")
    save(gset, file = "GSE36278_450K.RData")
}
```

Adjust row names in the matrix to be the same as the probes. 

```r
mat = exprs(gset[[1]])
colnames(mat) = phenoData(gset[[1]])@data$title
mat = mat[rownames(getAnnotation(probe, what = "Locations")), ]
```


`probe` contains locations of probes and also information whether the CpG sites overlap
with SNPs. Here we remove probes that are on sex chromosomes and probes that overlap with SNPs.


```r
l = getAnnotation(probe, what = "Locations")$chr %in% paste0("chr", 1:22) &
    is.na(getAnnotation(probe, what = "SNPs.137CommonSingle")$Probe_rs)
mat = mat[l, ]
```

Extract the CpG annotations (i.e. CpG Islands, shores, seas).

```r
cgi_anno = getAnnotation(probe, "Islands.UCSC")$Relation_to_Island[l]
```

Extract the matrix for tumor samples. Also modify column names for the tumor
samples to be consistent with the phenotype data which we will read later.

```r
mat1 = as.matrix(mat[, grep("GBM", colnames(mat))])   # tumor samples
colnames(mat1) = gsub("GBM", "dkfz", colnames(mat1))
```

Read the annotations for samples. Phenotype data (`450K_annotation.txt`) is from [Sturm et al., 2012](http://dx.doi.org/10.1016/j.ccr.2012.08.024), supplementary table S1.

```r
phenotype = read.table("450K_annotation.txt", header = TRUE, sep = "\t", row.names = 1,
    check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
phenotype = phenotype[colnames(mat1), 1:2]
colnames(phenotype) = c("dkfz_subtype", "tcga_subtype")
```

Assign random values for `NA`.

```r
mat1[is.na(mat1)] = runif(sum(is.na(mat1)))
```

A subset of matrix which corresponds to the current CpG annotation.
The value of `cgi` can be `island`/`shore`/`sea`/`all`.

```r
cgi = "all" # island/shore/sea/all
cgi_regexp = cgi
if(cgi == "all") {
    cgi_regexp = ".*"
} else if(cgi %in% c("shelf", "sea")) {
    cgi_regexp = "shelf|sea"
}
mat1 = mat1[grepl(cgi_regexp, cgi_anno, ignore.case = TRUE), , drop = FALSE]
```

Define the colors for the annotations.

```r
anno_col = list(
    dkfz_subtype = structure(names = c("IDH", "K27", "G34", "RTK I PDGFRA", "Mesenchymal", "RTK II Classic"), brewer.pal(6, "Set1")),
    tcga_subtype = structure(names = c("G-CIMP+", "Cluster #2", "Cluster #3"), brewer.pal(3, "Set1"))
)
```

Perform the consensus partitioning:

```r
set.seed(123)
rl = run_all_consensus_partition_methods(
    mat1, 
    top_n = seq(min(2000, round(nrow(data) * 0.1)), min(10000, round(nrow(data) * 0.5)), length.out = 5),
    max_k = 10,
    scale_rows = FALSE, 
    anno = phenotype, 
    anno_col = anno_col, 
    mc.cores = 4)

saveRDS(rl, file = qq("GBM_450K_cgi_@{cgi}_subgroup.rds"))

cola_opt(group_diff = 0.1)
cola_report(rl, output_dir = qq("GBM_450K_cgi_@{cgi}_subgroup_cola_report"), mc.cores = 4)
```
