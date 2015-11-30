# sleuth

In the wild, a sleuth is a pack of bears. In the context of RNA-Seq, a sleuth is
a pack of kallisto. The job of sleuth is to perform aggregate analysis of many
related samples at once. Currently, the main role of sleuth is to do
differential expression analysis at the transcript level. Sleuth differs from most other differential
expression tools by modeling the technical variance due to the transcript abundance
estimation along with the biological variability between samples. Most methods
only model the biological variability.

To explain how to use __sleuth__ we provide an example based on the data in the "Cuffdiff2 paper":

* [Differential analysis of gene regulation at transcript resolution with RNA-seq](http://www.nature.com/nbt/journal/v31/n1/full/nbt.2450.html)	by Cole Trapnell, David G Henderickson, Martin Savageau, Loyal Goff, John L Rinn and Lior Pachter, Nature Biotechnology __31__, 46--53 (2013).

The human fibroblast RNA-Seq data for the paper is available on GEO at accession [GSE37704](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37704). The samples to be analyzed are the six samples LFB_scramble_hiseq_repA, LFB_scramble_hiseq_repB, LFB_scramble_hiseq_repC, LFB_HOXA1KD_hiseq_repA, LFB_HOXA1KD_hiseq_repA, and LFB_HOXA1KD_hiseq_repC. These are three biological replicates in each of two conditions (scramble and HoxA1 knockdown) that will be compared with __sleuth__.

## preliminaries

Download the kallisto qaunt output to your home directory:
```{r}
wget http://de.iplantcollaborative.org/dl/d/9171661C-3746-4DD7-87D0-078B76D28CFC/kallisto_qaunt_human_RNAseq_output.zip

```
Start up RStudio and navigate to `kallisto_qaunt_human_RNAseq_output` subdirectory in the directory we've been
working in.

```{r}
setwd('~/analysis/bears_iplant/R')
```

First, let's install `sleuth` and `biomaRt`, a tool that we will use later for
getting the gene names:

```{r}
source('http://bioconductor.org/biocLite.R')
biocLite("devtools")
biocLite("pachterlab/sleuth")
biocLite("biomaRt")
```

Open a new file if you would like to type the commands and add them as we go
along, or you can simply open `analysis.R` and follow along. You can execute a
line in Rstudio using `ctrl + enter`.

Next, load sleuth:

```{r}
library('sleuth')
```

Though not required, I also suggest loading a package called `cowplot` which makes the `ggplot`
default much more aesthetically pleasing:

```{r}
# install.packages('cowplot')
# if it isn't installed
library('cowplot')
```

Let's also set the base working directory:

```{r}
base_dir <- '..'
```

From here on, all the commands will be in R unless otherwise specified.

## preparing your data

The main requirements of sleuth are:

- sample to covariates mapping
- output from kallisto

### sample to covariate mapping

The sample to covariate mapping is simply a table that describes the experiment. possibly the most
challenging part of sleuth is organizing your data into a way that it can be
read and easily. The only real requirement is that there is at least one column
labeled 'sample'. The remaining columns are important as they describe your
experiment, but the column names can be pretty much any valid string.

Our data is pretty simple in that there is only one covariate here: the
experimental condition.

This is what the file looks like (from the terminal):

```{sh}
cat metadata/sample_info.tsv
sample  condition
SRR493366       scramble
SRR493367       scramble
SRR493368       scramble
SRR493369       HOXA1KD
SRR493370       HOXA1KD
SRR493371       HOXA1KD
```

Let's load this file in R:

```{r}
s2c <- read.table(file.path(base_dir, 'metadata', 'sample_info.tsv'),
  header = TRUE, stringsAsFactors = FALSE)
```

### locating kallisto output

Next, we have to tell `sleuth` where all the kallisto output is. If you've
organized your data like we did in the snake file, this is quite easy. Sleuth is
simply expecting a character array pointing to all the directories.

Next get the list of sample IDs with

```{r}
sample_id <- dir(file.path(base_dir, "results", "paired"))
```

The result can be displayed by typing
```{r}
sample_id
## [1] "SRR493366" "SRR493367" "SRR493368" "SRR493369" "SRR493370" "SRR493371"
```

In the box above, lines beginning with ## show the output of the command (in
what follows we include the output that should appear with each command).

A list of paths to the kallisto results indexed by the sample IDs is collated with

```{r}
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results",
  "paired", id, "kallisto"))
```

### getting gene names

Since the gene names are not automatically in the annotation, we need to get
them from elsewhere. `biomaRt` provide a good way of getting the mappings for
Ensembl annotations.

```{r}
# get the gene names using biomaRt
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
```

### fitting the model

Now the “sleuth object” can be constructed. This requires three commands that
(1) load the kallisto processed data into the object (2) estimate parameters for
the sleuth response error measurement model and (3) perform differential analyis
(testing). On a laptop the three steps should take about 2 minutes altogether.

First type

```{r}
so <- sleuth_prep(kal_dirs, s2c, ~ condition, target_mapping = t2g)
## reading in kallisto results
## ......
## normalizing est_counts
## 42193 targets passed the filter
## normalizing tpm
## normalizing bootstrap samples
```

then

```{r}
so <- sleuth_fit(so)
## summarizing bootstraps
## fitting measurement error models
## shrinkage estimation
## computing variance of betas
```

and finally

```{r}
so <- sleuth_test(so, which_beta = 'conditionscramble')
```

In general, one can see the possible tests that could be performed using the
`which_beta` parameter in `sleuth_test` and examining the coefficients:

```{r}
models(so)
## [  full  ]
## formula:  ~condition
## coefficients:
##  (Intercept)
##      conditionscramble
## tests:
##  conditionscramble
```

### interactive analysis

Sleuth provides many different ways visualize your data. Most visualizations are
prefixed by `plot_`. While this is true, we think the best way to analyze your
data is using sleuth live. Sleuth live gives you an interactive visualization
along with all the differential expression analysis together. You can execute
sleuth live with:

```{r}
sleuth_live(so)
```

Let's chat about what sort of things to look out for.
