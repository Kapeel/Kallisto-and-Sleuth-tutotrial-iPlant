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

create a folder:
```{r}
mkdir kallisto_qaunt_human_RNAseq_output
cd kallisto_qaunt_human_RNAseq_output
```

Download the kallisto qaunt output to your home directory and unzip:
```{r}
wget http://de.iplantcollaborative.org/dl/d/9171661C-3746-4DD7-87D0-078B76D28CFC/kallisto_qaunt_human_RNAseq_output.zip
unzip kallisto_qaunt_human_RNAseq_output.zip
```


Start up RStudio and navigate to `kallisto_qaunt_human_RNAseq_output`  directory we be
working .

```{r}
setwd('~/Desktop/kallisto_qaunt_human_RNAseq_output')
```

First, let's install `sleuth` and `biomaRt`, a tool that we will use later for
getting the gene names:

```{r}
source('http://bioconductor.org/biocLite.R')
biocLite("devtools")
biocLite("pachterlab/sleuth")
biocLite("biomaRt")
```

You can execute a line in Rstudio using `ctrl + enter`.

Next, load sleuth:

```{r}
library('sleuth')
```

Though not required, I also suggest loading a package called `cowplot` which makes the `ggplot`
default much more aesthetically pleasing:

```{r}
install.packages('cowplot')
library('cowplot')
```

Let's also set the base working directory:

```{r}
base_dir <- '~/Desktop/kallisto_qaunt_human_RNAseq_output'
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
more hiseq_info.txt
run_accession experiment_accession spots condition sequencer sample
SRR493366 SRX145662 15117833 scramble hiseq A
SRR493367 SRX145663 17433672 scramble hiseq B
SRR493368 SRX145664 21830449 scramble hiseq C
SRR493369 SRX145665 17916102 HOXA1KD hiseq A
SRR493370 SRX145666 20141813 HOXA1KD hiseq B
SRR493371 SRX145667 23544153 HOXA1KD hiseq C

```


### locating kallisto output

Next, we have to tell `sleuth` where all the kallisto output is. If you've
organized your data like we did in the snake file, this is quite easy. Sleuth is
simply expecting a character array pointing to all the directories.

Next get the list of sample IDs with

```{r}
sample_id <- dir(file.path(base_dir,"results"))
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
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results", id, "kallisto"))
kal_dirs
```
```{r}
##                                                                SRR493366 
## "~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493366/kallisto" 
##                                                                SRR493367 
## "~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493367/kallisto" 
##                                                                SRR493368 
## "~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493368/kallisto" 
##                                                                SRR493369 
## "~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493369/kallisto" 
##                                                                SRR493370 
## "~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493370/kallisto" 
##                                                                SRR493371 
## "~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493371/kallisto"
```
The next step is to load an auxillary table that describes the experimental design and the relationship between the kallisto directories and the samples:

```{r}
s2c <- read.table(file.path(base_dir, "hiseq_info.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c
```
```{r}
##      sample condition
## 1 SRR493366  scramble
## 2 SRR493367  scramble
## 3 SRR493368  scramble
## 4 SRR493369   HOXA1KD
## 5 SRR493370   HOXA1KD
## 6 SRR493371   HOXA1KD
```
Now, we must enter the directories into a column in the table describing the experiment. This column must be labeled path, otherwise sleuth will throw an error. This is to ensure that the user can check which samples correspond to which kallisto runs

```{r}
s2c <- dplyr::mutate(s2c, path = kal_dirs)
```
The user should check whether or not the order is correct. In this case, the kallisto output is correctly matched with the sample identifiers.
```{r}
print(s2c)
##      sample condition
## 1 SRR493366  scramble
## 2 SRR493367  scramble
## 3 SRR493368  scramble
## 4 SRR493369   HOXA1KD
## 5 SRR493370   HOXA1KD
## 6 SRR493371   HOXA1KD
##                                                                     path
## 1 ~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493366/kallisto
## 2 ~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493367/kallisto
## 3 ~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493368/kallisto
## 4 ~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493369/kallisto
## 5 ~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493370/kallisto
## 6 ~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493371/kallisto
```

### fitting the model

Now the “sleuth object” can be constructed. This requires three commands that
(1) load the kallisto processed data into the object (2) estimate parameters for
the sleuth response error measurement model and (3) perform differential analyis
(testing). On a laptop the three steps should take about 2 minutes altogether.

First type

```{r}
so <- sleuth_prep(s2c, ~ condition)
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
### getting gene names

At this point the sleuth object constructed from the kallisto runs has information about the data, the experimental design, the kallisto estimates, the model fit, and the testing. In other words it contains the entire analysis of the data. There is, however, one piece of information that can be useful to add in, but that is optional. In reading the kallisto output sleuth has no information about genes, but this can be added allowing for searching and analysis by gene instead of transcript.

Since the example was constructed with the ENSEMBL human transcriptome, we will add gene names from ENSEMBL using biomaRt (there are other ways to do this as well):

First, install biomaRt with
```{r}
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
```
Then collect gene names with 
```{r}
# get the gene names using biomaRt
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
```
and add them into the sleuth table with
```{r}
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)
```
```{r}
## reading in kallisto results
## ......
## normalizing est_counts
## 50844 targets passed the filter
## normalizing tpm
## merging in metadata
## normalizing bootstrap samples
## summarizing bootstraps
```
```{r}
so <- sleuth_fit(so)
## shrinkage estimation
## computing variance of betas
so <- sleuth_wt(so, which_beta = 'conditionscramble')
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
