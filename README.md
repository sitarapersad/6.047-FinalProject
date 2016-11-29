# 6.047-FinalProject

## Timeline 

Edit as you see fit.

### November 29
1. Install LDSC tools for score regression.

Reasonably chill process using https://github.com/bulik/ldsc. There are some useful tutorials there, as well. Our primary interest is in estimating genetic correlation, as per the following: https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation.

2. Download SNP data from the website

Done, uploading to git. Be careful in the future, since these data files are so big, only git add . subfolders which are modified, to save time.

### December 2

1. Write scripts to compute LD gene heritability using the aforementioned tools

I'm guessing a bash script will be more appropriate for this, although we'll have to see.

Apparently the data looks like this initially:

snpid hg18chr bp a1 a2 or se pval info ngt CEUaf

rs3131972   1   742584  A   G   1   0.0966  0.9991  0.702   0   0.16055

We're interested in knowing chromosome (col[1]) and base pair (col[2]) for analysing chromosome by chromosome, and then for binary partitioning.


### December 3

1. Write code to implement binary search of a given chromosome.

2. Test code on one chromosome for two selected diseases.

### December 4 -  6

1. Run code on all the combinations to uncover all the regions of interest.

We need a way to verify these regions compared with known regions of interest.

## Random Notes for Now

IMPT: How do we force LDSC to compute the LD score over a particular chromosome region?
We can do this by parsing the text file. This might be the job for a particular python script. Perhaps this is somehow build into LDSC.
Parsing via text files will make this a slightly time-consuming process to create and destroy these large text files repeatedly.


How do we parse the log file of outputs to obtain the genetic correlation?

IMPT: How do we select a threshold beyond which a region no longer explains significant genetic correlation?

How do we verify that regions we find are good regions?
Is there any previous work that identifies genes/ regions that are biologically interesting?



## Writeup
https://www.overleaf.com/7238811cscdzymnmskg

## Citations
Neale et al., Meta-analysis of genome-wide association studies of attention-
deficit/hyperactivity disorder. J. Am. Acad. Child Adolesc. Psychiatry 49, 884
(Sep, 2010).

Psychiatric GWAS Consortium Bipolar Disorder Working Group, Large-scale
genome-wide association analysis of bipolar disorder identifies a new
susceptibility locus near ODZ4. Nature Genetics 43, 977 (2011).

Major Depressive Disorder Working Group of the PGC, A mega-analysis of
genome-wide association studies for major depressive disorder. Molecular
Psychiatry. (In press).

Schizophrenia PGC, Genome-wide association study of schizophrenia identifies
five novel loci. Nature Genetics 43, 969 (2011).
