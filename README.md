# ArchIE
ArchIE (ARCHaic Introgression Explorer) is a method for reference-free inference of archaic local ancestry. It uses coalescent simulations to train a model to distinguish between archaic and non-archaic windows. The process is as follows:

1. run `bash create_training.sh` to simulate training data and calculate features (from the directory)
2. run `bash create_test.sh` to simulate test data (from the directory)
3. run `Rscript train.R` to train a logistic regression model (from the directory)
4. run `python data/calc_stats_window_data.py -s data.snp -i data.ind -a data.geno -r ref.geno -c 1 -b $start -e $end -w 50000 -z 50000 | gzip > output.gz`

The last step will calculate features in 50KB windows with a 50KB step size (i.e. nonoverlapping) from `$start` to `$end` on chromosome 1.

You can find out more about the options for `calc_stats_window_data.py` using the -h flag:

  python data/calc_stats_window_data.py -h
  usage: calc_stats_window_data.py [-h] -s SNP -i IND -a ADMIX [-r REFERENCE] -c
                                   CHROM -b BEGIN -e END -w WINDOW -z STEP

  Script to calculate stats for ArchIE. Prints stats to stdout.

  optional arguments:
    -h, --help            show this help message and exit
    -s SNP, --snp SNP     .snp file with positions (eigenstrat format)
    -i IND, --ind IND     .ind file of IDS (eigenstrat format)
    -a ADMIX, --admix ADMIX
                          .geno file for admixed/test population (eigenstrat
                          format)
    -r REFERENCE, --reference REFERENCE
                          .geno file for reference population (eigenstrat
                          format).
    -c CHROM, --chrom CHROM
                          Chromosome
    -b BEGIN, --begin BEGIN
                          Beginning position
    -e END, --end END     End position
    -w WINDOW, --window WINDOW
                          Window length (eg. 50000)
    -z STEP, --step STEP  Step size for the window (eg. 10000)

# Requirements

ArchIE requires the following dependencies:

- Python 3
- Numpy (1.13.0)
- Scipy (0.19.0)
- Scikit Learn (0.19.1)
- R (3.4.0)

Note that other versions may work -- these are the ones I used.

# Data formats

ArchIE uses the eigenstrat format, which contains `.snp` (list of polymorphic positions), `.geno` (genotypes), and `.ind` (individual IDs) files. ArchIE expects phased genotypes. You can convert from VCF format to eigenstrat using Iain Mathieson's [genetic data converter](https://github.com/mathii/gdc). To calculate features on your data, use the `calc_stats_window_data.py` script in the data folder.


# Applying to your data

When running ArchIE, you can specify the number of haploid genomes to simulate (n), the window size, or the demographic history (in `ms.sh`). You must simulate the same number of haploid genomes as are in your data set. E.g if I have 15 diploid individuals with phased genomes, I need to simulate a dataset with 30 haploid genomes.  

# Features

We use several features to predict local ancestry. When the number of haplotypes simulated is 100 (the default) this is what the columns describe:

- 1-100 The individual frequency spectrum (sample size dependent)
- 101-201 - the distribution of distances between haplotypes (sample size dependent)
- 202 - Mean of the distribution of distances between haplotypes
- 203 - Variance of the distribution of distances between haplotypes
- 204 - Skew of the distribution of distances between haplotypes
- 205 - Kurtosis of the distribution of distances between haplotypes
- 206 - Minimum distance to the reference population
- 207 - [S*](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0020105)
- 208 - Number of SNPs private to the target population
- 209-211 one hot encoded label for [archaic, not archaic, in-between]. 'In-between' is chosen if haplotypes fall between 30-70% archaic. This can be modified in the `calc_stats_ms.py` script
- 212 - Proportion of the haplotype that is archaic ranges [0-1]

The order of the haplotypes remains the same if the sample size is changed, but the column numbers will change.

# Note about `ms`

We use a customized version of `ms` that allows us to track the ancestry of each SNP. Check around line 218 in ms.c for the modifications. `ms` was compiled on Mac OSX, if you run into problems you can recompile it with `gcc -o ms ms.c streec.c rand1.c -lm`
