# ArchIE
ArchIE (ARCHaic Introgression Explorer) is a method for reference-free inference of archaic local ancestry. It uses coalescent simulations to train a model to distinguish between archaic and non-archaic windows. The process is as follows:

1. run `bash create_training.sh` to simulate training data and calculate features (from the directory)
2. run `bash create_test.sh` to simulate test data (from the directory)
3. run `Rscript train.R` to train a logistic regression model (from the directory)

# Data formats

ArchIE uses the eigenstrat format, which contains `.snp` (list of polymorphic positions), `.geno` (genotypes), and `.ind` (individual IDs) files. ArchIE expects phased genotypes. You can convert from VCF format to eigenstrat using Iain Mathieson's [genetic data converter](https://github.com/mathii/gdc). To calculate features on your data, use the `calc_stats_window_data.py` script in the data folder.


# Applying to your data

When running ArchIE, you can specify the number of haploid genomes to simulate ($n$), the window size, or the demographic history (in `ms.sh`). You must simulate the same number of haploid genomes as are in your data set. E.g if I have 15 diploid individuals with phased genomes, I need to simulate a dataset with 30 haploid genomes.  

# Features

# Note about `ms`

We use a customized version of `ms` that allows us to track the ancestry of each base. Check around line 218 in ms.c for the modifications. `ms` was compiled on Mac OSX, if you run into problems you can recompile it with `gcc -o ms ms.c streec.c rand1.c -lm`
