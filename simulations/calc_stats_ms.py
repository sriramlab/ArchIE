import argparse
import sys
import math
import csv
import copy
import scipy

import numpy as np
import scipy.stats as sp
import sklearn.cluster
import sklearn.metrics.pairwise

import warnings
warnings.filterwarnings("ignore")

def min_dist_ref(genotypes, ref_genotypes, focal_idx):
    f = np.array(genotypes)[:,focal_idx]
    t_r = np.transpose(np.array(ref_genotypes))
    d = []
    for r in t_r:
        r_np = np.array(r)
        arr = np.array([f, r_np])
        d.append(np.max(sklearn.metrics.pairwise.pairwise_distances(arr)))
    return(np.min(d))

def N_ton(genotypes, n_samples, focal_idx):
    N_ton_vec = [0] * (n_samples + 1)
    for g in genotypes:
        freq = g.count(1)
        if g[focal_idx] == 1:
            N_ton_vec[freq] = N_ton_vec[freq] + 1
    return(N_ton_vec)

def distance_vector(genotypes, focal_idx):
    dist = sklearn.metrics.pairwise.pairwise_distances(np.transpose(genotypes))
    focal_dist = dist[focal_idx]
    focal_dist.sort()
    m = np.mean(focal_dist)
    var = np.var(focal_dist)
    skew = sp.skew(focal_dist)
    kurtosis = sp.kurtosis(focal_dist)
    return(np.ndarray.tolist(focal_dist) + [m] + [var] + [skew] + [kurtosis])


def score(twosites, twopos):
    length_factor = np.absolute(np.diff(twopos))
    if length_factor < 10 or length_factor > 50000:
        return(-np.inf)

    nderiv = np.sum(twosites, 0)
    n1dot = nderiv[0]
    ndot1 = nderiv[1]
    if n1dot <= 1 or ndot1 <=1: # no singletons
        return(-np.inf)

    d = []
    for chrom in twosites:
        d.append(np.abs(chrom[0]-chrom[1]))
    distance = np.sum(d)

    if distance < 1:
        return(length_factor[0]+5000)
    if distance <= 5:
        return(-10000)
    return(-np.inf)

def s_star_ind(s, pos, focal_idx):
    # remove SNPs from the matrix where focal_idx is ancestral
    focal_hap = s[focal_idx]
    to_remove = [] # list of indexes
    for idx,site in enumerate(focal_hap):
        if site == 0:
            to_remove.append(idx)

    new_s = np.delete(s, to_remove, axis=1)
    s = new_s
    nsites = s.shape[1]

    if(nsites == 0):
        return(0)
    best_score = np.zeros((nsites, nsites))
    best_at_this_point = ["NA"] * nsites
    best_at_this_point[0] = 0
    for i in range(1, nsites):
        best_at_this_point[i] = 0
        for j in range(0, i):
            best_score[i,j] = np.maximum(0, best_at_this_point[j] + score(s[:,[j,i]], pos[[j,i],]))
            if best_score[i,j] > best_at_this_point[i]:
                best_at_this_point[i] = best_score[i,j]

    return(np.amax(best_score))

def num_private(clean_geno, focal_idx):
    """"input is a genotype matrix with (rows=samples), (columns=snps) and SNPs in reference removed"""
    focal_hap = clean_geno[focal_idx]
    return(np.sum(focal_hap))


def parse_snp(snp):
    """Returns a list of SNPs"""
    pos = []
    with open(snp, 'r') as f:
        for line in f:
            pos.append(int(line.split("\t")[3]))
    return(pos)

def parse_geno(geno_f):
    """Returns a list of lists of genotypes (2D array) (sample x snp)"""
    geno = []
    with open(geno_f, 'r') as f:
            for line in f:
                ll = list(line)[:-1]
                li = [int(i) for i in ll]
                geno.append(li)
    return(geno)

def parse_anc(anc_f, n_samples):
    """Returns a list n_samples long of the number of archaic bases per sample"""
    anc = []
    with open(anc_f, 'r') as f:
        for line in f:
            anc.append(list(line)[:-1])
    return(anc)

def label(bases, snp_list, n_sites, arch_thresh, not_arch_thresh):
    archaic_anc = 0
    prev_pos = 0
    previous_anc = 0
    for idx,base in enumerate(bases):
        s = snp_list[idx]
        if base == "1" and previous_anc == "1":
            archaic_anc = archaic_anc + (s-previous_pos)
        previous_pos = s
        previous_anc = base
    prop = archaic_anc/n_sites
    if prop > arch_thresh:
        return([1,0,0, prop])
    elif prop < not_arch_thresh:
        return([0,1,0, prop])
    else:
        return([0,0,1, prop])


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script to calculate stats for AINN (archaic introgression neural network). Prints stats to stdout.")
    parser.add_argument("-s", "--snp", action="store", required=True, help=".snp file with positions (eigenstrat format)")
    parser.add_argument("-a", "--admix", action="store", required=True, help=".geno file for admixed/test population (eigenstrat format)")
    parser.add_argument("-r", "--reference", action="store", required=False, help=".geno file for reference population (eigenstrat format).")
    parser.add_argument("--anc", action="store", required=True, help=".anc file to output a label ")
    parser.add_argument("-n", "--nsites", action="store", required=True, help="Number of sites in the data including homozygous ancestral")
    args = parser.parse_args()

    try:
        mutation_positions = parse_snp(args.snp)
    except FileNotFoundError:
        print("Error reading SNP file. Check file path.", file=sys.stderr)
        sys.exit(1)

    try:
        genotypes = parse_geno(args.admix)
    except FileNotFoundError:
        print("Error reading admix .geno file. Check file path.", file=sys.stderr)
        sys.exit(1)

    try:
        ref_geno = parse_geno(args.reference)
    except FileNotFoundError:
        print("Error reading reference .geno file. Check file path.", file=sys.stderr)
        sys.exit(1)

    n_samples = len(genotypes[0])

    try:
        arch = parse_anc(args.anc, n_samples)
    except FileNotFoundError:
        print("Error reading anc file. Check file path.", file=sys.stderr)
        sys.exit(1)

    n_sites = int(args.nsites)

    ## set up S* stuff -- remove mutations found in reference set
    t_ref = list(map(list, zip(*ref_geno)))
    t_geno =  list(map(list, zip(*genotypes)))
    pos_to_remove = set() # contains indexes to remove
    s_star_haps = []
    for idx, hap in enumerate(t_ref):
        for jdx, site in enumerate(hap):
            if site == 1:
                pos_to_remove.add(jdx)

    for idx, hap in enumerate(t_geno):
        s_star_haps.append([v for i, v in enumerate(hap) if i not in pos_to_remove])

    #individual level
    for focal_idx in range(0, n_samples):
        calc_N_ton = N_ton(genotypes, n_samples, focal_idx)
        dist = distance_vector(genotypes, focal_idx)
        min_d = [min_dist_ref(genotypes, ref_geno, focal_idx)]
        ss = [s_star_ind(np.array(s_star_haps), np.array(mutation_positions), focal_idx)]
        n_priv = [num_private(np.array(s_star_haps), focal_idx)]
        focal_arch = [row[focal_idx] for row in arch ]
        lab = label(focal_arch, mutation_positions, n_sites, 0.7, 0.3)

        output = calc_N_ton + dist + min_d + ss + n_priv + lab
        print(*output, sep="\t") # print stats to standard out
