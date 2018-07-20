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

def window_kmeans(mutations, window_size):
    mutations = mutations
    win_size = window_size

    num_windows = int(n_sites)/win_size
    win_end = np.linspace(mutation_positions[0], mutation_positions[0]+int(n_sites), num_windows, endpoint=True)
    win_k1 = np.zeros(num_windows)
    win_k2 = np.zeros(num_windows)

    for idx,window in enumerate(win_end):
        current_win = []
        for jdx,mutation in enumerate(mutations):
            if mutation < win_end[idx] and mutation > win_end[idx-1]:
                current_win.append(genotypes[jdx])
        try:
            t_genotypes = list(map(list, zip(*current_win)))
            win_k1[idx] = sklearn.cluster.KMeans(n_clusters=1).fit(t_genotypes).inertia_
            win_k2[idx] = sklearn.cluster.KMeans(n_clusters=2).fit(t_genotypes).inertia_
        except Exception as e:
            #print(str(e), t_genotypes, mutation)
            continue

    return(list(win_k1) + list(win_k2))

def cSFS(genotypes, n_samples, focal_idx):
    focal = [0, 0]
    others = [0] * (2*n_samples)
    geno = copy.deepcopy(genotypes)

    for row in geno:
        freq = row.count(1)
        focal_row = row[focal_idx]
        others_row = row
        del others_row[focal_idx]

        focal_freq = focal_row
        try:
            focal[focal_freq] = focal[focal_freq] + 1
        except IndexError:
            pass

        others_freq = others_row.count(1)
        others[others_freq] = others[others_freq] + 1

    #sFocal = [i/sum(focal) for i in focal]
    #sOthers = [i/sum(others) for i in others]
    return(focal + others)

def SFS(genotypes, n_samples):
    SFS = [0] * (2*n_samples)
    geno = copy.deepcopy(genotypes)
    for row in geno:
        freq = row.count(1)
        try:
            SFS[freq] = SFS[freq] + 1
        except IndexError:
            pass #fixed in outgroup sample, but not in anc ref

    sSFS = [i/sum(SFS) for i in SFS]
    sSFS.pop(0) # remove zeroth entry
    return(sSFS)

def pi(genotypes, n_sites, n_samples):
    site_pi = []
    for site in genotypes:
        p = site.count(0)
        q = site.count(1)
        site_pi.append(2*p*q/n_samples)
    return(sum(site_pi)/n_sites)

def S(genotypes, n_sites):
    return(len(genotypes)/n_sites)

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

def window_pi(mutations, window_size):
    mutations = mutations
    win_size = window_size

    num_windows = int(n_sites)/win_size
    win_end = np.linspace(mutation_positions[0], mutation_positions[0]+int(n_sites), num_windows, endpoint=True)
    win_pi = np.zeros(num_windows)

    for idx,window in enumerate(win_end):
        for jdx,mutation in enumerate(mutations):
            if mutation < win_end[idx] and mutation > win_end[idx-1]:
                p = genotypes[jdx].count(0)/len(genotypes[jdx])
                q = genotypes[jdx].count(1)/len(genotypes[jdx])
                win_pi[idx] += 2*p*q

    return([round(pi/win_size, 3) for pi in win_pi])

def LOO_kmeans(genotypes, focal_idx):
    # leave one out kmeans (remove focal individual and see which K it belongs to)
    geno = copy.deepcopy(genotypes)
    t_genotypes = list(map(list, zip(*geno)))

    t_genotypes2 = list(map(list, zip(*genotypes)))
    focal = t_genotypes2[focal_idx]

    k1 = sklearn.cluster.KMeans(n_clusters=1).fit(t_genotypes)
    k2 = sklearn.cluster.KMeans(n_clusters=2).fit(t_genotypes)
    _all = scipy.stats.mode(k2.predict(t_genotypes2))[0]
    foc = k2.predict(focal)[0]
    if foc == _all:
        return(1) # 1 if same cluster as most genotypes
    else:
        return(0) # 0 if different from most genotypes
    return([k2.predict(t_genotypes2),k2.predict(focal)])

def kmeans(genotypes):
    t_genotypes = list(map(list, zip(*genotypes)))
    k1 = sklearn.cluster.KMeans(n_clusters=1).fit(t_genotypes)
    k2 = sklearn.cluster.KMeans(n_clusters=2).fit(t_genotypes)
    return([k1.inertia_, k2.inertia_])

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
    #return([prop])

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

    # calc_pi = [pi(genotypes, n_sites, n_samples)]
    # calc_S = [S(genotypes, n_sites)]
    # calc_kmeans = kmeans(genotypes)
    # calc_SFS = SFS(genotypes, n_samples)
    # w_pi_1000 = window_pi(mutation_positions, 1000)
    # w_pi_5000 = window_pi(mutation_positions, 5000)
    # w_pi_10000 = window_pi(mutation_positions, 10000)
    #
    # w_k_1000 = window_kmeans(mutation_positions, 1000)
    # w_k_5000 = window_kmeans(mutation_positions, 5000)
    # w_k_10000 = window_kmeans(mutation_positions, 10000)

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
        # calc_cSFS = cSFS(genotypes, n_sites, focal_idx)
        # calc_loo_kmeans = [LOO_kmeans(genotypes, focal_idx)]
        calc_N_ton = N_ton(genotypes, n_samples, focal_idx)
        dist = distance_vector(genotypes, focal_idx)
        min_d = [min_dist_ref(genotypes, ref_geno, focal_idx)]
        ss = [s_star_ind(np.array(s_star_haps), np.array(mutation_positions), focal_idx)]
        n_priv = [num_private(np.array(s_star_haps), focal_idx)]
        focal_arch = [row[focal_idx] for row in arch ]
        lab = label(focal_arch, mutation_positions, n_sites, 0.7, 0.3)
        # output=lab
        #output = calc_SFS + calc_cSFS + calc_pi + calc_S + calc_kmeans + calc_loo_kmeans + w_pi_1000 + w_pi_5000 + w_pi_10000 + w_k_1000 + w_k_5000 + w_k_10000 + calc_N_ton + dist + min_d + ss + lab
        output = calc_N_ton + dist + min_d + ss + n_priv + lab
        print(*output, sep="\t") # print stats to standard out
