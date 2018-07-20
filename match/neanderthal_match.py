from sys import argv
import sys
import datetime
import gzip
from bisect import bisect, bisect_left

calls = argv[1]
freq = argv[2]
vcf = argv[3]
CHROM = argv[4]

nean_dict = {}
nean_pos = []
nean_geno = []
print("["+str(datetime.datetime.now().time())+'] Parsing Neanderthal map...', file=sys.stderr)
with open(freq, 'r') as f:
    for line in f:
        c = int(float(line.split(" ")[1].rstrip())*2) # convert freq to counts, missing is now -2
        p = int(line.split(" ")[0].split(":")[1])
        nean_dict[p] = c
        nean_pos.append(p)
        nean_geno.append(c)

samples = []
human_pos = []
human_geno = []
print("["+str(datetime.datetime.now().time())+'] Parsing human VCF...', file=sys.stderr)
with gzip.open(vcf, 'rt') as f:
    for line in f:
        if line[:2]=="##":				  # Comment line
            continue
        elif line[:6]=="#CHROM":
            samples=line.split()[9:]
        else:
            b = line.split()
            if "," in b[4]:
                continue
            if len(b[3])!=1 or len(b[4])!=1:
                continue
            else:
                g1 = b[9:]
                g2 = [i.split("|") for i in g1]
                g3 = [int(item) for sublist in g2 for item in sublist]
                human_geno.append(g3)
                human_pos.append(int(b[1]))

t_human_geno = list(map(list, zip(*human_geno)))
print("["+str(datetime.datetime.now().time())+'] Calculating statistics...', file=sys.stderr)
with open(calls, 'r') as f:
    for j,line in enumerate(f):
        b = line.split()
        sample = b[2]
        sample_idx = samples.index(sample)
        if j%2==0:
            sample_geno = t_human_geno[sample_idx*2]
        else:
            sample_geno = t_human_geno[sample_idx*2+1]

        start = int(b[4])
        end = int(b[5])
        try:
            h_start_idx = bisect(human_pos, start)
            h_end_idx = bisect_left(human_pos, end)
            n_start_idx = bisect(nean_pos, start)
            n_end_idx = bisect_left(nean_pos, end)
        except StopIteration:
            print(*[start, end, sample, b[0], b[1],0, 0, 0, 0], sep="\t")
            continue
        region_h_pos = human_pos[h_start_idx:h_end_idx+1]
        region_h_geno = sample_geno[h_start_idx:h_end_idx+1]
        region_n_pos = nean_pos[n_start_idx:n_end_idx+1]
        region_n_geno = nean_geno[n_start_idx:n_end_idx+1]


        num_shared = 0
        for i,p in enumerate(region_h_pos):
            try:
                n_g = nean_dict[p]
            except KeyError:
                continue #if there is no neanderthal variant, should skip
            h_g = region_h_geno[i]
            if n_g >= 1 and h_g == 1: #nean is unphased so can be 2
                num_shared += 1

        num_nean = region_n_geno.count(1) + region_n_geno.count(2) - num_shared
        num_human = region_h_geno.count(1)
        try:
            statistic = num_shared/(num_nean+num_human)
        except ZeroDivisionError:
            statistic = 0
        print(*[CHROM, start, end, sample, b[0], b[1], num_nean, num_human, num_shared, statistic], sep="\t")
