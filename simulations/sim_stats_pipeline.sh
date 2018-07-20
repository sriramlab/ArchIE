### Test the stat distributions for admix/not
set -e
SUFF=$RANDOM

f=0.02
for i in `seq 1 10000`;
do
  echo "Simulation " ${i}
  bash ms.sh 50000 ${f}
  python calc_stats_ms.py -s out.snp -a out.ADMIXED.geno -r out.1.geno --anc out.ADMIXED.anc -n 50000 > sims/${i}.txt
done

cat sims/* > test_stats.txt
rm sims/*
