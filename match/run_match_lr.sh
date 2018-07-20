set -e
for i in `seq 1 22`
do
echo "CEU Chromosome $i"
python -u scripts/neanderthal_match.py ../1000G_predict_LR/CEU/CEU.chr$i-prediction.txt data/neanderthal/neander.chr-$i.freq data/vcf/CEU.chr$i.vcf.gz $i > results_LR/CEU-CHR_$i-match-stat_LR.txt
echo "YRI Chromosome $i"
python -u scripts/neanderthal_match.py ../1000G_predict_LR/YRI/YRI.chr$i-prediction.txt data/neanderthal/neander.chr-$i.freq data/vcf/YRI.chr$i.vcf.gz $i > results_LR/YRI-CHR_$i-match-stat_LR.txt
done
