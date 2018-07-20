set -e
L=$1
F=$2
T=`python -c "print(5e-4*${L})"`
R=`python -c "print(4e-4*${L})"`
A=`python -c "print(1-${F})"`
../msmodified/ms 202 1 -T -t ${T} -r ${R} ${L} -I 4 100 100 1 1 g  -en 0 1 1  -es 0.05 1 ${A} -ej 0.05 5 3  -ej 0.0625 2 1 -en 0.15 3 0.01 -en 0.153 3 1 -ej 0.175 4 3 -ej 0.3 3 1| tail -n1
