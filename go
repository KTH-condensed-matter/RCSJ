#!/bin/sh

for RC in 0.1 0.5 1 2 5 10; do

dir=RC=$RC
mkdir -p $dir
cd $dir

ln -s ../rcsj .

echo RC=$RC

for T in 0.01 0.1 0.5 ; do

dir=T=$T
mkdir -p $dir
cd $dir

echo $dir
sed -e 's/%RC%/'$RC'/g' -e 's/%T%/'$T'/g' ../../xrcsj.tmpl > xrcsj
chmod u+x xrcsj

# Start the job:
# submit job
# qsub -t 1 xrcsj
# ./xrcsj &

cd ..
done

cd ..
done
