# RCSJ
Simulation of Resistively and Capacitively Shunted Josephson junction array.

## Get started:

    mkdir snspd
    cd snspd
    git clone https://github.com/jlidmar/RCSJ.git src

    cd src
    make
    make install

    cd ../run

    nice ./xiv

or

    sbatch -t 60 -a 1 slurm r xxiv

## TODO:

- Normal juncion not fully consistent as it is now.  vgap is modified not not further down.