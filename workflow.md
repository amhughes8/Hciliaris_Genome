https://github.com/amhughes8/Hciliaris_Genome
https://github.com/remygatins/Holacanthus_passer-ONT-Illumina-Genome-Assembly

ssh hughes.annab@login.discovery.neu.edu
cd /work/gatins/hci_genome

# get out of login node
srun --partition=short --nodes=1 --cpus-per-task=1 --pty /bin/bash
module load anaconda3/2022.05 discovery

# creating environment -- took FOREVER, don't run again! just needed to set up environment
conda create --prefix=/work/gatins/hci_genome/env python=3.11 anaconda

# activate environment:
source activate /work/gatins/hci_genome/env

#setting up environments:

##Porechop v0.2.4
conda install bioconda::porechop
### this install was corrupted somehow, so I had to install from source
git clone https://github.com/rrwick/Porechop.git
cd Porechop
python3 setup.py install
porechop -h
### also didn't work because I received a permissions error due to gcc?...
module load cmake/3.27.9
module load gcc/9.2.0
git clone https://github.com/rrwick/Porechop.git
cd Porechop
make
./porechop-runner.py -h
### this worked! I really don't get why UGH but it's installed so yay

##seqkit v2.10.0
conda install bioconda::seqkit

##Hifiasm-0.25.0-r726
git clone https://github.com/chhylp123/hifiasm
cd hifiasm && make
### this I installed before reinstalling Porechop (above) and it had worked without module loading cmake which is why I'm confused about needing that for porechop

##Nanostat
### so annoying to download but eventually worked
module load gcc/9.2.0
pip install nanoplot

##Nanoplot
pip install nanoplot

##Flye
git clone https://github.com/fenderglass/Flye
cd Flye
python setup.py install

##Medaka
### the user must provide several binaries: samtools, minimap2, tabix, and bgzip
module load samtools/1.9
module load minimap2/2.26
### need HTSlib to download bgzip
git clone https://github.com/samtools/htslib.git
git submodule update --init --recursive #ensure all submodules are up to date
module load cmake/3.27.9
module load gcc/9.2.0
make; make install #I kept getting permission errors when trying to install, eventually I had to alter my .bashrc file to look as follows to direct the installation only within directories i had downloaded myself:
export PATH="/work/gatins/hci_genome/htslib:$PATH"
export LD_LIBRARY_PATH="/work/gatins/hci_genome/htslib:$LD_LIBRARY_PATH"
export PKG_CONFIG_PATH="/work/gatins/hci_genome/htslib/config.mk.in:$PKG_CONFIG_PATH"

##Busco
### Busco gave me errors every time I tried to download it, so I created a new environment for it and it worked
conda create --prefix=/work/gatins/hci_genome/busco
source activate /work/gatins/hci_genome/busco
conda install -c conda-forge -c bioconda busco=5.8.2

export LD_LIBRARY_PATH="/work/gatins/hci_genome/htslib/lib:$LD_LIBRARY_PATH"

export PKG_CONFIG_PATH="work/gatins/hci_genome/htslib/lib/pkgconfig:$PKG_CONFIG_PATH"