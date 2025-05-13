ssh hughes.annab@login.discovery.neu.edu
cd /work/gatins/hci_genome

# get out of login node
srun --partition=short --nodes=1 --cpus-per-task=5 --time=48:00:00 --pty /bin/bash
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
### could not get anything to work except the conda command which created a medaka environment:
conda create -n medaka -c conda-forge -c nanoporetech -c bioconda medaka

##Busco
### Busco gave me errors every time I tried to download it, so I created a new environment for it and it worked
conda create --prefix=/work/gatins/hci_genome/busco
source activate /work/gatins/hci_genome/busco
conda install -c conda-forge -c bioconda busco=5.8.2

##Kraken
conda install bioconda::kraken2
### this actually installed the beta version... looking for the newest version
module load cmake/3.27.9
module load gcc/9.2.0
git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
make
./install_kraken2.sh /work/gatins/hci_genome/kraken2
cp /work/gatins/hci_genome/kraken2/kraken2{,-build,-inspect} /home/hughes.annab

##Blobtools
pip install "blobtoolkit[full]"
curl -L https://github.com/blobtoolkit/blobtoolkit/releases/download/4.1.5/blobtoolkit-api-linux > blobtoolkit-api
curl -L https://github.com/blobtoolkit/blobtoolkit/releases/download/4.1.5/blobtoolkit-viewer-linux > blobtoolkit-viewer
BTK_API_PORT=8880 BTK_PORT=8881 BTK_FILE_PATH=/work/gatins/hci_genome/processing/blobtools ./blobtoolkit-api
BTK_API_PORT=8880 BTK_PORT=8881 ./blobtoolkit-viewer

##RepeatModeler
module load singularity
curl -sSLO https://github.com/Dfam-consortium/TETools/raw/master/dfam-tetools.sh
chmod +x dfam-tetools.sh
./dfam-tetools.sh
