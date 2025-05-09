# Use Northeastern's large file transfer node
ssh hughes.annab@xfer.discovery.neu.edu

# For individual files, use wget
wget https://s3.amazonaws.com/gtl-public-data/miten/remy/03_24_25_R10_HCI_CUR_092401_GENOME.5mC_5hmC.sup.dorado.0.9.1.bam /work/gatins/hci_genome/bams
wget https://s3.amazonaws.com/gtl-public-data/miten/remy/03_24_25_R10_HCI_CUR_092401_GENOME_SS.5mC_5hmC.sup.dorado.0.9.1.bam /work/gatins/hci_genome/bams

# For recursive download, i.e. a full directory, use AWS CLI from command line
module load anaconda3/2022.05
source activate /work/gatins/hci_genome/env
pip install awscli

##flow cell 1
aws s3 cp s3://gtl-public-data/miten/remy/03_24_25_R10_HCI_CUR_092401_GENOME/03_24_25_R10_HCI_CUR_092401_GENOME/ /work/gatins/hci_genome/raw/ --recursive --no-sign-request
##flow cell 2 (size-selected)
aws s3 cp s3://gtl-public-data/miten/remy/03_24_25_R10_HCI_CUR_092401_GENOME/03_24_25_R10_HCI_CUR_092401_GENOME_SS/ /work/gatins/hci_genome/raw/ --recursive --no-sign-request

