```bash

# Install Perl & dependencies
alces gridware install main/apps/perl/5.20.2
module load apps/perl

# These are needed by XML::Simple
sudo -s
yum install expat
yum install expat-devel
exit

cpan XML::Simple


# Get Cluster Flow
wget https://github.com/ewels/clusterflow/archive/v0.5.tar.gz
tar xvzf v0.5.tar.gz
rm v0.5.tar.gz
cd clusterflow-0.5/
export PATH=$PATH:$(pwd)

# Configure
cat > clusterflow.config
# paste contents of clusterflow_aws.config

# Optional - depends on how customised the main config file is
cf --setup

# Download some test data
wget https://export.uppmax.uu.se/b2013064/test-data/ngi-rna_test_set.tar.bz2
tar xvjf ngi-rna_test_set.tar.bz2
rm ngi-rna_test_set.tar.bz2
cd ngi-rna_test_set/
pwd # copy path

# Set up a reference genome
cf --add_genome
# (all users, Yeast, test, test, paste, y,y,y,y,y,.... enter)
cf --genomes

# ================================
# Cluster Flow - available genomes
# ================================
# 
# -----------------------------------------------------------------------------------------------------------------
#  /home/phil/clusterflow/clusterflow-master/genomes.config
#  Name           Type      Species                  Assembly       Path
# -----------------------------------------------------------------------------------------------------------------
#  test           fasta     Yeast                    test           /home/phil/genomes/ngi-rna_test_set
#  test           gtf       Yeast                    test           /home/phil/genomes/ngi-rna_test_set/genes.gtf
#  test           star      Yeast                    test           /home/phil/genomes/ngi-rna_test_set/star
# -----------------------------------------------------------------------------------------------------------------

# Install and load required software
alces gridware install main/apps/fastqc/0.11.3
alces gridware install main/apps/trimgalore/0.4.2
alces gridware install main/apps/samtools/1.4
alces gridware install main/apps/star/2.5.2a
module load apps/fastqc/0.11.3
module load apps/trimgalore/0.4.2
module load apps/samtools/1.4
module load apps/star/2.5.2a

# Run test run
cf --genome test fastq_star *.fastq.gz

```