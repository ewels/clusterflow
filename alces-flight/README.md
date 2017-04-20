# Cluster Flow on Alces Flight

[Alces Flight](http://alces-flight.com/) provides a super easy-to-use
HPC environment running on Amazon Web Services that will be familiar to many
coming from traditional academic cluster environments.

Alces Flight comes with loads of pre-installed software, many of which are
bioinformatics related. As such, it's perfectly suited for running Cluster Flow.

We are [currently working](https://community.alces-flight.com/t/module-request-cluster-flow/303/13)
to get Cluster Flow set up as an available package on Alces Flight to make it
a two-command setup procedure. But in the mean time, here is a run through
of how to use Cluster Flow on Alces Flight.

## Step 1: Create an Alces Flight Cluster

First, we need to launch an [Alces Flight](http://alces-flight.com) system.

1. Go to the Community Edition Alces Flight [AWS Market Page](https://aws.amazon.com/marketplace/pp/B01GC9E3OG?ref=_ptnr_www)
2. Select a nearby region (eg. `EU (Ireland)` for me) and choose _Personal HPC compute cluster_ from the _Delivery Method_ dropdown. Click _Continue_.
3. Log in to AWS
4. Leave the defaults and hit _Launch with CloudFormation Console_
5. Go through the wizard selecting your options.
    * We won't do any work on the login node, so it can be very small (and cheap). Select `other` and `t2.micro-1C-1GB`
        * I think that if you're on the [AWS Free Tier](https://aws.amazon.com/free/) you can then leave this head node running for a year for free? I'm not certain though.
    * The compute instance node type needs to be bigger. How much bigger depends on what analysis you're going to run. `memory-8C-60GB.medium-r3.2xlarge` should be fine for most things.
    * We highly recommend using [Spot Pricing](https://aws.amazon.com/ec2/spot/pricing/) - you run a small risk of the analysis failing because of getting outbid, but it should cost far far less than the AWS list pricing (typically around 80% savings).
        * If you set a high bid then you're almost certain not to get knocked off, but will pay very little most of the time.
        * You can see the above link for the current spot pricing. At the time of writing, `r3.2xlarge` spot pricing is `$0.1482 per Hour`. This is in contrast to `$0.741 per Hour` for the [on-demand pricing](https://aws.amazon.com/ec2/pricing/on-demand/). So entering `0.741` or more into the `Spot price` field should be fine for ensuring that you won't be knocked off. I'm often lazy and just put `1` (you only ever pay the current bid price).
    * We recommend setting the _Autoscaling policy_ to `enabled`. This means that new compute nodes will be created and destroyed according to how many jobs are currently sitting in the queue. Cool stuff!
        * If you set the _Initial compute nodes_ note that it'll take a few minutes after you launch a Cluster Flow pipeline for anything to actually run, as a new node will need to be created.
    * Be careful about the _Login node system volume size_. If you have a lot of data then you'll probably want to mount and run on a separate file system (eg. EFS is more expensive, but doesn't have a size limit).
6. Hit Create!
7. In the CloudFormation _Stacks_ page that you're taken to, click the refresh button in the top right until your new stack appears.
8. Select your stack to get the details tabs to show below.
9. Once the status turns green (`CREATE_COMPLETE`) and your Cluster is ready, SSH in using a terminal.
    * Click the _Outputs_ tab and copy the `AccessIP` field. Note also the `Username`.
    * Use these to SSH in. eg: `ssh -i "something@keypair.pem" Username@AccessIP`
    * Note that if you didn't set up your keypairs properly in the setup wizard (or access IP addresses) then you will have a bad time.

## Step 2: Install Cluster Flow

> The team behind Alces Flight are currently working towards releasing Cluster Flow
> as a pre-configured app. Hopefully, this will essentially mean that all of the following
> section can be achieved with a single `alces gridware install main/apps/clusterflow` command
> in the near future.

#### Load Perl and install some dependencies
First, we need to load a Perl installation and add a couple of extra dependencies that
are used by Cluster Flow:
```bash
alces gridware install main/apps/perl/5.20.2
module load apps/perl

# Install the required XML::Simple Perl package
module load services/pdsh
pdsh -g cluster 'sudo yum -y -e0 install expat-devel'
cpan XML::Simple

# This is used by Cluster Flow to load environment modules
sudo ln -s /opt/clusterware/opt/modules/bin/modulecmd /usr/local/bin
```

#### Install Cluster Flow
Next we need to get hold of a copy of Cluster Flow itself.

```bash
wget https://github.com/ewels/clusterflow/archive/master.zip
unzip master.zip
rm master.zip
cd clusterflow-master/
export PATH=$PATH:$(pwd)
mv alces-flight/clusterflow.config clusterflow.config
cd ~
```

#### Configure Cluster Flow
To enter your e-mail address for notification e-mails and to
set other config options, run the setup wizard:

```bash
cf --setup
```

#### Install required Alces Flight Gridware Software
Now we will install all of the programs used by the Cluster Flow modules.
Once installed, they're available to load as environment modules (done by
Cluster Flow automatically when a pipeline is run).

```bash
# Install and load required software
alces gridware install main/apps/bedtools/2.25.0
alces gridware install main/apps/star/2.5.2a
alces gridware install main/apps/bismark/0.16.3
alces gridware install main/apps/bowtie/1.1.0
alces gridware install main/apps/bowtie2/2.3.0
alces gridware install main/apps/bwa/0.7.15
alces gridware install main/apps/cutadapt/1.8
alces gridware install main/apps/fastqscreen/0.6.1
alces gridware install main/apps/fastqc/0.11.3
alces gridware install main/apps/htseq/0.6.1p1
alces gridware install main/apps/kallisto/0.43.0
alces gridware install main/apps/phantompeakqualtools/1.1
alces gridware install main/apps/picard/2.6.0
alces gridware install main/apps/samtools/1.4
alces gridware install main/apps/sra/2.3.5-2
alces gridware install main/apps/subread/1.5.0-p3
alces gridware install main/apps/tophat/2.1.0
alces gridware install main/apps/trimgalore/0.4.2

# Volatile packages
sed -i 's?^# - /opt/clusterware/var/lib/gridware/repos/volatile? - /opt/clusterware/var/lib/gridware/repos/volatile?g' /opt/gridware/etc/gridware.yml
alces gridware update
alces gridware install volatile/apps/multiqc/0.9

# These currently fail: "ERROR: Unable to download source."
alces gridware install volatile/apps/preseq/1.0.0
alces gridware install volatile/apps/rseqc/2.4

# Missing packages:
# - deepTools
# - hicup
# - hisat2
# - salmon
```

## Step 3: Set up your reference genomes
### Provided references
We have created a publicly accessibly S3 bucket containing the reference genomes and indices required
for a number of reference genomes, which you can use when running Cluster Flow. This is based on the
[illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html) resource,
uncompressed and with a few extra indices added (STAR, Bismark). Pulling your required references from this
is the fastest and easiest way to run alignments on Alces Flight using Cluster Flow. Note that the S3 bucket
is set to use _Requester Pays_ policy. We pay to keep the data hosted, but you will pay for any fees associated
with accessing the resource.

> **NB:** This is working out to be more expensive that we had hoped. It may well have been taken down by
> the time that you read this. We hope that Amazon will provide a similar public dataset resource in the near
> future and will update this documentation if / when that happens.

First, pull your required data from the S3 bucket. What you will need depends on what species you want to
align to and what tools you will be using. You can find the available species and references in the bundled
[`s3-igenomes.txt`](s3-igenomes.txt) file. Each reference listed contains the following files
_(may vary a little - Human Ensembly GRCh37 shown below)_:

```
Annotation/README.txt
Annotation/Genes/
Annotation/SmallRNA/
Annotation/Variation/
Sequence/AbundantSequences/
Sequence/BismarkIndex/
Sequence/Bowtie2Index/
Sequence/BowtieIndex/
Sequence/BWAIndex/
Sequence/Chromosomes/
Sequence/STARIndex/
Sequence/WholeGenomeFasta/
```

You now need to sync the required files to your local disk space. For example, if you intend to run
the `fastq_star` pipeline to align Human RNA-seq reads, you will need the following commands:
```bash
# Make a directory to hold your references
mkdir references

# Configure s3 access
s3cmd --configure
# Use User access key from AWS IAM Users console

# Sync the s3 bucket data to your local disk
s3cmd --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/ references/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/
s3cmd --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/STARIndex/ references/Homo_sapiens/Ensembl/GRCh37/Sequence/STARIndex/

# Tell Cluster Flow about your reference indices
cf --add_genome
# ..use the wizard to search this folder and find references automatically
```

> **NB:** We are currently talking with Amazon to make the above references into a public resource.
> If this goes through as planned, it may be possible to pre-configure the Alces Fight installation
> of Cluster Flow with ready-to-use references. This will mean that you can skip the above section - neat!

### Adding a reference manually
If the above resource is no longer available or doesn't contain the reference you need, you will need to fetch your
reference genome manually. For testing purposes, you can use the below dataset:

```bash
# Download some test data
wget https://export.uppmax.uu.se/b2013064/test-data/ngi-rna_test_set.tar.bz2
tar xvjf ngi-rna_test_set.tar.bz2
rm ngi-rna_test_set.tar.bz2
cd ngi-rna_test_set/
pwd # copy path

# Set up a reference genome
cf --add_genome
# (all users, Yeast, test, test, paste, y,y,y,y,y,.... enter)
```

### Testing added genomes
If everything has worked properly, running `cf --genomes` should show you something like this:

```
================================
Cluster Flow - available genomes
================================

----------------------------------------------------------------------------------------------------------------------
 /home/alces/clusterflow-master/genomes.config
 Name           Type      Species     Assembly       Path
----------------------------------------------------------------------------------------------------------------------
 GRCh37         fasta     Human       GRCh37         references/Homo_sapiens/Ensembl/GRCh37/Sequence/STARIndex
 GRCh37         gtf       Human       GRCh37         references/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf
 GRCh37         star      Human       GRCh37         references/Homo_sapiens/Ensembl/GRCh37/Sequence/STARIndex
----------------------------------------------------------------------------------------------------------------------
```

## Step 4: Get your raw data and run Cluster Flow
Everything should now be installed, configured and ready to run! You just need some data to process.
You can sync data from your local computer to the cluster as follows:

```bash
mkdir data
```

Then on your local machine (example only):
```
scp -i "something@keypair.pem" *.fastq.gz Username@AccessIP:/home/alces/data/
```

If you're just testing Cluster Flow out, you can use the `ngi-rna_test_set/` described above.

Now back on the cluster, you can kick off a Cluster Flow analysis run! For example:
```bash
cd data
cf --genome GRCh37 fastq_star data/*.fastq.gz
```

## Step 5: Download processed data
Hopefully if everything goes well, you will get an e-mail appear in your inbox
as normal saying that the analysis is complete. You can use the reverse `scp` above
or something else (sync to a s3 bucket?) to pull down your results files, and you're done!

For example (on your local machine):
```
mkdir my_results/
scp -i "something@keypair.pem" -r Username@AccessIP:/home/alces/data/* my_results/
```

> **NB:** This will re-download your raw data. Remember to move this out of your download data
> directory first (or to use `rsync` / something more clever than `scp`).

## Step 6: Shut down the cluster
Done except for one last thing that is: **you must shut down your cluster!**. If you don't,
then it will keep running for ages and cost you loads of money.

> TODO: Investigate whether running the head node alone on a micro instance
> costs any money on the free tier. Would be great if that could be left alive,
> just kicking off compute nodes when required and killing them when things are
> done..

1. Go back to your AWS CloudFormation console page, with the _Stacks_ page showing your Alces Flight cluster stack.
2. Select the row with your stack and click the _Actions_ dropdown. Click _Delete Stack_.
3. Hope that you downloaded everything you needed and wait for everything to be killed!



