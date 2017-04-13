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

This is the stuff that will hopefully be redundant if/when Cluster Flow comes bundled with
the Alces Flight Gridware.

#### Load Perl and install some dependencies
```bash
alces gridware install main/apps/perl/5.20.2
module load apps/perl

sudo -s
# This is needed by XML::Simple
yum install expat-devel
# This is used by Cluster Flow
ln -s /opt/clusterware/opt/modules/bin/modulecmd /usr/local/bin
exit

cpan XML::Simple
```

#### Install Cluster Flow

```bash
wget https://github.com/ewels/clusterflow/archive/v0.5.tar.gz
tar xvzf v0.5.tar.gz
rm v0.5.tar.gz
cd clusterflow-0.5/
export PATH=$PATH:$(pwd)
```

#### Configure Cluster Flow
Find the [Alces Flight Cluster Flow config file](https://raw.githubusercontent.com/ewels/clusterflow/alces-flight/clusterflow_aws.config)
and copy the contents.

```bash
cat > clusterflow.config
# paste contents of clusterflow_aws.config
# ctrl-D to complete
cd ../

# Run setup to add your e-mail address and preferences
cf --setup
```

#### Install required Alces Flight Gridware Software
```bash
# Install and load required software
alces gridware install main/apps/fastqc/0.11.3
alces gridware install main/apps/trimgalore/0.4.2
alces gridware install main/apps/samtools/1.4
alces gridware install main/apps/star/2.5.2a

# .. plus everything else. Exhaustive list coming soon.
```

## Step 3: Set up your reference genomes

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
```

## Step 4: Get your raw data and run Cluster Flow

```bash
# Run test run
cf --genome test fastq_star *.fastq.gz
```