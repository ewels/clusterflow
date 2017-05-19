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
    * If you would like to use the [AWS-iGenomes](https://github.com/ewels/AWS-iGenomes) reference genomes, we highly recommend selecting `EU (Ireland)`. The iGenomes are stored in that region, so transfers will be free if your cluster is also there.
3. Log in to AWS
4. Leave the defaults and hit _Launch with CloudFormation Console_
5. Go through the wizard selecting your options.
    * The first page can be left as it is, just click _Next_.
    * **Access and security**
        * Ensure that you have a keypair selected that you'll be use - otherwise you won't be able to log into your cluster! See the [AWS docs](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html) for more information.
        * Set _Access network address_ to an IP range that you can access. If in doubt, specify `0.0.0.0/0` so that anyone can access the login node (they'll still need SSH authentication using the keypair above).
    * **Configuration and customization**
        * Under _Additional features to enable_ enter `clusterflow`. This tells Alces Flight to configure a Cluster Flow installation for your cluster setup.
    * **Front end node configuration**
        * We won't do any work on the login node, so it can be very small (and cheap). Select `other` and `t2.micro-1C-1GB`
        * I _think_ that if you're on the [AWS Free Tier](https://aws.amazon.com/free/) you can then leave this head node running for a year for free? Do this at your own risk though!
    * **Compute estate**
        * The compute instance node type needs to be bigger. How much bigger depends on what analysis you're going to run. `memory-8C-60GB.medium-r3.2xlarge` should be fine for most things.
        * We highly recommend using [Spot Pricing](https://aws.amazon.com/ec2/spot/pricing/) - you run a small risk of the analysis failing because of getting outbid, but it should cost far far less than the AWS list pricing (typically around 80% savings). Enter a value into _Spot price_ to enable this.
            * If you set a high bid then you're almost certain not to get knocked off, but will pay very little most of the time.
            * You can see the above link for the current spot pricing. At the time of writing, `r3.2xlarge` spot pricing is `$0.1482 per Hour`. This is in contrast to `$0.741 per Hour` for the [on-demand pricing](https://aws.amazon.com/ec2/pricing/on-demand/). So entering `0.741` or more into the `Spot price` field should be fine for ensuring that you won't be knocked off. I'm often lazy and just put `1` (you only ever pay the current bid price).
        * We recommend setting the _Autoscaling policy_ to `enabled`. This means that new compute nodes will be created and destroyed according to how many jobs are currently sitting in the queue. Cool stuff!
            * If you set the _Initial compute nodes_ to `0` then all of the expensive compute nodes should be automatically shut down when the pipeline completes (leaving only your cheap head node running).
            * Note that there is a slight lag for the autoscaling policy, so it'll take a few minutes after you launch a Cluster Flow pipeline for anything to actually run, whilst you wait for a new node to be created.
    * **Disks and storage**
        * Be careful about the _Login node system volume size_. If you have a lot of data then you'll probably want to mount and run on a separate file system (eg. EFS is more expensive, but doesn't have a size limit).
    * Click through the rest of the wizard - defaults should be fine for the rest of the pages.
    * Hit Create!
6. In the CloudFormation _Stacks_ page that you're taken to, click the refresh button in the top right until your new stack appears.
7. Select your stack to get the details tabs to show below.
    * If you can't see the tabs below, click the little square icons in the bottom right to bring it up.
8. Once the status turns green (`CREATE_COMPLETE`) and your Cluster is ready, SSH in using a terminal.
    * Click the _Outputs_ tab and copy the `AccessIP` field. Note also the `Username`.
    * Use these to SSH in. eg: `ssh -i "something@keypair.pem" Username@AccessIP`
    * Note that if you didn't set up your keypairs properly in the setup wizard (or access IP addresses) then you won't be able to get in and need to start again :(

## Step 2: Configure Cluster Flow

1. _Optional:_ If you created your cluster with 0 compute nodes, we need to create one so that Gridware can properly configure Cluster Flow.
    * Running the `qhost` command should show just `global` and no compute nodes.
    * Launch any job with the scheduler. For example:
        ```
        echo "sleep 1800" | qsub -N sleep
        ```
    * You'll probably get a message like _warning: no suitable queues_. That's fine - it's just telling you that no compute nodes are available.
    * Wait (~10 minutes) until you have some compute nodes running - check using the command `qhost`. For example:
        ```
        HOSTNAME                ARCH         NCPU  LOAD  MEMTOT  MEMUSE  SWAPTO  SWAPUS
        -------------------------------------------------------------------------------
        global                  -               -     -       -       -       -       -
        flight-108              linux-x64       8  1.23   59.6G  619.1M     0.0     0.0
        ```
2. Load the Cluster Flow Alces Flight Gridware module.
    * Load using the `module load` command:
        ```
        module load apps/clusterflow
        ```
    * You should get a notice saying something like this:
        ```
        Cluster Flow has been configured:

          Job scheduler: SGE
                  Cores:   8
                 Memory:   59G
        ```
    * If _Cores_ and _Memory_ don't show a value, you probably didn't have any compute nodes running when you loaded  the module. Things will probably break. You probably need to start again. Sorry.
3. Configure Cluster Flow with your personal setup (eg. e-mail address):
    * Run the Cluster Flow command line setup wizard:
        ```
        cf --setup
        ```
    * There's no need to add anything to `/home/alces/.bash_profile`, so you can say _no_ to all of those prompts.

## Step 3: Set up your reference genomes
### AWS iGenomes on S3
To help running analyses on AWS, we have created an S3 bucket containing the reference genomes and indices required for a number of reference genomes, which you can use when running Cluster Flow. This is based on the [illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html) resource, uncompressed and with a few extra indices added (STAR, Bismark).

Pulling your required references from this is the fastest and easiest way to run alignments on Alces Flight using Cluster Flow. Note that the S3 bucket is set to use _Requester Pays_ policy. We pay to keep the data hosted, but you will pay for any fees associated with accessing the resource. The S3 bucket is in _EU (Ireland)_ - if your cluster is in the same region then any transfers should be free. Our hosting costs for this are kindly paid by Amazon through a research grant.

1. First, you need to configure the `s3cmd` on your cluster with your authentication.
    * Run the following command and follow the prompts:
        ```
        s3cmd --configure
        ```
2. Next, run the AWS-iGenomes command line wizard to sync the required genomes:
    ```bash
    curl -fsSL https://ewels.github.io/AWS-iGenomes/aws-igenomes.sh | bash
    ```
3. Once downloaded, configure Cluster Flow so that it knows about these genomes:
    ```bash
    cf --add_genome
    ```

You can find more information and documentation about this resource at https://github.com/ewels/AWS-iGenomes

### Adding a reference manually
If AWS-iGenomes doesn't have the reference you need, you will need to fetch your
reference genome manually. For example, you can upload it from your local computer as follows:

```bash
scp -r -i "something@keypair.pem" /my/reference/ Username@AccessIP:/home/alces/references/
```

Then on the AF cluster, tell Cluster Flow about your reference:
```bash
cf --add_genome
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

Then on your local machine: _(example only)_
```
scp -i "something@keypair.pem" *.fastq.gz Username@AccessIP:/home/alces/data/
```

Now back on the cluster, you can kick off a Cluster Flow analysis run! For example:
```bash
cf --genome GRCh37 fastq_star data/*.fastq.gz
```

## Step 5: Download processed data
Hopefully if everything goes well, you will get an e-mail appear in your inbox
as normal saying that the analysis is complete. You can download the results using
`scp` to pull down your results files, and you're done!

For example (on your local machine):
```
mkdir my_results/
scp -i "something@keypair.pem" -r Username@AccessIP:/home/alces/data/* my_results/
```

> **NB:** This will re-download your raw data. Remember to move this out of your download data
> directory first (or to use `rsync` / something more clever than `scp`).

## Step 6: Shut down the cluster
Remember - when you're finished **you must shut down your cluster!**. If you don't,
then it will keep running for ages and cost you loads of money.

1. Go back to your AWS CloudFormation console page, with the _Stacks_ page showing your Alces Flight cluster stack.
2. Select the row with your stack and click the _Actions_ dropdown. Click _Delete Stack_.
3. Hope that you downloaded everything you needed and wait for everything to be killed.

> If you're running the tiny head node, have no compute nodes running and are on
> the AWS free tier, then it may not cost you any money to leave the head node
> running, which could be convenient. But try this at your own risk!

