## Installation
### Environment Module
If you have Cluster Flow installed on your cluster as an environment module,
you can load it using:
```
module load clusterflow
```

### Manual Installation
If you want to install Cluster Flow yourself manually, follow these steps:

1. [Download Cluster Flow](https://github.com/ewels/clusterflow/releases)
```
wget https://github.com/ewels/clusterflow/archive/v0.4.tar.gz
```
2. Extract the files
```
tar -C clusterflow -zxvf v0.4.tar.gz
```
3. Create & configure the site-wide configuration file
```
cd clusterflow
cp clusterflow.config.example clusterflow.config
nano clusterflow.config
```
You must specify your environment (running locally, or the type of cluster).
See the comments in the config file for explanation of the other options.

## Configuration
Once Cluster Flow has been set up site-wide, you need to configure it for your
personal use:

	cf --setup

This will launch a wizard to write a config file for you, with details such
as e-mail address and notification settings.

## Adding Reference Genomes
Most analysis pipelines need a reference genome. This can exist in a central
location or in your personal setup (or both).

> If you're using the Swedish UPPMAX cluster, please see
> [these instructions](https://github.com/ewels/clusterflow-uppmax).

You can add your reference genome paths with the following wizard:

	cf --add_genome

## Do a test run!
That should be it! If you added some bash aliases such as `qs`, you should
log out and log in again. Then you can try running a test run:

    cf --genome GRCh37 sra_bowtie ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX031/SRX031398/SRR1068378/SRR1068378.sra

This will download [SRR1068378](http://www.ncbi.nlm.nih.gov/sra/?term=SRR1068378)
_(Human H3K4me3 ChIP-Seq data)_, convert to FastQ, run FastQC, Trim Galore! and align with bowtie.

See the [usage instructions](usage) for more details on how to get started
using Cluster Flow.
