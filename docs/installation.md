## Requirements
Cluster Flow is designed to work with a computing cluster. It currently supports
the Sun GRIDEngine, LSF and SLURM job managers (not PBS, Torque or others).

If you don't have a cluster with a supported manager, you can run Cluster Flow on
any command-line machine in `local` mode. This writes a bash script and runs it as
a job in the background.

To run analyses, you will also need the required tools to be installed. Cluster Flow
is designed to work with the environment module system and load tools as required, but
if software is available on the `PATH` it can work without this.

Cluster Flow itself is written in Perl. It has minimal dependencies, all of which are
core Perl packages.

## Environment Module
If you are a user on a HPC cluster, you may already have Cluster Flow installed
on your cluster as an environment module. If so, you may be able to load it using:
```
module load clusterflow
```

## Manual Installation
Cluster Flow is a collection of stand-alone scripts, mostly written in Perl.

1. Download Cluster Flow (see the [releases page](https://github.com/ewels/clusterflow/releases))
```bash
wget https://github.com/ewels/clusterflow/archive/v0.5.tar.gz
```
2. Extract the files
```bash
tar xvzf v0.5.tar.gz
```
3. Create & configure the site-wide configuration file
```bash
cd clusterflow-0.5
cp clusterflow.config.example clusterflow.config
vi clusterflow.config
```

You must specify your environment in the config file (`@cluster_environment`:
`local`, `GRIDEngine`, `SLURM` or `LSF`), most other things are optional.

The `cf` executable must be in your system `PATH`, so that you can run it easily
from any directory. Ensure that you run the Configuration Wizard (described below)
so that this config is created in your `~/.bashrc` file.

If you prefer, you can symlink the `cf` executable to `~/bin`.

## Configuration Wizard
Once Cluster Flow has been set up site-wide, you need to configure it for your
personal use:
```bash
cf --setup
```

This will launch a wizard to write a config file for you, with details such
as e-mail address and notification settings.

## Adding Reference Genomes
Most analysis pipelines need a reference genome. This can exist in a central
location or in your personal setup (or both).

> If you're using the Swedish UPPMAX cluster, please see
> [these instructions](https://github.com/ewels/clusterflow-uppmax).

You can add your reference genome paths with the following wizard:
```bash
cf --add_genome
```

## Do a test run!
That should be it! Log out of your session and in again to activate any new
bash settings. Then try launching a test run:
```bash
cf --genome GRCh37 sra_bowtie ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX031/SRX031398/SRR1068378/SRR1068378.sra
```

This will download [SRR1068378](http://www.ncbi.nlm.nih.gov/sra/?term=SRR1068378)
_(Human H3K4me3 ChIP-Seq data)_, convert to FastQ, run FastQC, Trim Galore! and align with bowtie.

