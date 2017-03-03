# Cluster Flow Installation

For full installation instructions, please see the
[`docs/installation.md`](docs/installation.md) file or read
online at [http://clusterflow.io/docs/#installation](http://clusterflow.io/docs/#installation)

## Quick Start
Cluster Flow is written in Perl. To install,
[download from GitHub](https://github.com/ewels/clusterflow/releases))
and move to the desired installation location. Extract the files:

```
tar -C clusterflow -zxvf v0.4.tar.gz
```

Before running Cluster Flow, you need to configure it to work with your
system. Copy the example `clusterflow.config.example` config file to
a new file called `clusterflow.config` (in the main CF script directory).
Edit it as appropriate, ensuring that you specify `@cluster_environment`.

Now add the main `cf` executable to the `PATH`. You can do this by symlinking
the executable to your `bin` directory or adding a bash alias to your
`~/.bashrc` file:

```
alias cf='/path/to/cf'
```

Now that the central Cluster Flow installation is complete, you can set
up your personal config using the bundled setup tool:

```
cf --setup
```

Finally, use the reference genome tool to create a config file describing paths
to your reference genome paths:

```
cf --add_genome
```

