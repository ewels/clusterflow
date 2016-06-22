---
Using Cluster Flow:
  Installation: installation.md
  General Usage: usage.md
  Command Line Reference: cl_reference.md
  Troubleshooting: troubleshooting.md
Coding with Cluster Flow:
  Writing new pipelines and modules: writing_pipelines_modules.md
---

# Welcome!
## Cluster Flow v0.4 Documentation

Cluster Flow is simple package
designed to run bioinformatics pipelines. It is operated through a single command
`cf`, which can be used to launch, configure, monitor and cancel pipelines.

When you run Cluster Flow, you choose a pipeline. This has a list of modules
in order. Each module is a wrapper around a bioinformatics tool.
When the pipeline has finished, a notification e-mail is sent to you with
status messages from the log.

## Tutorial Videos

[Usage / Installation Tutorial](http://youtu.be/b2g_zQiz9ys) | [Advanced Tutorial](http://youtu.be/aBHOcsA2M6w)
-------------------------------------------------------------|-------------------------------------------------------
<iframe width="300" height="169" src="https://www.youtube.com/embed/b2g_zQiz9ys?rel=0&amp;showinfo=0" frameborder="0" allowfullscreen></iframe> | <iframe width="300" height="169" src="https://www.youtube.com/embed/aBHOcsA2M6w?rel=0&amp;showinfo=0" frameborder="0" allowfullscreen></iframe>

See the Installation page for setup instructions.

## Contributing to Cluster Flow
If you write a module or pipeline which could be of use to others, it would be
great to merge those changes back into the core Cluster Flow project.

For instructions on how best to do this, please see the
[contributing instructions](https://github.com/ewels/clusterflow/blob/master/CONTRIBUTING.md).

