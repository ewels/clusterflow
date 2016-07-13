## Pipeline syntax
All pipelines conform to a standard syntax. The name of the pipeline is
given by the filename, which should end in `.config`. The top of the file
should contain a title and description surrounded by `/*` and `*/`

Variables can be set using the same `@key value` syntax as in
`clusterflow.config` files.

Modules are described using `#` prefixes. Tab indentation denotes dependencies
between modules. Syntax is `#module_name parameters`, where there can be any
number of space separated parameters which will be passed on to the module
at run time.

### Example pipeline
Here is an example pipeline, which requires a genome path and uses
three modules:

```
/*
Example Pipeline
================
This pipeline is an example of running three modules which depend on
each other. Module 2 is run with a parameter that modifies its behaviour.
This block of text is used when cf --help example_pipeline is run
*/
#module1
       #module2
       #module2 parameter
             #module3
```

Remember to run `dos2unix` on your pipeline before you run it, if you're
working on a windows machine.

## Run files
Cluster Flow works by creating `.run` files for each batch of input files.
These are a copy of the pipeline file, with filenames appended for each
step of the pipeline. These files are used by subsequent steps in the
pipeline to know which input files to use.

Inspecting run files is a quick way to see exactly what analysis was
done in a directory.
