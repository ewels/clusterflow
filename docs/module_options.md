Many modules can have their default behaviour modified through the use of
Cluster Flow `--params`. These are described below.

See the documentation about [Module Paramters](https://ewels.github.io/clusterflow/docs/#module-parameters)
for more information about how to specify these options.

## BEDTools intersectNeg
### `blacklistFile`
Use to define a blacklist file (overrides any set as a genome reference).
```bash
cf --params blacklistFile="/path/to/file"
```

## Bismark align
### `pbat`
Use the Bismark `--pbat` flag.
```bash
cf --params pbat
```

### `unmapped`
Save the unmapped reads to a file (Bismark `--unmapped` flag).
```bash
cf --params unmapped
```

### `bt1`
Align with Bowtie1 instead of Bowtie2 (default).
```bash
cf --params bt1
```

### `single_cell`
Use the `--non_directional` Bismark flag.
```bash
cf --params single_cell
```

### `subsample`
Only align the first 1000000 reads.
```bash
cf --params subsample
```

## Bowtie 1
### `mirna`
Use alignment paramters suitable for miRNA alignment against miRBase references,
instead of the standard Bowtie1 command.
Uses `-n 0 -l 15 -e 99999 -k 200` bowtie flags, instead of default `-m 1 --strata`.
```bash
cf --params mirna
```

## CF merge files
### `regex`
Override any merge regex set in the Cluster Flow configuration and use this instead.
```bash
cf --params regex="/REMOVE_([KEEP]+).fastq.gz/"
```

## deeptools bamCoverage
### `fragmentLength`
Set the fragment length to use for bamCoverage, instead of taking from the
phantompeaktools cross correlation analysis or using the default (200).
```bash
cf --params fragmentLength=120
```

## deeptools bamFingerprint
### `fragmentLength`
Set the fragment length to use for bamCoverage, instead of taking from the
phantompeaktools cross correlation analysis or using the default (200).
```bash
cf --params fragmentLength=120
```

## FastQ Screen
### `fastq_screen_config`
Use a specific FastQ Screen config file (with `--conf` FastQ Screen flag).
```bash
cf --params fastq_screen_config="/path/to/config"
```

## FastQC
### `nogroup`
Use the `--nogroup` option with FastQC to prevent automatic grouping of
base pair positions in plots. You can end up with some very large plots if
you have long reads!
```bash
cf --params nogroup
```

## featureCounts
### `stranded`
Set the `-s 1` flag for featureCounts.
```bash
cf --params stranded
```

### `stranded_rev`
Set the `-s 2` flag for featureCounts.
```bash
cf --params stranded_rev
```

### `id_tag`
Specify the tag to use for counting in the GTF file. If not specified,
module tries to guess by looking for a field called `gene_id` or `ID`.
```bash
cf --params id_tag="Gene"
```

## HiCUP
### `longest`
The longest fragment to accept (HiCUP parameter `--longest`). Default: `800`
```bash
cf --params longest=900
```
### `shortest`
The shortest fragment to accept (HiCUP parameter `--shortest`). Default: `100`
```bash
cf --params shortest=50
```
### `re1`
The restriction enzyme recognition pattern to use. Default: `"A^AGCTT,HindIII"`
```bash
cf --params re1="A^GATCT,BglII"
```

## HTSeq Counts
### `stranded`
Set the `-s yes` flag for HTSeq Counts. Default is to set `-s no`
```bash
cf --params stranded
```

### `stranded_rev`
Set the `-s reverse` flag for HTSeq Counts. Default is to set `-s no`.
```bash
cf --params stranded_rev
```

### `id_tag`
Specify the tag to use for counting in the GTF file. If not specified,
module tries to guess by looking for a field called `gene_id` or `ID`.
```bash
cf --params id_tag="Gene"
```

## Kallisto
### `estFragmentLength`
Specify the estimated fragment length (Kallisto `--fragment-length` option). Default: 200.
```bash
cf --params estFragmentLength=300
```
### `est_sd`
Specify the fragment length standard deviation (Kallisto `--sd` option). Default: 20.
```bash
cf --params est_sd=30
```

## MultiQC
### `template`
Specify the MultiQC template to use. Default: `default`
```bash
cf --params template=geo
```


## RSeQC (all modules)
### `keep_intermediate`
Do not delete the R files used to generate the PDF figures. Useful when running
downstream tools such as MultiQC, that use these intermediate files.
```bash
cf --params keep_intermediate
```

## Samtools sort + index
### `byname`
Sort by name instead of position (`-n` flag).
```bash
cf --params byname
```
### `forcesort`
Don't skip the sorting step, even if the file already seems to be sorted.
```bash
cf --params forcesort
```


## STAR
### `LoadAndRemove`
Load and remove genome index (`--genomeLoad LoadAndRemove`). Default: `NoSharedMemory`.
```bash
cf --params LoadAndRemove
```
### `LoadAndKeep`
Load and keep genome index (`--genomeLoad LoadAndKeep`). Default: `NoSharedMemory`.
```bash
cf --params LoadAndKeep
```
### `outSAMattributes`
Specify SAM attributes (`--outSAMattributes [attr]`). Default: `Standard`.
```bash
cf --params outSAMattributes="attr"
```


## TrimGalore!
### `min_readlength`
Minimum read length for trimming to run. If the first file in each run group
has reads less than this length, trimming will be skipped. Default: 50
```bash
cf --params min_readlength=30
```

### `force_trim`
Force TrimGalore! to run, even if reads are below minimum read length.
```bash
cf --params force_trim
```

### `q_cutoff`
Specify quality for trimming low-quality ends from reads in addition to adapter removal.
Default Phred score: 20.
```bash
cf --params q_cutoff=10
```

### `stringency`
Number of bases of overlap with adapter sequence required to trim a sequence.
Default: 1
```bash
cf --params stringency=3
```

### `adapter`
Specify an adapter sequence to trim. Default: Auto-detect (_Illumina universal_,
_Nextera transposase_ or _Illumina small RNA adapter_).
```bash
cf --params adapter=ATACAGCTAGCAGTAC
```

### `RRBS`
Specifies that the input file was an MspI digested RRBS sample.
```bash
cf --params RRBS
```

### `nofastqc`
Do not run FastQC after trimming is complete.
```
cf --params nofastqc
```

### Specific trimming
To remove a custom number of bases from reads after adapter removal,
the following parameters can be set:

* `cf --params clip_r1=<int>`
  * Remove bp from the 5' end of read 1 (or single-end reads).
* `cf --params clip_r2=<int>`
  * Remove bp from the 5' end of read 2 (paired-end only).
* `cf --params three_prime_clip_r1=<int>`
  * Remove bp from the 3' end of read 1 _AFTER_ adapter/quality trimming has been performed.
* `cf --params three_prime_clip_r2=<int>`
  * Remove bp from the 3' end of read 2 _AFTER_ adapter/quality trimming has been performed.

The following params are presets which are easier to remember and use:

* `cf --params trim=<int>`
  * Trim from 5' of R1 and R2. Equivalent to `clip_r1=<int> clip_r2=<int>`.
* `cf --params pbat`
  * `clip_r1 6`
  * `clip_r2 6`
* `cf --params ATAC`
  * `clip_r1 4`
  * `clip_r2 4`
* `cf --params single_cell`
  * `clip_r1 9`
  * `clip_r2 9`
* `cf --params epignome`
  * `clip_r1 7`
  * `clip_r2 7`
  * `three_prime_clip_r1 7`
  * `three_prime_clip_r2 7`
* `cf --params accel`
  * `clip_r1 10`
  * `clip_r2 15`
  * `three_prime_clip_r1 10`
  * `three_prime_clip_r2 10`
* `cf --params cegx`
  * `clip_r1 6`
  * `clip_r2 6`
  * `three_prime_clip_r1 2`
  * `three_prime_clip_r2 2`






