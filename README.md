# Context 

Mappability can be challenging in repetitive regions when analyzing NGS-derived data. Local **duplication of genomic** regions hampers the alignment of reads in those parts of the genome, since it **prevents unambiguously** **alingning reads**. Therefore, those reads are usually assigned a random location (among the possible ones) and a low _mapping quality_. Since reads are assigned a random location, **mutations** **are usually either** **missed** (mutated reads are dispersed among all possible alignment regions) or, even worse, **called multiple times** leading to the rise of false positives (enough reads to call a variant in multiple alignment regions). For this reason, **low mapping quality regions are often ignored** in variant calling pipelines.

__Armadillo__ is a somatic variant calling pipeline that aims at repetitive regions for tumor-control paired samples, in order to recover mutations that are systematically missed. In order to avoid both false positives and false negatives, it first __finds all candidate regions__ for a repetitive region of interest (RROI). Then, it extracts all reads coming from all regions and aligns them at the RROI sequence only instead of the whole genome, forcing all reads to be aligned together. This way, dilution of variant-supporting reads is prevented, as well as the false positives caused by mistakenly call a mutation at multiple copies. 
Finally, the variant calling step is carried out. After usual __heuristic filters__, we use __context-specific filters__ to prevent false positives related the stacking step and, lastly, run a __CNN model to classify the mutations__.

## Armadillo data-prep

Armadillo data-prep is a small tool to make easy building the reference file needed by [armadillo](https://platform.batchx.io/uniovi/tools/labxa%2Farmadillo), our main tool. This tool was thought to be pipelined with [getFasta](https://platform.batchx.io/batchx/tools/bioinformatics%252Fbedtools%252Fgetfasta) and [BLAT](https://platform.batchx.io/kent-informatics/tools/blat%2Falign/). This way, with a simple BED file and a FASTA file, the reference file needed by Armadillo is built easily.
 
More information on the pipeline can be retrieved in our [GitHub repository](https://github.com/xa-lab/armadillo) or in [our paper](https://www.nature.com/articles/s41525-022-00292-2).

# Inputs

## Required inputs

1. `blast8`: Input file. It's the BLAT output in blast8 format for all regions needed. 
2. `reference`: Fasta file of species genome.

## Optional inputs

1. `referenceIdx`: Reference's index (fai file). Recommended for shorter runs
2. `identity`: Minimum identity for a hit to be considered a copy of the regions of interest (ROIs)
3. `lendiff`: Maximum %length difference allowed between each hit and the input sequence.
4. `mlen`: Minimum length of a ROI.
7. `outputName`: Basename of the output file


# Outputs

1. `armadillo_data`: Armadillo reference file in gz format.

# Example 

1. Copy bams and ROIs to batchX
```bash
bx cp rois.bed hg19.fa rois bx://test/
```

> A list of precomputed rois of interest is available in our github repo for [hg19](https://github.com/pbousquets/armadillo-batchx/tree/master/lib/armadillo_data_hg19/rois) or [here](https://github.com/pbousquets/armadillo-batchx/tree/master/lib/armadillo_data_hg38/rois) for hg38. 

3. Submit job

```bash
bx submit -v=1 -m=4000 uniovi@labxa/armadillo-dataprep:1.0.0 \
'{
  "blast8": "test",
  "reference": "/test/rois.bed",
  "tumorGenome": "/test/hg19.fa"}'
```

# Tools version

* [Samtools](https://www.htslib.org/doc/1.9/samtools.html) (built with v.1.9)
* [Python3](https://www.python.org) (built with v.3.7)

# Authors

* **Pablo Bousquets Muñoz** - bousquetspablo@uniovi.es
* **Ander Díaz Navarro**
* **Xose Antón Suárez Puente**
