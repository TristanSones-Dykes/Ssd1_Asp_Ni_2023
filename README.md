[![DOI](https://zenodo.org/badge/288707159.svg)](https://zenodo.org/badge/latestdoi/288707159)

# Ssd1_CRACanalysis_2020

This repository contains extended data for the manuscript:

> Yeast Ssd1 is a non-enzymatic member of the RNase II family with an alternative RNA recognition interface.
> Rosemary A. Bayne, Uma Jayachandran, Aleksandra Kasprowicz, Stefan Bresson, David Tollervey, Edward W. J. Wallace, and Atlanta G. Cook.
> Nucleic Acids Research, 2021 [doi: 10.1093/nar/gkab615](https://doi.org/10.1093/nar/gkab615)

This repository concentrates on analysis of CRAC data measuring the RNA-binding of yeast [Ssd1/YDR293C](https://www.yeastgenome.org/locus/S000002701), 2020.
The nearly-raw demultiplexed CRAC sequencing data is archived on Gene Expression Omnibus, accession [GSE159835](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159835).

This repository also contains some related analysis and supporting information that is not directly analysing the CRAC data, see the table of contents below.

Please address questions to Edward.Wallace@ed.ac.uk, at [the Wallace lab, University of Edinburgh](https://ewallace.github.io/).

Structural data is available on PDB, PDBid 7AM1.

# Overview - The pipeline for CRAC read processing

The main data here is CRAC, an Illumina short-read protocol related to HITS-CLIP, that detects RNA crosslinking sites for specific proteins. 
We analyse CRAC data using the python script, [src/CRAC_pipeline_SE_demult_dedup.py](src/CRAC_pipeline_SE_demult_dedup.py). 
This script uses the [ruffus](https://cgat-ruffus.readthedocs.io/) computational pipeline library to organise running of many 3rd-party tools.
This script was based on an earlier version by [Sander Granneman](http://sandergranneman.bio.ed.ac.uk/) for the paper:

> van Nues R, Schweikert G, de Leau E, Selega A, Langford A, Franklin R, Iosub I, Wadsworth P, Sanguinetti G, Granneman S. Kinetic CRAC uncovers a role for Nab3 in determining gene expression profiles during stress. Nat Commun. 2017 Apr 11;8(1):12. doi: 10.1038/s41467-017-00025-5. [PMID: 28400552](https://pubmed.ncbi.nlm.nih.gov/28400552/); PMCID: PMC5432031.

The script makes extensive use of his [pyCRAC software](https://git.ecdf.ed.ac.uk/sgrannem/pycrac):

> Webb S, Hector RD, Kudla G, Granneman S. PAR-CLIP data indicate that Nrd1-Nab3-dependent transcription termination regulates expression of hundreds of protein coding genes in yeast. Genome Biol. 2014 Jan 7;15(1):R8. doi: 10.1186/gb-2014-15-1-r8. [PMID: 24393166](https://pubmed.ncbi.nlm.nih.gov/24393166/); PMCID: PMC4053934.

This repository is "as self-contained as possible". Most files are kept locally, with the exception of the 
* .fastq data (very large, some is available on GEO)
* .novoindex aligner index (large, created from the genome sequence in `input_annotation`

The pipeline reads as inputs:

* one multiplexed .fastq file with all raw reads (available on request because GEO will not accept multiplexed files).
* annotation files, in `input_annotation`
* barcode files, in `input_barcodes`

The pipeline writes as outputs:

* demultiplexed reads in .fastq format (stored on GEO)
* various intermediate files of trimmed or unsorted reads (very large, not stored)
* reads aligned to the genome and sorted, in .bam format (very large, not stored in repository)
* read counts for each annotated transcript, in a .gff-like format
* pileups of local read density on the genome in .bedgraph format
* locations and hit counts for peaks in the CRAC data, in a .gff-like format

To run the pipeline as used here and in our paper, use this code:

```
python src/CRAC_pipeline_SE_demult_dedup.py \
  -f /homes/ewallac2/mount/datastore/wallace_rna/bigdata/2019/Ssd1CRAC-2019-01-15/20190114_Ssd1_CRAC.fastq  \
  -c input_annotation/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.lengths \
  --gtf input_annotation/Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.quotefix.gtf \
  --novoindex /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.novoindex \
  -p 5 \
  --name Ssd1_CRAC_demult_dedup_20190114_all \
  -a input_barcodes/3primeadapter.fasta \
  -b input_barcodes/Ssd1CRACSampleBarcodes.txt \
  --genelist input_annotation/Ssd1TargetGeneNamesOnly.txt \
  --genometab input_annotation/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.fa.tab \
  --transcriptgff input_annotation/abundant_verified_full-ORF_ypd_plus_other_fixed_UTR_length_transcripts.gff
```

Unfortunately, due to a bug that we have not been able to fix, this pipeline script needs to be run *twice* for all tasks to complete. We think the bug arises from the `ruffus` pipeline not waiting for `novoalign` to finish running, so that the first time the script runs all steps before genomic alignment, but those post-alignment return errors. Then the second time the pipeline runs, the post-alignment tasks run to completion (#overlyhonestmethods).


# Contents of this repository

Data and code is organised into subdirectories, whose contents are briefly mentioned in this section.
Each directory has its own `README.md` file that describes the contents in more detail.

We start with the contents that are directly part of to the CRAC data analysis.

## src

Source code that runs the pipeline, `CRAC_pipeline_SE_demult_dedup.py`.

## input_annotation

Annotation files (genome sequence, transcript maps) needed to run CRAC pipeline for data from S288C yeast.

## input_barcodes

Sequencing adapters and barcodes needed to run CRAC pipeline.

## Ssd1_CRAC_demult_dedup_20190114_all

This contains all the outputs generated by `src/CRAC_pipeline_SE_demult_dedup.py`.

## input_targets

Published data on Ssd1 targets from other studies.
Because these data are relatively large and publicly available, we did not directly include them in the repository.

They can be downloaded from the links provided in `input_targets/README.md`.

## rmarkdown

Analysis scripts in Rmarkdown (.Rmd) format, that produce many figures in the manuscript.
These run *after* the CRAC pipeline has completed.

## figure_out

Figures that are used for the manuscript, generated by scripts in `rmarkdown`.

## results

Tabular data and other results generated by scripts in directory `rmarkdown`. 


# Contents continued, that are not directly part of the CRAC analysis

## motif_conservation

Motif conservation analysis in alternative fungi.
This analysis takes the `CNYTCNYT` motif found via CRAC analysis,as calculate by code in `rmarkdown/peakanalysis_Ssd1.html`.
It then searches for that in 5'UTR regulatory sequences in _C. albicans_, _S. pombe_, and _A. fumigatus_, and reports gene ontology analysis.

## motifs_SUNfamily

Analysis of S. cerevisiae SUN family (SIM1/UTH1/NCA3/SUN4) homologs in ascomycetes in _C. albicans_, _S. pombe_, _A. fumigatus_, and _N. crassa_, to look for conservation of Ssd1 regulation by `CNYTCNYT` and `CCAACT` motifs. 

## strain_plasmid_info

Information on yeast strain construction for Ssd1 analysis.
Including construction of the Ssd1-HTP strain that we used for the CRAC analysis, and other yeast strains used in the paper but not directly for CRAC analysis.

It focuses on:
* Ssd1 locus maps detailing sequence and edits
* Plasmid maps for guide RNA plasmids for CRISPR
