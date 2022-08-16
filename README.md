# Annelid functional genomics reveal the origins of bilaterian life cycles - Supplementary Data
This GitHub repository contains the supplementary data to the publication cited below.

## How to use and cite these files
All files are made publicly available and can be used for further research and other applications. However, if you use these resources in your work, we kindly ask you to cite our original publication.
> **Annelid functional genomics reveal the origins of bilaterian life cycles.**
> Francisco M. Martín-Zamora§, Yan Liang§, Kero Guynes, Allan M. Carrillo-Baltodano, Billie E. Davies, Rory D. Donellan, Yongkai Tan, Giacomo Moggioli, Océane Seudre, Martin Tran, Kate Mortimer, Nicholas M. Luscombe, Andreas Hejnol, Ferdinand Marlétaz, José M. Martín-Duran.
> bioRxiv 2022.02.05.479245; doi: https://doi.org/10.1101/2022.02.05.479245

## Author contact
- [José M. Martín-Durán](mailto:chema.martin@qmul.ac.uk) (senior co-author)
- [Ferdinand Marlétaz](mailto:f.marletaz@ucl.ac.uk) (senior co-author)
- [Francisco M. Martín-Zamora](mailto:f.m.martinzamora@qmul.ac.uk) (co-first author)
- [Yan Liang](mailto:y.liang@qmul.ac.uk) (co-first author)
- Lab website: [martinduranlab.com](https://www.martinduranlab.com)

## How to download these files
The files inside this repository are all compressed `.tar.gz` files stored in [GitHub Large File Storage (LFS)](https://git-lfs.github.com/) which will need to be decompressed. To view them, select the `.tar.gz` file you are interested in (e.g. [05-Chordin_orthology_assignment.tar.gz](05-Chordin_orthology_assignment.tar.gz)), click on *Download* and the file will be downloaded to your device.

> :warning: **Some files are very heavy and download might take time to begin:** :warning: Be patient, they are correctly stored in [GitHub Large File Storage (LFS)](https://git-lfs.github.com/) and your download will eventually start.

## Index of contents
- [01-Owenia_fusiformis_genome_assemblies.tar.gz](01-Owenia_fusiformis_genome_assemblies.tar.gz) contains the unmasked genome files in `.fasta`/`.fa` format of the two genome assemblies versions we report in the publication:
  - **Scaffold_level_assembly_unmasked.fa**, which is the first high-quality genome assembly that we created, upon which we have based our RNA-seq and ATAC-seq analyses, and all other analyses presented in the original publication (see above).
  - **Chromosome_level_assembly_unmasked.fa**, which is the chromosome-level genome assembly after HiC chromosome-scale scaffolding.
- [02-Genome_annotations.tar.gz](02-Genome_annotations.tar.gz) contains the genome annotation files with all gene models in `.gff3` format and the funcional annotation reports in `.xlsx` format for both *O. fusiformis* and *C. teleta*. Again, as described above for [01-Owenia_fusiformis_genome_assemblies.tar.gz](01-Owenia_fusiformis_genome_assemblies.tar.gz), there are two different annotation files for *O. fusiformis*, which match each of the genome asembly files.
  - **Capitella_teleta_genome_annotation.gff3** contains the new updated gene models of *C. teleta*.
  - **Capitella_teleta_functional_annotation.xlsx** is the functional annotation report of *C. teleta*.
  - **Owenia_fusiformis_genome_annotation_scaffold_assembly.gff3** contains the gene models of *O. fusiformis* corresponding to **Scaffold_level_assembly_unmasked.fa**.
  - **Owenia_fusiformis_genome_annotation_chromosome_assembly.gff3** contains the gene models of *O. fusiformis* corresponding to **Chromosome_level_assembly_unmasked.fa**.
  - **Owenia_fusiformis_functional_annotation_scaffold_assembly.xlsx** is the functional annotation report of *O. fusiformis* corresponding to **Owenia_fusiformis_genome_annotation_scaffold_assembly.gff3**.
  - **Owenia_fusiformis_functional_annotation_chromosome_assembly.xlsx** is the functional annotation report of *O. fusiformis* corresponding to **Owenia_fusiformis_genome_annotation_chromosome_assembly.gff3**.
- [03-Transposable_elements.tar.gz](03-Transposable_elements.tar.gz) contains the repetitive and transposable elements repertoires of both *O. fusiformis* and *C. teleta* in `.fasta`/`.fa` and `.gff3` formats, which were used to mask the unmasked genome files.
  - **Capitella_teleta_transposable_elements.fa** contains the sequences of the repetitive and transposable elements of *C. teleta*.
  - **Owenia_fusiformis_transposable_elements.fa** contains the sequences of the repetitive and transposable elements of *O. fusiformis*.
  - **Capitella_teleta_transposable_elements.gff3** contains the repetitive and transposable elements models of *C. teleta*.
  - **Owenia_fusiformis_transposable_elements.gff3**, contains the repetitive and transposable elements models of *O. fusiformis*, corresponding to  **Scaffold_level_assembly_unmasked.fa**.
- [04-Orthogroups.tar.gz](04-Orthogroups.tar.gz) contains the classification by orthogroups of the non-redundant genome-based gene models of the analysed metazoan lineages in `.tsv` format, from which single copy orthologs can be parsed. There are two different files, which were used for different analyses of the publication:
  - **Orthogroups_analysis1.tsv** contains the orthogroup classification of 22 metazoan lineages. This analysis corresponds to the gene family evolution analysis used for the results reported in Fig. 1 and Extended Data Fig. 1. 
  - **Orthogroups_analysis1.tsv** contains the orthogroup classification of 27 metazoan lineages. This analysis corresponds to the gene family evolution analysis used for the results reported in the rest of the publication.
- [05-Chordin_orthology_assignment.tar.gz](05-Chordin_orthology_assignment.tar.gz) contains the multiple sequence alignment (MSA) files in `.fasta`/`.fa` format used for the phylogenetic inference on Chordin sequences.
  - **Chordin_mother_alignment.fa** contains the "mother" MSA of 15 curated Chordin and Chordin-like protein sequences, an outgroup (human BMP endothelial regulator, BMPER), and the Chordin and Chordin-like sequences from *O. fusiformis*.
  - **Chordin_daughter_alignment_1.fa** is the "daughter" MSA #1, which contains all sequences from **Chordin_mother_alignment.fa** and all candidate annelid sequences that have a 10-residue or longer properly aligned fragment in the CHRD domain region of the MSA.
  - **Chordin_daughter_alignment_2.fa** is the "daughter" MSA #2, which contains all sequences from **Chordin_mother_alignment.fa** and all candidate annelid sequences that have a 10-residue or longer properly aligned fragment in any of the vWFC domain regions of the MSA.
- [06-Hox_orthology_assignment.tar.gz](06-Hox_orthology_assignment.tar.gz) contains the multiple sequence alignment (MSA) files in `.fasta`/`.fa` format used for the phylogenetic inference on *Hox* sequences. 
  - **Hox_alignment1.fa** contains the alignment with gaps as hyphens(`-`).
  - **Hox_alignment2.fa** contains the alignment with gaps as hyphens(` `).
- [07-Annelid_RNAseq.tar.gz](07-Annelid_RNAseq.tar.gz) contains the gene expression matrices in transcripts-per-million (TPM) and after DESeq2 normalisation of the developmental RNA-seq time courses of *O. fusiformis*, *C. teleta*, and *D. gyrociliatus* containing all transcripts (i.e., all isoforms per gene model), both split by biological replicates and averaged by developmental stage. For each species:
  - **<Species_name>_all_transcripts_TPM_replicates.txt** is the gene expression matrix of all transcripts in TPM split by replicates.
  - **<Species_name>_all_transcripts_TPM_average.txt** is the gene expression matrix of all transcripts in TPM averaged by developmental stage.
  - **<Species_name>_all_transcripts_DESeq2_replicates.txt** is the gene expression matrix of all transcripts in DESeq2 normalised values split by replicates.
  - **<Species_name>_all_transcripts_DESeq2_average.txt** is the gene expression matrix of all transcripts in DESeq2 normalised values averaged by developmental stage.
- [08-Metazoan_RNAseq.tar.gz](08-Metazoan_RNAseq.tar.gz) contains the gene expression matrices in transcripts-per-million (TPM) and after DESeq2 normalisation of the developmental RNA-seq time courses on the non-redundant genome-based gene models (i.e., with only one isoform, the longest, per gene model) of *A. queenslandica*, *B. lanceolatum*, *C. elegans*, *C. teleta*, *C. hemisphaerica*, *C. gigas*, *D. rerio*, *D. gyrociliatus*, *D. melanogaster*, *N. vectensis*, *O. fusiformis*, and *S. purpuratus*, both split by biological replicates and averaged by developmental stage. For each species:
  - **<Species_name>_longest_isoform_TPM_replicates.txt** is the gene expression matrix of the longest isoform of each gene in TPM split by replicates.
  - **<Species_name>_longest_isoform_TPM_average.txt** is the gene expression matrix of the longest isoform of each gene in TPM averaged by developmental stage.
  - **<Species_name>_longest_isoform_DESeq2_replicates.txt** is the gene expression matrix of the longest isoform of each gene in DESeq2 normalised values split by replicates.
  - **<Species_name>_longest_isoform_DESeq2_average.txt** is the gene expression matrix  of the longest isoform of each gene in DESeq2 normalised values averaged by developmental stage.
- [09-Metazoan_Hox_expression_RNAseq.tar.gz](09-Metazoan_Hox_expression_RNAseq.tar.gz) contains the gene expression matrices as described above for [08-Metazoan_Hox_expression_RNAseq.tar.gz](08-Metazoan_Hox_expression_RNAseq.tar.gz) but with the added *Hox* sequences for the species whose gene models were lacking any. These are *B. lanceolatum*, *C. teleta*, *D. rerio*, and *D. melanogaster*. For each stage, the file structure is the same as above, except for the fact there is not any TPM averaged gene expression matrices. 
- [10-Owenia_fusiformis_WGCNA.tar.gz](10-Owenia_fusiformis_WGCNA.tar.gz) and [11-Capitella_teleta_WGCNA.tar.gz](11-Capitella_teleta_WGCNA.tar.gz) contain the nodes and edges files of the full WGCNA network, as well as the nodes and edges files of the random 30% subset of the networks used for visualisation, for both *O. fusiformis* and *C. teleta*. For each `.tar.gz` file:
  - **<Species_name>_WGCNA_nodes.txt** is the nodes file of the full WGCNA network.
  - **<Species_name>_WGCNA_edges.txt** is the edges file of the full WGCNA network.
  - **<Species_name>_WGCNA_30%_nodes.txt** is the nodes file of the random subset of 30% nodes of the WGCNA network.
  - **<Species_name>_WGCNA_30%_edges.txt** is the edges file of the edges corresponding to the random subset of 30% nodes of the WGCNA network.
- [12-Adult_tissues_RNAseq.tar.gz](12-Adult_tissues_RNAseq.tar.gz) contains the gene expression matrices in transcripts-per-million (TPM) and after DESeq2 normalisation of the adult tissue RNA-seq samples of *O. fusiformis*.
  - **Adult_tissues_RNAseq_TPM.txt** is the gene expression matrix in TPM.
  - **Adult_tissues_RNAseq_DESeq2.txt** is the gene expression matrix in DESeq2 normalised values.
- [13-Owenia_fusiformis_ATACseq_coverage.tar.gz](13-Owenia_fusiformis_ATACseq_coverage.tar.gz) and [14-Capitella_teleta_ATACseq_coverage.tar.gz](14-Capitella_teleta_ATACseq_coverage.tar.gz) contain the ATAC-seq stage-specific coverage files in `.bigwig`/`.bw` format for both *O. fusiformis* and *C. teleta*.
- [15-Owenia_fusiformis_ATACseq_peaks.tar.gz](15-Owenia_fusiformis_ATACseq_peaks.tar.gz) and [16-Capitella_teleta_ATACseq_peaks.tar.gz](16-Capitella_teleta_ATACseq_peaks.tar.gz) which contains the ATAC-seq stage-specific peak sets and the consensus peak sets in `.bed` format for both *O. fusiformis* and *C. teleta*.

