# Supplementary Information
This GitHub repository contains the supplementary files of the publication: XXXXX.

## How to use and cite these files
This repository belongs to Francisco M. Martín-Zamora, Yan Liang (as co-first authors of the associated publication), José María Martín-Durán (as senior co-author), and all other authors of the publication. All files are made publicly available and can be used for further research and other applications. However, if you use these resources in your work, we kindly ask you to cite our original publication: XXXXX.

## Author contact
- [José María Martín-Durán](mailto:chema.martin@qmul.ac.uk) (senior co-author)
- [Ferdinand Marlétaz](mailto:f.marletaz@ucl.ac.uk) (senior co-author)
- [Francisco M. Martín-Zamora](mailto:f.m.martinzamora@qmul.ac.uk) (co-first author)
- [Yan Liang](mailto:y.liang@qmul.ac.uk) (co-first author)
- Lab website: [martinduranlab.com](https://www.martinduranlab.com)

## How to download the files
The files inside this repository are of two types. Most of them are compressed `.tar.gz` files stored in [GitHub Large File Storage (LFS)](https://git-lfs.github.com/) which will need to be decompressed, whereas there are two additional files in `.xlsx` format which were uploaded straight to GitHub. To view them, select the `.tar.gz` or `.xslx` file you are interested in (e.g. [Chordin_alignments.tar.gz](Chordin_alignments.tar.gz)), click on *View Raw* or *Download*, and the file will be downloaded to your device. 

## Index of contents
- [ATACseq_coverage.bw.tar.gz](ATACseq_coverage.bw.tar.gz), which contains the ATAC-seq stage-specific coverage files in `.bigwig`/`bw` format.
- [ATACseq_peaks.bed.tar.gz](ATACseq_peaks.bed.tar.gz), which contains the ATAC-seq stage-specific peak set and the consensus peak set in `.bed` format.
- [Chordin_alignments.tar.gz](Chordin_alignments.tar.gz), which contains the multiple sequence alignment (MSA) files in `.fasta`/`.fa` format which were used for the phylogenetic analyses on Chordin. In particular:
  - **Chordin_SuppFile1.fa** is the "mother" MSA, which contains 15 curated Chordin and Chordin-like protein sequences, an outgroup (human BMP endothelial regulator), and the Chordin and Chordin-like sequences from *Owenia fusiformis*.
  - **Chordin_SuppFile2.fa** is the "daughter" MSA #1, which contains all sequences from **Chordin_SuppFile1.fa** and all candidate annelid sequences that contained a 10-residue or longer properly aligned fragment in the CHRD domain region of the MSA.
  - **Chordin_SuppFile3.fa** is the "daughter" MSA #2, which contains all sequences from **Chordin_Suppfile1.fa** and all candidate annelid sequences that contained a 10-residue or longer properly aligned fragment in any of the vWFC domain regions of the MSA.
- [Genome_files.tar.gz](Genome_files.tar.gz), which contains the unmasked genome files in `.fasta`/`.fa` format. There are two different versions:
  - **Genome_scaffold_level_unmasked.fa**, which is the first high-quality genome assembly that we created, upon which we have based our RNA-seq and ATAC-seq analyses, and all other analyses presented in the original publication (see above).
  - **Genome_chromosome_level_unmasked.fa**, which is the chromosome-level genome assembly after chromosome-scale scaffolding.
- [Genome_annotation.tar.gz](Genome_annotation.tar.gz), which contains the genome annotation files in `.gff3` format with all gene models. Again, as described above for [Genome_files.tar.gz](Genome_files.tar.gz), there are two different annotation files, which match each of the genome files.
  - **Genome_annotation_scaffold_level.gff3**, which is the annotation corresponding to **Genome_scaffold_level_unmasked.fa**.
  - **Genome_annotation_chromosome_level.gff3**, which is the annotation corresponding to **Genome_chromosome_level_unmasked.fa**.
- [Hox_alignments.tar.gz](Hox_alignments.tar.gz), which contains the multiple sequence alignment (MSA) files in `.fasta`/`.fa` format which were used for the phylogenetic analysis on Hox genes. The only difference **Hox_SuppFile2.fa** with **Hox_SuppFile1.fa** is that the latter contains gaps as hyphens(`-`), whereas the former denotes them as blank spaces (` `).
- [Owenia_annotation_chromosome_TrinoPanther.xlsx](Owenia_annotation_chromosome_TrinoPanther.xlsx), which contains an Excel spreadsheet in `.xslx` format containing the genome functional annotation of the chromosome-level genome annotation, i.e., **Genome_annotation_chromosome_level.gff3**. Likewise [Owenia_annotation_scaffold_TrinoPanther.xlsx](Owenia_annotation_scaffold_TrinoPanther.xlsx) contains the same for the chromosome-level genome annotation, i.e., **Genome_annotation_scaffold_level.gff3**.
- [Transposable_elements_annotation.tar.gz](Transposable_elements_annotation.tar.gz) contains two different files:
  - **Transposable_elements_SuppFile1.fa** contains the sequence of the annotated transposable elements in `.fasta`/`.fa` format.
  - **Transposable_elements_SuppFile2.fa** contains the scaffold-level annotation of transposable elements in `.gff3` format.
- [WGCNA_files.tar.gz](WGCNA_files.tar.gz), which contains the edges and nodes files corresponding to the gene co-expression network analysis we performed.
  - **WGCNA_SuppFile1_nodes.txt** and **WGCNA_SuppFile2_edges.txt** contain the nodes and edges, respectively, of the whole network.
  - **WGCNA_SuppFile3_30%\_nodes.txt** and **WGCNA_SuppFile4_30%\_edges.txt** contain the nodes and edges, respectively, of the random 30% subset we used for visualisation purposes.
