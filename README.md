# R_code_Master_thesis

### Div:
- Code for saving data is commented out
- Code for printing to console and plotting is not commented out
- Some of the code is from or inspired by Lars-Gustav Snipen’s lectures (BIN310, NMBU)




### Folders and files for each Rmarkdown file:

#### downloading_Refseq_genomes.Rmd
-	Sub_data/prokaryotes.csv

#### fastq_check_and_filter_assembled_contigs.Rmd
-	Sub_data/Staphylococcus_haemolyticus_LMGT4071_R1_subset.fq
-	Sub_data/Staphylococcus_haemolyticus_LMGT4071_R2_subset.fq
-	Sub_data/wt_contigs_subset.fa

#### gene_annotation_pan_genome_tree.Rmd
-	Sub_data/Prokka_results_GFF/*
-	Sub_data/gene_presence_absence.Rtab
-	Sub_data/refseq_assembly_tbl.txt
-	Sub_data/fasta_files/*
-	Sub_data/fasta_files/wt_contigs_filtered.fasta
-	Sub_data/Hmmer_results_txt/*
-	Sub_data/LORFs_fasta/*
-	Sub_data/long_list_all_CDS_prokka_filtered_LORFs.fasta
-	Sub_data/IPScan_results_GFF/*

#### Kallisto_index.Rmd
-	Sub_data/all.fasta.tbl_subset.RData
-	Sub_data/Super.GFF.tbl_kallisto_index_subset.RData
-	Sub_data/blastn_LMGT4071genesvsrefgenes_270322_subset.txt

#### DEA_Sh_DESeq.Rmd
-	Sub_data/meta_data_genes.RData
-	Sub_data/meta_data_samples.RData
-	Sub_data/Kallisto_results/*
-	Sub_data/annotation_blast_results.RData

#### operon_clustering.Rmd
-	Sub_data/num_reads_filtered_RNA_reads.txt
-	Sub_data/fastqc_sequence_length_distribution_plot.tsv
-	Sub_data/SPAdes.contigs_subset.fasta
-	Sub_data/carefulSPAdes.contigs_subset.fasta
-	Sub_data/metaSPAdes.contigs_subset.fasta
-	Sub_data/rnaSPAdes.contigs_subset.fasta
-	Sub_data/GFF_info_LMGT4071_genes.RData
-	Sub_data/Blast_genes_against_rnaseq_contigs/*

#### selected_DEgenes.Rmd
-	Sub_data/operon_groups_tbl.RData
-	Sub_data/GFF_info_LMGT4071_genes.RData
-	Sub_data/meta_data_genes.RData
-	Sub_data/Tables_with_DE_genes/*
-	Sub_data/gene_seqs_LMGT4071_subset.fasta

organizing_res_prokka_interproscan.Rmd
-	Sub_data/subset_all_prokka.RData
-	Sub_data/subset_all.LORFs.RData
-	Sub_data/subset_interproscan.RData
-	Sub_data/SuperGFF_fixed_seqid_290322.RData

