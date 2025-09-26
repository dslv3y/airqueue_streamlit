```mermaid
graph TD
    tree[tree]
    tree --> tree__
    tree__[.]
    tree --> tree_Coverage
    tree_Coverage[Coverage]
    tree_Coverage --> tree_Coverage_CoverageAnalysis_md
    tree_Coverage_CoverageAnalysis_md[CoverageAnalysis.md]
    tree_Coverage --> tree_Coverage_window_based_coverage_md
    tree_Coverage_window_based_coverage_md[window_based_coverage.md]
    tree --> tree_Epigenetics
    tree_Epigenetics[Epigenetics]
    tree_Epigenetics --> tree_Epigenetics_ChIP_seq
    tree_Epigenetics_ChIP_seq[ChIP_seq]
    tree_Epigenetics --> tree_Epigenetics_ChIPseq_Geneprofile_md
    tree_Epigenetics_ChIPseq_Geneprofile_md[ChIPseq_Geneprofile.md]
    tree_Epigenetics --> tree_Epigenetics_PeakCallingMacs2_md
    tree_Epigenetics_PeakCallingMacs2_md[PeakCallingMacs2.md]
    tree_Epigenetics --> tree_Epigenetics_README_md
    tree_Epigenetics_README_md[README.md]
    tree --> tree_Evolutionary_theory
    tree_Evolutionary_theory[Evolutionary_theory]
    tree_Evolutionary_theory --> tree_Evolutionary_theory_Foundations_of_evolutionary_biology_PhD_course
    tree_Evolutionary_theory_Foundations_of_evolutionary_biology_PhD_course[Foundations_of_evolutionary_biology_PhD_course]
    tree_Evolutionary_theory --> tree_Evolutionary_theory_Evolutionary_quantitative_genetics_md
    tree_Evolutionary_theory_Evolutionary_quantitative_genetics_md[Evolutionary_quantitative_genetics.md]
    tree_Evolutionary_theory --> tree_Evolutionary_theory_Genetic_drift_md
    tree_Evolutionary_theory_Genetic_drift_md[Genetic_drift.md]
    tree_Evolutionary_theory --> tree_Evolutionary_theory_Interactions_between_selection_and_drift_md
    tree_Evolutionary_theory_Interactions_between_selection_and_drift_md[Interactions_between_selection_and_drift.md]
    tree_Evolutionary_theory --> tree_Evolutionary_theory_Pop_structure_and_nonrandom_mating_md
    tree_Evolutionary_theory_Pop_structure_and_nonrandom_mating_md[Pop_structure_and_nonrandom_mating.md]
    tree_Evolutionary_theory --> tree_Evolutionary_theory_README_md
    tree_Evolutionary_theory_README_md[README.md]
    tree_Evolutionary_theory --> tree_Evolutionary_theory_Single_locus_natural_selection_md
    tree_Evolutionary_theory_Single_locus_natural_selection_md[Single_locus_natural_selection.md]
    tree_Evolutionary_theory --> tree_Evolutionary_theory_Two_locus_evolution_md
    tree_Evolutionary_theory_Two_locus_evolution_md[Two_locus_evolution.md]
    tree --> tree_Gene_expression_analysis
    tree_Gene_expression_analysis[Gene_expression_analysis]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_Bulk
    tree_Gene_expression_analysis_Bulk[Bulk]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_DESeq2_md
    tree_Gene_expression_analysis_DESeq2_md[DESeq2.md]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_Differential_splicing_analysis
    tree_Gene_expression_analysis_Differential_splicing_analysis[Differential_splicing_analysis]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_Differential_splicing_analysis_md
    tree_Gene_expression_analysis_Differential_splicing_analysis_md[Differential_splicing_analysis.md]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_Pie_AS_events_R
    tree_Gene_expression_analysis_Pie_AS_events_R[Pie_AS_events.R]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_Volcano_MA_plot_brain_R
    tree_Gene_expression_analysis_Volcano_MA_plot_brain_R[Volcano_MA_plot-brain.R]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_Volcano_MA_plot_R
    tree_Gene_expression_analysis_Volcano_MA_plot_R[Volcano_MA_plot.R]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_split_file_R
    tree_Gene_expression_analysis_split_file_R[split_file.R]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_Expression_analysis_200bp_windows_no_annotation_md
    tree_Gene_expression_analysis_Expression_analysis_200bp_windows_no_annotation_md[Expression_analysis_200bp_windows_no_annotation.md]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_KALLISTO_md
    tree_Gene_expression_analysis_KALLISTO_md[KALLISTO.md]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_KallistoSleuth_md
    tree_Gene_expression_analysis_KallistoSleuth_md[KallistoSleuth.md]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_NormalizerDE_md
    tree_Gene_expression_analysis_NormalizerDE_md[NormalizerDE.md]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_README_md
    tree_Gene_expression_analysis_README_md[README.md]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_SLEUTH_md
    tree_Gene_expression_analysis_SLEUTH_md[SLEUTH.md]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_StarDEseq2_md
    tree_Gene_expression_analysis_StarDEseq2_md[StarDEseq2.md]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_Transposable_elements
    tree_Gene_expression_analysis_Transposable_elements[Transposable_elements]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_TE_expression_md
    tree_Gene_expression_analysis_TE_expression_md[TE_expression.md]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_WGCNA_singlesample_md
    tree_Gene_expression_analysis_WGCNA_singlesample_md[WGCNA_singlesample.md]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_edgeR_DE_AS_md
    tree_Gene_expression_analysis_edgeR_DE_AS_md[edgeR_DE_AS.md]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_normalisation_expression_ipynb
    tree_Gene_expression_analysis_normalisation_expression_ipynb[normalisation_expression.ipynb]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_Single_cell
    tree_Gene_expression_analysis_Single_cell[Single_cell]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_DropSeq_cookbook_md
    tree_Gene_expression_analysis_DropSeq_cookbook_md[DropSeq_cookbook.md]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_Dropseq_10x_Salmon_md
    tree_Gene_expression_analysis_Dropseq_10x_Salmon_md[Dropseq_10x_Salmon.md]
    tree_Gene_expression_analysis --> tree_Gene_expression_analysis_pseudobulk_deseq2_analysis_for_single_cell_data_md
    tree_Gene_expression_analysis_pseudobulk_deseq2_analysis_for_single_cell_data_md[pseudobulk_deseq2_analysis_for_single_cell_data.md]
    tree --> tree_General_bioinformatics
    tree_General_bioinformatics[General_bioinformatics]
    tree_General_bioinformatics --> tree_General_bioinformatics_BasicSLURM_md
    tree_General_bioinformatics_BasicSLURM_md[BasicSLURM.md]
    tree_General_bioinformatics --> tree_General_bioinformatics_Instructions_md
    tree_General_bioinformatics_Instructions_md[Instructions.md]
    tree_General_bioinformatics --> tree_General_bioinformatics_README_md
    tree_General_bioinformatics_README_md[README.md]
    tree_General_bioinformatics --> tree_General_bioinformatics_SRAdownload_md
    tree_General_bioinformatics_SRAdownload_md[SRAdownload.md]
    tree_General_bioinformatics --> tree_General_bioinformatics_fasta_manipulations_md
    tree_General_bioinformatics_fasta_manipulations_md[fasta_manipulations.md]
    tree_General_bioinformatics --> tree_General_bioinformatics_no_password_ssh_md
    tree_General_bioinformatics_no_password_ssh_md[no_password_ssh.md]
    tree_General_bioinformatics --> tree_General_bioinformatics_tips_and_tricks_md
    tree_General_bioinformatics_tips_and_tricks_md[tips_and_tricks.md]
    tree --> tree_Genome_or_transcriptome_assembly
    tree_Genome_or_transcriptome_assembly[Genome_or_transcriptome_assembly]
    tree_Genome_or_transcriptome_assembly --> tree_Genome_or_transcriptome_assembly_BUSCO_md
    tree_Genome_or_transcriptome_assembly_BUSCO_md[BUSCO.md]
    tree_Genome_or_transcriptome_assembly --> tree_Genome_or_transcriptome_assembly_Genome_assembly
    tree_Genome_or_transcriptome_assembly_Genome_assembly[Genome_assembly]
    tree_Genome_or_transcriptome_assembly --> tree_Genome_or_transcriptome_assembly_Megahit_SOAPfusion_md
    tree_Genome_or_transcriptome_assembly_Megahit_SOAPfusion_md[Megahit_SOAPfusion.md]
    tree_Genome_or_transcriptome_assembly --> tree_Genome_or_transcriptome_assembly_README_md
    tree_Genome_or_transcriptome_assembly_README_md[README.md]
    tree_Genome_or_transcriptome_assembly --> tree_Genome_or_transcriptome_assembly_SOAPdenovo_GenomeAssembly_md
    tree_Genome_or_transcriptome_assembly_SOAPdenovo_GenomeAssembly_md[SOAPdenovo_GenomeAssembly.md]
    tree_Genome_or_transcriptome_assembly --> tree_Genome_or_transcriptome_assembly_GetORF_md
    tree_Genome_or_transcriptome_assembly_GetORF_md[GetORF.md]
    tree_Genome_or_transcriptome_assembly --> tree_Genome_or_transcriptome_assembly_IsoCollapse_md
    tree_Genome_or_transcriptome_assembly_IsoCollapse_md[IsoCollapse.md]
    tree_Genome_or_transcriptome_assembly --> tree_Genome_or_transcriptome_assembly_orf_to_strand_md
    tree_Genome_or_transcriptome_assembly_orf_to_strand_md[orf_to_strand.md]
    tree_Genome_or_transcriptome_assembly --> tree_Genome_or_transcriptome_assembly_ragtag_md
    tree_Genome_or_transcriptome_assembly_ragtag_md[ragtag.md]
    tree_Genome_or_transcriptome_assembly --> tree_Genome_or_transcriptome_assembly_transcriptome_assembly
    tree_Genome_or_transcriptome_assembly_transcriptome_assembly[transcriptome_assembly]
    tree_Genome_or_transcriptome_assembly --> tree_Genome_or_transcriptome_assembly_SOAPdenovo_trans_md
    tree_Genome_or_transcriptome_assembly_SOAPdenovo_trans_md[SOAPdenovo-trans.md]
    tree_Genome_or_transcriptome_assembly --> tree_Genome_or_transcriptome_assembly_de_novo_md
    tree_Genome_or_transcriptome_assembly_de_novo_md[de_novo.md]
    tree_Genome_or_transcriptome_assembly --> tree_Genome_or_transcriptome_assembly_trinity_md
    tree_Genome_or_transcriptome_assembly_trinity_md[trinity.md]
    tree --> tree_Miscellaneous
    tree_Miscellaneous[Miscellaneous]
    tree_Miscellaneous --> tree_Miscellaneous_README_md
    tree_Miscellaneous_README_md[README.md]
    tree --> tree_Orthology_and_phylogenomics
    tree_Orthology_and_phylogenomics[Orthology_and_phylogenomics]
    tree_Orthology_and_phylogenomics --> tree_Orthology_and_phylogenomics_Dn_ds
    tree_Orthology_and_phylogenomics_Dn_ds[Dn_ds]
    tree_Orthology_and_phylogenomics --> tree_Orthology_and_phylogenomics_README_md
    tree_Orthology_and_phylogenomics_README_md[README.md]
    tree_Orthology_and_phylogenomics --> tree_Orthology_and_phylogenomics_dNdS_md
    tree_Orthology_and_phylogenomics_dNdS_md[dNdS.md]
    tree_Orthology_and_phylogenomics --> tree_Orthology_and_phylogenomics_Phylogenetic_trees
    tree_Orthology_and_phylogenomics_Phylogenetic_trees[Phylogenetic_trees]
    tree_Orthology_and_phylogenomics --> tree_Orthology_and_phylogenomics_DNA_IQtree_md
    tree_Orthology_and_phylogenomics_DNA_IQtree_md[DNA_IQtree.md]
    tree_Orthology_and_phylogenomics --> tree_Orthology_and_phylogenomics_Protein_IQtree_md
    tree_Orthology_and_phylogenomics_Protein_IQtree_md[Protein_IQtree.md]
    tree --> tree_Population_and_quantitative_genomics
    tree_Population_and_quantitative_genomics[Population_and_quantitative_genomics]
    tree_Population_and_quantitative_genomics --> tree_Population_and_quantitative_genomics_Filtering
    tree_Population_and_quantitative_genomics_Filtering[Filtering]
    tree_Population_and_quantitative_genomics --> tree_Population_and_quantitative_genomics_LDpruning_plink2_md
    tree_Population_and_quantitative_genomics_LDpruning_plink2_md[LDpruning_plink2.md]
    tree_Population_and_quantitative_genomics --> tree_Population_and_quantitative_genomics_README_md
    tree_Population_and_quantitative_genomics_README_md[README.md]
    tree_Population_and_quantitative_genomics --> tree_Population_and_quantitative_genomics_Sample_and_individual_filtering_plink2_md
    tree_Population_and_quantitative_genomics_Sample_and_individual_filtering_plink2_md[Sample_and_individual_filtering_plink2.md]
    tree_Population_and_quantitative_genomics --> tree_Population_and_quantitative_genomics_Nonsynonymous_synonymous_diversity
    tree_Population_and_quantitative_genomics_Nonsynonymous_synonymous_diversity[Nonsynonymous_synonymous_diversity]
    tree_Population_and_quantitative_genomics --> tree_Population_and_quantitative_genomics_PiAPiS_SNPgenie_RNAseq_md
    tree_Population_and_quantitative_genomics_PiAPiS_SNPgenie_RNAseq_md[PiAPiS_SNPgenie_RNAseq.md]
    tree_Population_and_quantitative_genomics --> tree_Population_and_quantitative_genomics_Population_structure
    tree_Population_and_quantitative_genomics_Population_structure[Population_structure]
    tree_Population_and_quantitative_genomics --> tree_Population_and_quantitative_genomics_Admixture_md
    tree_Population_and_quantitative_genomics_Admixture_md[Admixture.md]
    tree_Population_and_quantitative_genomics --> tree_Population_and_quantitative_genomics_Extract_PCs_plink2_md
    tree_Population_and_quantitative_genomics_Extract_PCs_plink2_md[Extract_PCs_plink2.md]
    tree_Population_and_quantitative_genomics --> tree_Population_and_quantitative_genomics_SNP_calling
    tree_Population_and_quantitative_genomics_SNP_calling[SNP_calling]
    tree_Population_and_quantitative_genomics --> tree_Population_and_quantitative_genomics_CallSNPs_artMAP_md
    tree_Population_and_quantitative_genomics_CallSNPs_artMAP_md[CallSNPs_artMAP.md]
    tree --> tree_Proteomics
    tree_Proteomics[Proteomics]
    tree_Proteomics --> tree_Proteomics_README_md
    tree_Proteomics_README_md[README.md]
    tree --> tree_Quality_control_of_raw_reads
    tree_Quality_control_of_raw_reads[Quality_control_of_raw_reads]
    tree_Quality_control_of_raw_reads --> tree_Quality_control_of_raw_reads_FASTQC_md
    tree_Quality_control_of_raw_reads_FASTQC_md[FASTQC.md]
    tree_Quality_control_of_raw_reads --> tree_Quality_control_of_raw_reads_FastQC_qualitycontrol_md
    tree_Quality_control_of_raw_reads_FastQC_qualitycontrol_md[FastQC_qualitycontrol.md]
    tree_Quality_control_of_raw_reads --> tree_Quality_control_of_raw_reads_README_md
    tree_Quality_control_of_raw_reads_README_md[README.md]
    tree_Quality_control_of_raw_reads --> tree_Quality_control_of_raw_reads_TRIMMOMATIC_SE_PE_md
    tree_Quality_control_of_raw_reads_TRIMMOMATIC_SE_PE_md[TRIMMOMATIC_SE_PE.md]
    tree_Quality_control_of_raw_reads --> tree_Quality_control_of_raw_reads_Trimgalore_md
    tree_Quality_control_of_raw_reads_Trimgalore_md[Trimgalore.md]
    tree_Quality_control_of_raw_reads --> tree_Quality_control_of_raw_reads_Trimmomatic_md
    tree_Quality_control_of_raw_reads_Trimmomatic_md[Trimmomatic.md]
    tree --> tree_README_md
    tree_README_md[README.md]
    tree --> tree_Read_mapping
    tree_Read_mapping[Read_mapping]
    tree_Read_mapping --> tree_Read_mapping_README_md
    tree_Read_mapping_README_md[README.md]
    tree --> tree_Small_RNA
    tree_Small_RNA[Small_RNA]
    tree_Small_RNA --> tree_Small_RNA_Readme_md
    tree_Small_RNA_Readme_md[Readme.md]
    tree --> tree_tree_md
    tree_tree_md[tree.md]
    tree --> tree_28directories_84files
    tree_28directories_84files[28directories,84files]
