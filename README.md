```mermaid
graph TD
    root[root]
    root --> root_
    root_[.]
    root --> root_Coverage
    root_Coverage[Coverage]
    root_Coverage --> root_Coverage_CoverageAnalysismd
    root_Coverage_CoverageAnalysismd[CoverageAnalysis.md]
    root_Coverage --> root_Coverage_window_based_coveragemd
    root_Coverage_window_based_coveragemd[window_based_coverage.md]
    root --> root_Epigenetics
    root_Epigenetics[Epigenetics]
    root_Epigenetics --> root_Epigenetics_ChIP_seq
    root_Epigenetics_ChIP_seq[ChIP_seq]
    root_Epigenetics --> root_Epigenetics_ ChIPseq_Geneprofilemd
    root_Epigenetics_ ChIPseq_Geneprofilemd[│   ChIPseq_Geneprofile.md]
    root_Epigenetics --> root_Epigenetics_ PeakCallingMacs2md
    root_Epigenetics_ PeakCallingMacs2md[│   PeakCallingMacs2.md]
    root_Epigenetics --> root_Epigenetics_ READMEmd
    root_Epigenetics_ READMEmd[│   README.md]
    root_Epigenetics --> root_Epigenetics_READMEmd
    root_Epigenetics_READMEmd[README.md]
    root --> root_Evolutionary_theory
    root_Evolutionary_theory[Evolutionary_theory]
    root_Evolutionary_theory --> root_Evolutionary_theory_Foundations_of_evolutionary_biology_PhD_course
    root_Evolutionary_theory_Foundations_of_evolutionary_biology_PhD_course[Foundations_of_evolutionary_biology_PhD_course]
    root_Evolutionary_theory --> root_Evolutionary_theory_ Evolutionary_quantitative_geneticsmd
    root_Evolutionary_theory_ Evolutionary_quantitative_geneticsmd[│   Evolutionary_quantitative_genetics.md]
    root_Evolutionary_theory --> root_Evolutionary_theory_ Genetic_driftmd
    root_Evolutionary_theory_ Genetic_driftmd[│   Genetic_drift.md]
    root_Evolutionary_theory --> root_Evolutionary_theory_ Interactions_between_selection_and_driftmd
    root_Evolutionary_theory_ Interactions_between_selection_and_driftmd[│   Interactions_between_selection_and_drift.md]
    root_Evolutionary_theory --> root_Evolutionary_theory_ Pop_structure_and_nonrandom_matingmd
    root_Evolutionary_theory_ Pop_structure_and_nonrandom_matingmd[│   Pop_structure_and_nonrandom_mating.md]
    root_Evolutionary_theory --> root_Evolutionary_theory_ READMEmd
    root_Evolutionary_theory_ READMEmd[│   README.md]
    root_Evolutionary_theory --> root_Evolutionary_theory_ Single_locus_natural_selectionmd
    root_Evolutionary_theory_ Single_locus_natural_selectionmd[│   Single_locus_natural_selection.md]
    root_Evolutionary_theory --> root_Evolutionary_theory_ Two_locus_evolutionmd
    root_Evolutionary_theory_ Two_locus_evolutionmd[│   Two_locus_evolution.md]
    root_Evolutionary_theory --> root_Evolutionary_theory_READMEmd
    root_Evolutionary_theory_READMEmd[README.md]
    root --> root_Gene_expression_analysis
    root_Gene_expression_analysis[Gene_expression_analysis]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_Bulk
    root_Gene_expression_analysis_Bulk[Bulk]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_ DESeq2md
    root_Gene_expression_analysis_ DESeq2md[│   DESeq2.md]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_ Differential_splicing_analysis
    root_Gene_expression_analysis_ Differential_splicing_analysis[│   Differential_splicing_analysis]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_  Differential_splicing_analysismd
    root_Gene_expression_analysis_  Differential_splicing_analysismd[│   │   Differential_splicing_analysis.md]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_  Pie_AS_eventsR
    root_Gene_expression_analysis_  Pie_AS_eventsR[│   │   Pie_AS_events.R]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_  Volcano_MA_plotbrainR
    root_Gene_expression_analysis_  Volcano_MA_plotbrainR[│   │   Volcano_MA_plot-brain.R]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_  Volcano_MA_plotR
    root_Gene_expression_analysis_  Volcano_MA_plotR[│   │   Volcano_MA_plot.R]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_  split_fileR
    root_Gene_expression_analysis_  split_fileR[│   │   split_file.R]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_ Expression_analysis_200bp_windows_no_annotationmd
    root_Gene_expression_analysis_ Expression_analysis_200bp_windows_no_annotationmd[│   Expression_analysis_200bp_windows_no_annotation.md]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_ KALLISTOmd
    root_Gene_expression_analysis_ KALLISTOmd[│   KALLISTO.md]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_ KallistoSleuthmd
    root_Gene_expression_analysis_ KallistoSleuthmd[│   KallistoSleuth.md]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_ NormalizerDEmd
    root_Gene_expression_analysis_ NormalizerDEmd[│   NormalizerDE.md]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_ READMEmd
    root_Gene_expression_analysis_ READMEmd[│   README.md]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_ SLEUTHmd
    root_Gene_expression_analysis_ SLEUTHmd[│   SLEUTH.md]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_ StarDEseq2md
    root_Gene_expression_analysis_ StarDEseq2md[│   StarDEseq2.md]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_ Transposable_elements
    root_Gene_expression_analysis_ Transposable_elements[│   Transposable_elements]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_  TE_expressionmd
    root_Gene_expression_analysis_  TE_expressionmd[│   │   TE_expression.md]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_ WGCNA_singlesamplemd
    root_Gene_expression_analysis_ WGCNA_singlesamplemd[│   WGCNA_singlesample.md]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_ edgeR_DE_ASmd
    root_Gene_expression_analysis_ edgeR_DE_ASmd[│   edgeR_DE_AS.md]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_ normalisation_expressionipynb
    root_Gene_expression_analysis_ normalisation_expressionipynb[│   normalisation_expression.ipynb]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_READMEmd
    root_Gene_expression_analysis_READMEmd[README.md]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_Single_cell
    root_Gene_expression_analysis_Single_cell[Single_cell]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_DropSeq_cookbookmd
    root_Gene_expression_analysis_DropSeq_cookbookmd[DropSeq_cookbook.md]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_Dropseq_10x_Salmonmd
    root_Gene_expression_analysis_Dropseq_10x_Salmonmd[Dropseq_10x_Salmon.md]
    root_Gene_expression_analysis --> root_Gene_expression_analysis_pseudobulk_deseq2_analysis_for_single_cell_datamd
    root_Gene_expression_analysis_pseudobulk_deseq2_analysis_for_single_cell_datamd[pseudobulk_deseq2_analysis_for_single_cell_data.md]
    root --> root_General_bioinformatics
    root_General_bioinformatics[General_bioinformatics]
    root_General_bioinformatics --> root_General_bioinformatics_BasicSLURMmd
    root_General_bioinformatics_BasicSLURMmd[BasicSLURM.md]
    root_General_bioinformatics --> root_General_bioinformatics_Instructionsmd
    root_General_bioinformatics_Instructionsmd[Instructions.md]
    root_General_bioinformatics --> root_General_bioinformatics_READMEmd
    root_General_bioinformatics_READMEmd[README.md]
    root_General_bioinformatics --> root_General_bioinformatics_SRAdownloadmd
    root_General_bioinformatics_SRAdownloadmd[SRAdownload.md]
    root_General_bioinformatics --> root_General_bioinformatics_fasta_manipulationsmd
    root_General_bioinformatics_fasta_manipulationsmd[fasta_manipulations.md]
    root_General_bioinformatics --> root_General_bioinformatics_no_password_sshmd
    root_General_bioinformatics_no_password_sshmd[no_password_ssh.md]
    root_General_bioinformatics --> root_General_bioinformatics_tips_and_tricksmd
    root_General_bioinformatics_tips_and_tricksmd[tips_and_tricks.md]
    root --> root_Genome_or_transcriptome_assembly
    root_Genome_or_transcriptome_assembly[Genome_or_transcriptome_assembly]
    root_Genome_or_transcriptome_assembly --> root_Genome_or_transcriptome_assembly_BUSCOmd
    root_Genome_or_transcriptome_assembly_BUSCOmd[BUSCO.md]
    root_Genome_or_transcriptome_assembly --> root_Genome_or_transcriptome_assembly_Genome_assembly
    root_Genome_or_transcriptome_assembly_Genome_assembly[Genome_assembly]
    root_Genome_or_transcriptome_assembly --> root_Genome_or_transcriptome_assembly_ Megahit_SOAPfusionmd
    root_Genome_or_transcriptome_assembly_ Megahit_SOAPfusionmd[│   Megahit_SOAPfusion.md]
    root_Genome_or_transcriptome_assembly --> root_Genome_or_transcriptome_assembly_ READMEmd
    root_Genome_or_transcriptome_assembly_ READMEmd[│   README.md]
    root_Genome_or_transcriptome_assembly --> root_Genome_or_transcriptome_assembly_ SOAPdenovo_GenomeAssemblymd
    root_Genome_or_transcriptome_assembly_ SOAPdenovo_GenomeAssemblymd[│   SOAPdenovo_GenomeAssembly.md]
    root_Genome_or_transcriptome_assembly --> root_Genome_or_transcriptome_assembly_GetORFmd
    root_Genome_or_transcriptome_assembly_GetORFmd[GetORF.md]
    root_Genome_or_transcriptome_assembly --> root_Genome_or_transcriptome_assembly_IsoCollapsemd
    root_Genome_or_transcriptome_assembly_IsoCollapsemd[IsoCollapse.md]
    root_Genome_or_transcriptome_assembly --> root_Genome_or_transcriptome_assembly_READMEmd
    root_Genome_or_transcriptome_assembly_READMEmd[README.md]
    root_Genome_or_transcriptome_assembly --> root_Genome_or_transcriptome_assembly_orf_to_strandmd
    root_Genome_or_transcriptome_assembly_orf_to_strandmd[orf_to_strand.md]
    root_Genome_or_transcriptome_assembly --> root_Genome_or_transcriptome_assembly_ragtagmd
    root_Genome_or_transcriptome_assembly_ragtagmd[ragtag.md]
    root_Genome_or_transcriptome_assembly --> root_Genome_or_transcriptome_assembly_transcriptome_assembly
    root_Genome_or_transcriptome_assembly_transcriptome_assembly[transcriptome_assembly]
    root_Genome_or_transcriptome_assembly --> root_Genome_or_transcriptome_assembly_SOAPdenovotransmd
    root_Genome_or_transcriptome_assembly_SOAPdenovotransmd[SOAPdenovo-trans.md]
    root_Genome_or_transcriptome_assembly --> root_Genome_or_transcriptome_assembly_de_novomd
    root_Genome_or_transcriptome_assembly_de_novomd[de_novo.md]
    root_Genome_or_transcriptome_assembly --> root_Genome_or_transcriptome_assembly_trinitymd
    root_Genome_or_transcriptome_assembly_trinitymd[trinity.md]
    root --> root_Miscellaneous
    root_Miscellaneous[Miscellaneous]
    root_Miscellaneous --> root_Miscellaneous_READMEmd
    root_Miscellaneous_READMEmd[README.md]
    root --> root_Orthology_and_phylogenomics
    root_Orthology_and_phylogenomics[Orthology_and_phylogenomics]
    root_Orthology_and_phylogenomics --> root_Orthology_and_phylogenomics_Dn_ds
    root_Orthology_and_phylogenomics_Dn_ds[Dn_ds]
    root_Orthology_and_phylogenomics --> root_Orthology_and_phylogenomics_ READMEmd
    root_Orthology_and_phylogenomics_ READMEmd[│   README.md]
    root_Orthology_and_phylogenomics --> root_Orthology_and_phylogenomics_ dNdSmd
    root_Orthology_and_phylogenomics_ dNdSmd[│   dNdS.md]
    root_Orthology_and_phylogenomics --> root_Orthology_and_phylogenomics_Phylogenetic_trees
    root_Orthology_and_phylogenomics_Phylogenetic_trees[Phylogenetic_trees]
    root_Orthology_and_phylogenomics --> root_Orthology_and_phylogenomics_ DNA_IQtreemd
    root_Orthology_and_phylogenomics_ DNA_IQtreemd[│   DNA_IQtree.md]
    root_Orthology_and_phylogenomics --> root_Orthology_and_phylogenomics_ Protein_IQtreemd
    root_Orthology_and_phylogenomics_ Protein_IQtreemd[│   Protein_IQtree.md]
    root_Orthology_and_phylogenomics --> root_Orthology_and_phylogenomics_READMEmd
    root_Orthology_and_phylogenomics_READMEmd[README.md]
    root --> root_Population_and_quantitative_genomics
    root_Population_and_quantitative_genomics[Population_and_quantitative_genomics]
    root_Population_and_quantitative_genomics --> root_Population_and_quantitative_genomics_Filtering
    root_Population_and_quantitative_genomics_Filtering[Filtering]
    root_Population_and_quantitative_genomics --> root_Population_and_quantitative_genomics_ LDpruning_plink2md
    root_Population_and_quantitative_genomics_ LDpruning_plink2md[│   LDpruning_plink2.md]
    root_Population_and_quantitative_genomics --> root_Population_and_quantitative_genomics_ READMEmd
    root_Population_and_quantitative_genomics_ READMEmd[│   README.md]
    root_Population_and_quantitative_genomics --> root_Population_and_quantitative_genomics_ Sample_and_individual_filtering_plink2md
    root_Population_and_quantitative_genomics_ Sample_and_individual_filtering_plink2md[│   Sample_and_individual_filtering_plink2.md]
    root_Population_and_quantitative_genomics --> root_Population_and_quantitative_genomics_Nonsynonymous_synonymous_diversity
    root_Population_and_quantitative_genomics_Nonsynonymous_synonymous_diversity[Nonsynonymous_synonymous_diversity]
    root_Population_and_quantitative_genomics --> root_Population_and_quantitative_genomics_ PiAPiS_SNPgenie_RNAseqmd
    root_Population_and_quantitative_genomics_ PiAPiS_SNPgenie_RNAseqmd[│   PiAPiS_SNPgenie_RNAseq.md]
    root_Population_and_quantitative_genomics --> root_Population_and_quantitative_genomics_Population_structure
    root_Population_and_quantitative_genomics_Population_structure[Population_structure]
    root_Population_and_quantitative_genomics --> root_Population_and_quantitative_genomics_ Admixturemd
    root_Population_and_quantitative_genomics_ Admixturemd[│   Admixture.md]
    root_Population_and_quantitative_genomics --> root_Population_and_quantitative_genomics_ Extract_PCs_plink2md
    root_Population_and_quantitative_genomics_ Extract_PCs_plink2md[│   Extract_PCs_plink2.md]
    root_Population_and_quantitative_genomics --> root_Population_and_quantitative_genomics_READMEmd
    root_Population_and_quantitative_genomics_READMEmd[README.md]
    root_Population_and_quantitative_genomics --> root_Population_and_quantitative_genomics_SNP_calling
    root_Population_and_quantitative_genomics_SNP_calling[SNP_calling]
    root_Population_and_quantitative_genomics --> root_Population_and_quantitative_genomics_CallSNPs_artMAPmd
    root_Population_and_quantitative_genomics_CallSNPs_artMAPmd[CallSNPs_artMAP.md]
    root --> root_Proteomics
    root_Proteomics[Proteomics]
    root_Proteomics --> root_Proteomics_READMEmd
    root_Proteomics_READMEmd[README.md]
    root --> root_Quality_control_of_raw_reads
    root_Quality_control_of_raw_reads[Quality_control_of_raw_reads]
    root_Quality_control_of_raw_reads --> root_Quality_control_of_raw_reads_FASTQCmd
    root_Quality_control_of_raw_reads_FASTQCmd[FASTQC.md]
    root_Quality_control_of_raw_reads --> root_Quality_control_of_raw_reads_FastQC_qualitycontrolmd
    root_Quality_control_of_raw_reads_FastQC_qualitycontrolmd[FastQC_qualitycontrol.md]
    root_Quality_control_of_raw_reads --> root_Quality_control_of_raw_reads_READMEmd
    root_Quality_control_of_raw_reads_READMEmd[README.md]
    root_Quality_control_of_raw_reads --> root_Quality_control_of_raw_reads_TRIMMOMATIC_SE_PEmd
    root_Quality_control_of_raw_reads_TRIMMOMATIC_SE_PEmd[TRIMMOMATIC_SE_PE.md]
    root_Quality_control_of_raw_reads --> root_Quality_control_of_raw_reads_Trimgaloremd
    root_Quality_control_of_raw_reads_Trimgaloremd[Trimgalore.md]
    root_Quality_control_of_raw_reads --> root_Quality_control_of_raw_reads_Trimmomaticmd
    root_Quality_control_of_raw_reads_Trimmomaticmd[Trimmomatic.md]
    root --> root_READMEmd
    root_READMEmd[README.md]
    root --> root_Read_mapping
    root_Read_mapping[Read_mapping]
    root_Read_mapping --> root_Read_mapping_READMEmd
    root_Read_mapping_READMEmd[README.md]
    root --> root_Small_RNA
    root_Small_RNA[Small_RNA]
    root_Small_RNA --> root_Small_RNA_Readmemd
    root_Small_RNA_Readmemd[Readme.md]
    root --> root_treemd
    root_treemd[tree.md]
    root --> root_28 directories, 84 files
    root_28 directories, 84 files[28 directories, 84 files]
