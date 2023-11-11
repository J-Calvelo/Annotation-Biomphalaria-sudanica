# Annotation-Biomphalaria-sudanica
Bioinformatic pipline utiliced in XXXXXX for the annotation of Biomphalaria sudanica.

## Publication_Bsudanica_annotation.sh
The script is subdevided into several blocks of analysis and/or file preparation in order to facilitate re-runs. Which blocks are ran at a time can be controlled by setting to true variables located from lines 7 to 52 (all written in all caps). Other important variables for the run are located betwee lines 59 to 113. Some will be discussed now while others more specific will be mentioned as they become relevant. Input files that are not third party databeses are included in: Other_Input_Files

### General variables:
1) species_name: Name for the species being annotated
2) genome_file: Full path to the genome assembly in fasta format. File available on the NCBI: SRRXXXXXX
3) work_dir: Full path to the working directory (needs to exists
4) threads: Number of computer cores available for the analysis. Used whenever possible.

### 1) CLEAN_RUN
Setting this variable to TRUE will remove all previous runs. Use only to start from a clean state:
#### Required software: 
- NA
#### Other Variables:
- NA
 
### 2) TEA_EARLGREY_HOT
Setting this variable to TRUE will run the earlGrey pipeline for repeat annotation
#### Required software: 
- EarlGrey (https://github.com/google/EarlGrey)
#### Other Variables: 
- earlgrey_instalation: Full path to EarlGrey

### 3) TEA_EARLGREY_SUMMARY
Generates a custom summary of EarlGrey results
#### Required software: 
- NA
#### Other Variables:
- NA

### 4) CLEAN_PACBIO_READS
Setting this variable to TRUE will run will retrieve full sequenced isoforms from the PACBIO CSS reads based on the presence of 5' and 3' isoforms
Setting CLEAN_PACBIO_READS_RM_TM to TRUE will also remove intermediary files.
#### Required software: 
- lima (https://github.com/PacificBiosciences/barcoding)
- isoseq3 (https://github.com/PacificBiosciences/IsoSeq)
- bam2fastq (https://github.com/jts/bam2fastq)
- seqkit (https://github.com/shenwei356/seqkit)
#### Other Variables:
- raw_pacbio_reads: Full path to the PACBIO CCS reads file (BAM Format). Data available on the NCBI: SRRXXXXXX
- pacbio_primers: Full path to a file with the adapter sequences of the PACBIO CCS readss in FASTA Format. File: Pacbio_adapters.fasta
  
### 5) CLEAN_ILLUMINA_READS 
Setting this variable to TRUE will remove adapter sequences and low quality bases from the paired end Illumina reads using trimmomatic (SLIDINGWINDOW:5:20 MINLEN:25)
#### Required software: 
- trimmomatic (https://github.com/usadellab/Trimmomatic)
- seqkit (https://github.com/shenwei356/seqkit)
#### Other Variables:
- raw_illumina_reads_dir: Full path to the folder with the Illumina read files. Data available on the NCBI: SRRXXXXXX
- raw_illumina_1: Identifier for the Left reads, located at the end of the file 
- raw_illumina_2: Identifier for the Left reads, located at the end of the file
- Illumina_adapters: Full path to a fasta file with Illumina adapters. File TruSeq3-PE-2.fa, included as trimmomatic's installation.
Note: Sample names are identified running the comand:
ls $raw_illumina_reads_dir"/"*$raw_illumina_1 | sed "s/.*\///" | sed "s/$raw_illumina_1//"

### 6) PACBIO_MAP
Setting this variable to TRUE will map the PACBIO reads to the genome
#### Required software: 
- minimap2 (https://github.com/lh3/minimap2)
- samtools (https://github.com/samtools/samtools)
#### Other Variables:
- NA

### 7) STAR_MAP
Setting this variable to TRUE will map the Illuimina reads to the genome
#### Required software: 
- STAR (https://github.com/alexdobin/STAR/tree/master)
#### Other Variables:
- star_index_num: Value recomended for the genome assembly size for --genomeChrBinNbits and --genomeSAindexNbases 

### 8) STRINGTIE_TRANS
Setting this variable to TRUE will combine the read mappings of Illumina and PACBIO reads to annotate the genome using Stringtie v2.2
Earlier versions of stringtie will not work
#### Required software: 
- stringtie (https://github.com/gpertea/stringtie)
- gffread (https://github.com/gpertea/gffread)
#### Other Variables:
- NA

Note: This block has the function kill_fai that will remove any index file of the genome that could cause problems with gffread
Command: rm $genome_file".fai"

### 9) FIND_MITOCHONDRIA 
Setting this variable to TRUE will try to identify the mitochondrial genome in the genome assembly based on blast searches against known mitochondrial gene sequences of other biomphalaria. 
Setting FIND_MITOCHONDRIA_RM_TM to TRUE will also remove intermediary files.
#### Required software: 
- BLAST+ (http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs)
#### Other Variables:
- mito_genes: Full path to a fasta file with known mitochondrial genes. File: Mitochondria_Coding_genes.fa

Note: Other blocks working with the mitochondria will assume that it's sequence was retrieved and annotated by means not included on this pipeline.

### 10) MAP_MITOCHONDRIA_PACBIO
Setting this variable to TRUE will map PACBIO reads to the mitochondrial genome and count reads asigned to the different features. 
#### Required software: 
- minimap2 (https://github.com/lh3/minimap2)
- samtools (https://github.com/samtools/samtools)
- htseq-count (https://github.com/htseq/htseq)
#### Other Variables:
- mito_genome: Full path to a fasta file with the conting/scaffold with the mitochondrial genome. File: Mitocondrial_genome.fa
- mito_gff: Full path to a GFF file with annotated mitochondrial genome. File: Mitocondrial_htseq.gff

### 11) MAP_MITOCHONDRIA_ILLUMINA
Setting this variable to TRUE will map paired Illumina reads to the mitochondrial genome and count reads asigned to the different features. 
#### Required software: 
- STAR (https://github.com/alexdobin/STAR/tree/master)
- samtools (https://github.com/samtools/samtools)
- htseq-count (https://github.com/htseq/htseq)
#### Other Variables:
- mito_genome: Full path to a fasta file with the conting/scaffold with the mitochondrial genome. File: Mitocondrial_genome.fa
- mito_htseq_illumina_gff: Full path to a GFF file with annotated mitochondrial genome (modified for convinience). File: Mitocondrial_htseq_Illumina.gff

Note: For convininence mito_htseq_illumina_gff has every feature identidied as a "gene" (third column of the gff) to better identify reads that overlap different features.

### 12) MAP_MITOCHONDRIA_PACBIO_NEW_ORIGIN
Setting this variable to TRUE will run the same analysis than MAP_MITOCHONDRIA_PACBIO  but over a genome file with a different origin coordinate
#### Required software: 
- minimap2 (https://github.com/lh3/minimap2)
- samtools (https://github.com/samtools/samtools)
- htseq-count (https://github.com/htseq/htseq)
#### Other Variables:
- new_ori_genome: Full path to a fasta file with the conting/scaffold with the mitochondrial genome. File: Mitocondrial_New_Origin_genome.fasta
- new_ori_gff: Full path to a GFF file with annotated mitochondrial genome. File: Mitocondrial_Mito_New_Origin_Htseq.gff

Note: This was done for quality control

### 13) MAP_MITOCHONDRIA_ILLUMINA_NEW_ORIGIN
Setting this variable to TRUE will run the same analysis than MAP_MITOCHONDRIA_ILLUMINA but over a genome file with a different origin coordinate
#### Required software: 
- STAR (https://github.com/alexdobin/STAR/tree/master)
- samtools (https://github.com/samtools/samtools)
- htseq-count (https://github.com/htseq/htseq)
#### Other Variables:
- new_ori_genome: Full path to a fasta file with the conting/scaffold with the mitochondrial genome. File: Mitocondrial_New_Origin_genome.fasta
- new_ori_htseq_illumina_gff: Full path to a GFF file with annotated mitochondrial genome. File: Mitocondrial_Mito_New_Origin_Htseq.gff

Note: This was done for quality control

### 14) ORF_PROT_SEQ_1
Setting this variable to TRUE will run the first half of Transdecoder pipeline for the identification of ORFs in the transcriptome predicted with STRINGTIE_TRANS
#### Required software: 
- Transdecoder (https://github.com/TransDecoder/TransDecoder)
- BLAST+ (http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs)
- hmmer (https://github.com/EddyRivasLab/hmmer/tree/master)
#### Other Variables:
- uniprot: Full path to a fasta file to reference proteins from Uniprot's Swiss-Prot database. Check main manuscript for details.
- transdecoder_evalue: Cut-off evalue used on the blastp search between the predicted proteins by TransDecoder.LongOrfs and the Uniprot reference
- pfam: Full path to hmmer frofiles from Pfam-A. Check manuscript for details.
- exclude_mitochondria: ID of the contig or scaffold that contains the mitochondrial genome. Any transcript assigned to it will not be analysed by trans-decoder

### 15) ORF_PROT_SEQ_2
Setting this variable to TRUE will run the second half of Transdecoder pipeline for the identification of ORFs in the transcriptome predicted with STRINGTIE_TRANS
#### Required software: 
- Transdecoder (https://github.com/TransDecoder/TransDecoder)
- gffread (https://github.com/gpertea/gffread)
#### Other Variables:
- NA
Note: This block has the function kill_fai that will remove any index file of the genome that could cause problems with gffread
Command: rm $genome_file".fai"
In addition it will remove temporary files in the working directory generated by TransDecoder.Predict and move the output with the commands:
  rm pipeliner.*
  mv PolyA_Transcripts_work* $work_dir"/TransDecoder/Predict"

### 16) BUSCO_TESTS
Setting this variable to TRUE will evaluate how complete are the genome assembly and the predicted transcriptome and proteome using BUSCO (-l mollusca). Which of the three is to be run can be adjusted with variables: BUSCO_GENOME, BUSCO_RNA and BUSCO_PROT_ALL
#### Required software: 
- Busco (http://busco.ezlab.org/)
#### Other Variables:
-NA

### 17) RUN_INTERPROT
Setting this variable to TRUE will annotate the predicted proteins by identifying protein domain signatures 
#### Required software: 
- interproscan (https://github.com/ebi-pf-team/interproscan)
#### Other Variables:
- interprot_path: Full path to the folder of the Interpro installation

### 18) RUN_EGGNOG
Setting this variable to TRUE will assign GO terms to the predicted proteins using eggnog mapper
#### Required software:
- eggnog-mapper (https://github.com/eggnogdb/eggnog-mapper/tree/master)
#### Other Variables:
- eggnog_data: Full path to the folder with eggnog-mapper reference data

### 19) RUN_TRNA
Setting this variable to TRUE will identify tRNA loci in the genome
#### Required software:
- tRNAscan-SE (https://github.com/UCSC-LoweLab/tRNAscan-SE)
#### Other Variables:
- NA

### 20) RUN_RIBOSOME
Setting this variable to TRUE will identify rRNA loci in the genome
#### Required software:
- barrnap (https://github.com/tseemann/barrnap)
#### Other Variables:
- NA

### 21) SEARCH_SIGNALS
Setting this variable to TRUE will identify signal peptides and mitochondrial taegeting peptides among the predicted proteins
#### Required software:
- signalp6 (https://github.com/fteufel/signalp-6.0)
- targetP (https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=targetp&version=2.0&packageversion=2.0&platform=Linux)
#### Other Variables:
- targetp_file. Full path to the TargetP installation folder 

### 22) SEARCH_SECRETOMEP
Setting this variable to TRUE will search for additional secreted proteins using SecretomeP
#### Required software:
- seqkit (https://github.com/shenwei356/seqkit)
- secretomeP (https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=secretomep&version=1.0&packageversion=1.0h&platform=Linux)
 #### Other Variables:
- secretomep: full path for the secretomep folder instalation. 

### 23) SIGNAL_RELOCATION_TABLE
Setting this variable to TRUE will reformat SecretomeP output in line with the rest of the pipeline
#### Required software:
- NA
#### Other Variables:
- NA

### 24) Run DeepTMHMM
Setting this variable to TRUE will run DeepTMHMM and identify transmembrane domains on the predicted proteins.
#### Required software:
- DeepTMHMM (https://dtu.biolib.com/DeepTMHMM)
#### Other Variables:
- NA

### 25) TABLE_TMHMM
Setting this variable to TRUE will reformat DeepTMHMM results, combine them with SEARCH_SIGNALS and SEARCH_SECRETOMEP and do a tentative clasification on the protein based on the presense/absence of signal peptides and transmembrane results.
#### Required software:
- DeepTMHMM (https://dtu.biolib.com/DeepTMHMM)
#### Other Variables:
- NA

### 26) SAM_IGSF_SEARCH
Setting this variable to TRUE will identify IgSF domains based on custom hmmer profiles. 
#### Required software:
- seqkit (https://github.com/shenwei356/seqkit)
- hmmer (https://github.com/EddyRivasLab/hmmer/tree/master)
#### Other Variables:
- igsf_hmmer_profile: Full path to a file with custom IgsF domains (IgSF1 and IgSF2). Originally designed in:  https://doi.org/10.1371/journal.pntd.0008780. File Bg_IgSF.hmm

### 27) PREPARE_ORTHOFINDER
Setting this variable to TRUE will prepare the input files to run Orthofinder. It requires the specification of several folders with input data. The longest reported isoform of each nuclear gene will be selected.
#### Required software:
- Emboss (http://emboss.open-bio.org)
- BLAST+ (http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs)
- seqkit (https://github.com/shenwei356/seqkit)
#### Other Variables:
- orthofinder_folder: Full path of the folder where the Orthofinder is going to be run
- other_species_orthofinder: Full path to a folder with the cDNA where the Orthofinder is going to be run. Files inside must contain their cDNA and end in ".cds"
- other_species_mitochondria: Full path to a file with mitochondrial sequences of the target species. File Other_Species_Mitochondria.fasta
- dash_ids: List separated by white species where Isoforms IDs fields are separated by "-" instead of "."

### 28) ORTHOFINDER
Setting this variable to TRUE will run Orthofinder with default parameters
#### Required software:
- Orthofinder (https://github.com/davidemms/OrthoFinder)
#### Other Variables:
- orthofinder_folder: Full path of the folder where the Orthofinder is going to be run

### 29) ORTHOFINDER_EGGNOG
Setting this variable to TRUE will run EGGNOG for all sequences anlaysed by Orthofinder
#### Required software:
- eggnog-mapper (https://github.com/eggnogdb/eggnog-mapper/tree/master)
- seqkit (https://github.com/shenwei356/seqkit)
#### Other Variables:
- orthofinder_folder: Full path of the folder where the Orthofinder is going to be run
- eggnog_data: Full path to the folder with eggnog-mapper reference data

### 30) ORTHOFINDER_EGGNOG_FOR_PHO
Setting this variable to TRUE will reformat ORTHOFINDER_EGGNOG and asign go terms to each HOG by assuming that if a GO term was assigned to one member it applies to the entire group.
#### Required software:
- NA
#### Other Variables:
- orthofinder_folder: Full path of the folder where the Orthofinder is going to be run

### 31) SEARCH_CREP_FREP_GREP
Setting this variable to TRUE will identify candidates members of the CREP, FREP, GREP and other related proteins based on the absence or presence of protein domains. And BLAST hits to known members these protein families.
#### Required software:
- seqkit (https://github.com/shenwei356/seqkit)
- BLAST+ (http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs)
#### Other Variables: 
- crep_frep_grep_blast_ident: Identity cutoff for the BLAST search against known members of these families.
- FBD_signatures: List separated by white spaces of interpro signatures associated with FBD domains
- C_lectin_signatures: List separated by white spaces of interpro signatures associated with c-lectin domains
- Galectin_signatures: List separated by white spaces of interpro signatures associated with galectin domains
- EGF_signatures_signatures: List separated by white spaces of interpro signatures associated with EGF domains
- Inmunoglobulin_signatures: List separated by white spaces of interpro signatures associated with Inmunoglobulin domains
- Dheilly_CREP_prot: Fasta file with reference CREP sequences taken from (doi: 10.1016/j.dci.2014.10.009). File: Dheilly_etal_2015_CREP.fa
- Dheilly_FREP_prot: Fasta file with reference FREP sequences taken from (doi: 10.1016/j.dci.2014.10.009). File: Dheilly_etal_2015_FREP.fa
- Dheilly_GREP_prot: Fasta file with reference GREP sequences taken from (doi: 10.1016/j.dci.2014.10.009). File:Dheilly_etal_2015_GREP.fa
- Lu_FREP_prot: Fasta file with reference FREP sequences taken from  (https://doi.org/10.1371/journal.pntd.0008780)

Note1: See main manuscript for the full description of the selecction procedure. Some of the sources of evidence were carried out as diagnostics.
Note2: The initial classification was then revised manually in a case by case.

### 32) CREP_FREP_PHYLOGENY
Setting this variable to TRUE will align CREP and FREP sequences, trim them automatically and callculate their phylogenetic relationships with iqtree. All sequences identified as Full will be included, with the exclusion of ones specified on a file (exclude_seq_manual)
#### Required software:
- seqkit (https://github.com/shenwei356/seqkit)
- Trimal (https://github.com/inab/trimal)
- iqtree (https://github.com/iqtree/iqtree2)
#### Other Variables: 
- exclude_seq_manual: Full path to a file with protein ids manually removed.

### 33) SUP_TAB_PAPER_CREP_FREP_GREP
Setting this variable to TRUE will summirize the results in a table 
#### Required software:
- NA
#### Other Variables: 
- NA

### 34) CAFE_EXPANSION_PREPARE
Setting this variable to TRUE will generate all input files for CAFE, from Orthofinder's results
#### Required software:
- make_ultrametric.py (one of orthofinde's utility scripts https://github.com/davidemms/OrthoFinder)
#### Other Variables: 
- orthofinder_folder: Full path of the folder where the Orthofinder is going to be run
- mini_ortho_min_Nspe: Minimun number of genes in the HOG to be considered.
- ultrametic_tree_root_age: Age of the tree's root in Million years
 
### 35) CAFE_EXPANSION_TEST_RUNS
Setting this variable to TRUE will conduct test runs for CAFE within a range of parameters and repetitions, or until one of them fails to compute.
#### Required software:
- cafe5 (https://github.com/hahnlab/CAFE5)
#### Other Variables: 
- orthofinder_folder: Full path of the folder where the Orthofinder is going to be run
- end_kvalue: maximun number of K to test
- max_cafe_runs_perK: Maximun number of repetitions per K to conduct
- max_number_iterations: Makimun number of itterations done by CAFE5

### 36) CAFE_EXPANSION_DEF_RUN
Setting this variable to TRUE will conduct one final run of CAFE5 that will be taken as definitive.
#### Required software:
- cafe5 (https://github.com/hahnlab/CAFE5)
- seqkit (https://github.com/shenwei356/seqkit)
#### Other Variables: 
- orthofinder_folder: Full path of the folder where the Orthofinder is going to be run
- max_number_iterations: Makimun number of itterations done by CAFE5
- Def_Kvalue: The best K-value identified in the test runs
- interest_nodes_file: Full path to a file with the node names separated by tabs: "Biomphalaria_glabrata_IM_GCA_025434175.1        <0>"

### 37) CAFE_EXPANSION_ANOT_EGGNOG
Setting this variable to TRUE will run EGGNOG for all sequences anlaysed by CAFE5.
#### Required software:
- eggnog-mapper (https://github.com/eggnogdb/eggnog-mapper/tree/master)
#### Other Variables:
- orthofinder_folder: Full path of the folder where the Orthofinder is going to be run
- eggnog_data: Full path to the folder with eggnog-mapper reference data

### 38) DATA_CROSS_POLYA_REPEATS
Setting this variable to TRUE will identify genes that overlap with repeated elements
#### Required software:
- NA
#### Other Variables:
- NA

### 39) SUMMARY
Setting this variable to TRUE will copy and move some result files for easier access
#### Required software:
- NA
#### Other Variables:
- NA

### 40) LOCATION_TABLE
Setting this variable to TRUE will compine the results of SignalP, TargetP, SecretomeP DeepTMHMM and InterproScan in a single table
#### Required software:
- NA
#### Other Variables:
- NA

### 41) GFF
Setting this variable to TRUE will combine the resulting GFFs of Stringtie, Transdecoder, barrnap and tRNAscan-SE in a single GFF file. Also adding the InterproScan annotation to each gene.
#### Required software:
- seqkit (https://github.com/shenwei356/seqkit)
#### Other Variables:
- NA

### 42) GFF_SORT
Setting this variable to TRUE will sort the resulting combined GFF. 
#### Required software:
- seqkit (https://github.com/shenwei356/seqkit)
#### Other Variables:
- NA
Note: It is slower than gffread and hardcodded to work with this pipeline of analysis, but it preserves the UTR predictions.

## Enriched_GO_Family_Expansion.R
Thuis script takes the annootations obtained with module ORTHOFINDER_EGGNOG_FOR_PHO and the lists of expanded/contracted proteins identified by CAFE5 and estimates enriched GO Terms.
