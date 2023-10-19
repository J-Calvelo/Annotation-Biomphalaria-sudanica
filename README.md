# Annotation-Biomphalaria-sudanica
Bioinformatic pipline utiliced in XXXXXX for the annotation of Biomphalaria sudanica.

The script is subdevided into several blocks of analysis and/or file preparation in order to facilitate re-runs. Which blocks are ran at a time can be controlled by setting to true variables located from lines 12 to 120 (all written in all caps). Other important variables for the run are located betwee lines 123 to 200. Some will be discussed now while others more specific will be mentioned as they become relevant:

### General variables:
1) species_name: Name for the species being annotated
2) genome_file: Full path to the genome assembly in fasta format
3) work_dir: Full path to the working directory (needs to exists
4) threads: Number of computer cores available for the analysis. Used whenever possible.

### 1) CLEAN_RUN
Setting this variable to TRUE will remove all previous runs. Us only to start from a clean state:
#### Required software: 
- NA
#### Other Variables:
- NA
 
### 2) TEA_EARLGREY_HOT
Setting this variable to TRUE will run the earlGrey pipeline for repeat annotation
#### Required software: 
- EarlGrey (https://github.com/google/EarlGrey)
#### Other Variables: 
earlgrey_instalation: Full path to EarlGrey

### 3) TEA_EARLGREY_SUMMARY
Generates a custom summary of EarlGrey results
#### Required software: 
- NA
#### Other Variables:
- NA

# 3) CLEAN_PACBIO_READS
Setting this variable to TRUE will run will retrieve full sequenced isoforms from the PACBIO CSS reads based on the presence of 5' and 3' isoforms
Setting CLEAN_PACBIO_READS_RM_TM to TRUE will also remove intermediary files.
#### Required software: 
- lima (https://github.com/PacificBiosciences/barcoding)
- isoseq3 (https://github.com/PacificBiosciences/IsoSeq)
- bam2fastq (https://github.com/jts/bam2fastq)
- seqkit (https://github.com/shenwei356/seqkit)
#### Other Variables:
- raw_pacbio_reads: Full path to the PACBIO CCS reads file (BAM Format)
- pacbio_primers: Full path to a file with the adapter sequences of the PACBIO CCS readss (FASTA Format)
  
# 4) CLEAN_ILLUMINA_READS 
Setting this variable to TRUE will remove adapter sequences and low quality bases from the paired end Illumina reads using trimmomatic (SLIDINGWINDOW:5:20 MINLEN:25)
#### Required software: 
- trimmomatic (https://github.com/usadellab/Trimmomatic)
- seqkit (https://github.com/shenwei356/seqkit)
#### Other Variables:
- raw_illumina_reads_dir: Full path to the folder with the Illumina read files 
- raw_illumina_1: Identifier for the Left reads, located at the end of the file 
- raw_illumina_2: Identifier for the Left reads, located at the end of the file
Note: Sample names are identified running the comand:
ls $raw_illumina_reads_dir"/"*$raw_illumina_1 | sed "s/.*\///" | sed "s/$raw_illumina_1//"

# 5) PACBIO_MAP
Setting this variable to TRUE will map the PACBIO reads to the genome
#### Required software: 
- minimap2 (https://github.com/lh3/minimap2)
- samtools (https://github.com/samtools/samtools)
#### Other Variables:
- NA

# 6) STAR_MAP
Setting this variable to TRUE will map the Illuimina reads to the genome
#### Required software: 
- STAR (https://github.com/alexdobin/STAR/tree/master)
#### Other Variables:
- NA

# 7) STRINGTIE_TRANS
Setting this variable to TRUE will combine the read mappings of Illumina and PACBIO reads to annotate the genome using Stringtie v2.2
Earlier versions of stringtie will not work
#### Required software: 
- stringtie (https://github.com/gpertea/stringtie)
- gffread (https://github.com/gpertea/gffread)
#### Other Variables:
- NA

Note: This block has the function kill_fai that will remove any index file of the genome that could cause problems with gffread
Command: rm $genome_file".fai"

# 8) FIND_MITOCHONDRIA 
Setting this variable to TRUE will try to identify the mitochondrial genome in the genome assembly based on blast searches against known mitochondrial gene sequences of other biomphalaria. 
Setting FIND_MITOCHONDRIA_RM_TM to TRUE will also remove intermediary files.
#### Required software: 
- BLAST+ (http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs)
#### Other Variables:
- mito_genes: Full path to a fasta file with known mitochondrial genes

Note: Other blocks working with the mitochondria will assume that it's sequence was retrieved and annotated by means not included on this pipeline.

# 9) MAP_MITOCHONDRIA_PACBIO
Setting this variable to TRUE will map PACBIO reads to the mitochondrial genome and count reads asigned to the different features. 
#### Required software: 
- minimap2 (https://github.com/lh3/minimap2)
- samtools (https://github.com/samtools/samtools)
- htseq-count (https://github.com/htseq/htseq)
#### Other Variables:
- mito_genome: Full path to a fasta file with the conting/scaffold with the mitochondrial genome
- mito_gff: Full path to a GFF file with annotated mitochondrial genome

# 10) MAP_MITOCHONDRIA_ILLUMINA
Setting this variable to TRUE will map paired Illumina reads to the mitochondrial genome and count reads asigned to the different features. 
#### Required software: 
- STAR (https://github.com/alexdobin/STAR/tree/master)
- samtools (https://github.com/samtools/samtools)
- htseq-count (https://github.com/htseq/htseq)
#### Other Variables:
- mito_genome: Full path to a fasta file with the conting/scaffold with the mitochondrial genome
- mito_htseq_illumina_gff: Full path to a GFF file with annotated mitochondrial genome

Note: For convininence mito_htseq_illumina_gff has every feature identidied as a "gene" (third column of the gff) to better identify reads that overlap different features.

# 11) MAP_MITOCHONDRIA_ILLUMINA_NEW_ORIGIN
Setting this variable to TRUE will run the same analysis than MAP_MITOCHONDRIA_ILLUMINA but over a genome file with a different origin coordinate
#### Required software: 
- STAR (https://github.com/alexdobin/STAR/tree/master)
- samtools (https://github.com/samtools/samtools)
- htseq-count (https://github.com/htseq/htseq)
#### Other Variables:
new_ori_genome: Full path to a fasta file with the conting/scaffold with the mitochondrial genome
new_ori_htseq_illumina_gff: Full path to a GFF file with annotated mitochondrial genome

Note: This was done for quality controlls

# 12) ORF_PROT_SEQ_1
Setting this variable to TRUE will run half of Transdecoder pipeline for the identification of ORFs in the transcriptome predicted with STRINGTIE_TRANS
#### Required software: 
- Transdecoder (https://github.com/TransDecoder/TransDecoder)
- 


ORF_PROT_SEQ_1=FALSE # runs the necessary predictions and tests




ORF_PROT_SEQ_2=FALSE # Predicts the Proteins

# 11) BUSCO_TESTS
Setting this variable to TRUE will evaluate how complete are the genome assembly and the predicted transcriptome and proteome using BUSCO. Which of the three is to be run can be adjusted with variables: BUSCO_GENOME, BUSCO_RNA and BUSCO_PROT_ALL
#### Required software: 
- Busco (http://busco.ezlab.org/)



# 13) Interprot # Requires the conda environment: conda activate Interprot
RUN_INTERPROT=FALSE

# 14) Orthofinder # requires the conda environment: conda activate Orthofinder
PREPARE_ORTHOFINDER=FALSE
RUN_ORTHOFINDER=FALSE

# 15) GO Terms with Eggnog # requires the conda environment: conda activate Eggnog
RUN_EGGNOG=FALSE

# 16) non coding search # Requires conda environment: conda activate RNA_NonCoding
RUN_TRNA=FALSE
RUN_RIBOSOME=FALSE

# 17) SignalP, TargetP  # Requires conda activate Signal
SEARCH_SIGNALS=FALSE
SEARCH_SECRETOMEP=FALSE # 1st Excludes Signalp and TargetP results. Then runs SecretomeP
SIGNAL_RELOCATION_TABLE=FALSE # SecretomeP is an asshole!! this part makes sense of the output (Re-formats a table and undoes an ID change)

# 18) Run DeepTMHMM
DeepTMHMM=FALSE
TABLE_TMHMM=FALSE

# 19) FREP Search # Serches FREP genes for BLAST and domain structure, and while it checks if there are some weird chimera like gene we should know about. Requires the basic conda environment
SAM_IGSF_SEARCH=FALSE
SEARCH_FREPS=FALSE
FREP_PHYLOGENY=FALSE
SUP_TAB_PAPER_FREP=TRUE

exclude_seq_manual=/home/amanda/PROYECTS/Project_Bsudanica/Bsudanica_Runs/Final_Genome_Annotation/FREP_like_Search/Best_Candidates_ReviewNeded/Manual_check_Missidentifications.txt # File with FREP and CREP candidates to be excluded because of anomalies in the predicted

# 20) Mini_Orthofinder (Orthofinder but just with a couple of species) # requires the conda environment: conda activate Orthofinder
PREPARE_MINI_ORTHOFINDER=FALSE
MINI_ORTHOFINDER=FALSE
MINI_ORTHOFINDER_EGGNOG=FALSE # Requires conca activate Eggnog
MINI_ORTHOFINDER_EGGNOG_FOR_PHO=FALSE
MINI_ORTHOFINDER_INTERPROT=FALSE # Requires conca activate Interprot

# 21) Protein expansion: CAFE  # requires the basic conda environment: conda activate
MINI_CAFE_EXPANSION_PREPARE=FALSE # Re-Pharses Orthofinder's results
MINI_CAFE_EXPANSION_TEST_RUNS=FALSE # Runs CAFE test runs for manual analysis
MINI_CAFE_EXPANSION_DEF_RUN=FALSE # Just runs the definitive run and summariced what families expanded or contracted
MINI_CAFE_EXPANSION_ANOT_INTERPROT=FALSE # Runs interpro for all non-Bsudanica species and store the results # requires conda interprot
MINI_CAFE_EXPANSION_ANOT_EGGNOG=FALSE # Runs Eggnog for all non-Bsudanica species and store the results # requires conda eggnog
MINI_CAFE_EXPANSION_ANOT_PREPARE_GO=FALSE

# 22) Data crossings
DATA_CROSS_INTERPROT_EGGNOG=FALSE # Registra que genes tienen las entradas de interprot identificadas como parte del sistema inmune
DATA_CROSS_INTERPROT_SIGNALP=FALSE
DATA_CROSS_POLYA_REPEATS=FALSE

# 23) Final Summary # Requires the basic conda environment: conda activate
FINAL_SUMMARY=FALSE
FINAL_Generate_GO_Lists=FALSE
FINAL_LOCATION_TABLE=FALSE
FINAL_GFF=FALSE
FINAL_GFF_SORT=FALSE

