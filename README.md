# Annotation-Biomphalaria-sudanica
Bioinformatic pipline utiliced in XXXXXX for the annotation of Biomphalaria sudanica.

The script is subdevided into several blocks of analysis and/or file preparation in order to facilitate re-runs. Which blocks are ran at a time can be controlled by setting to true variables located from lines 12 to 120 (all written in all caps):

### 1) CLEAN_RUN
Setting this variable to TRUE will remove all previous runs. Us only to start from a clean state:
#### Required software: 
- NA

### 2) TEA_EARLGREY_HOT
Setting this variable to TRUE will run the earlGrey pipeline for repeat annotation
Re

# Requires the conda enviroment earlGrey active: conda activate earlGrey. And the following additions to the PATH: PATH=$PATH:/home/amanda/programas/RECON-1.08/bin:/home/amanda/programas/rmblast-2.11.0/bin:/home/amanda/programas/RepeatScout-1.0.6:/home/amanda/programas/RepeatMasker:/home/amanda/programas/RepeatMasker/util:/home/amanda/programas/ucscTwoBitTools:/home/amanda/programas/RepeatModeler-2.0.2a:/home/amanda/programas/EarlGrey
  TEA_EARLGREY_SUMMARY=FALSE # my Summary for the repeats

# 3) Pacbio Cleaning
CLEAN_PACBIO_READS=FALSE # Requires the conda environment TAMA: conda activate TAMA
  CLEAN_PACBIO_READS_RM_TM=FALSE # Remove intermediary files

# 4) Illumina Reads Cleaning # Requires the conda environment Mixed_Transcriptome: conda activate Mixed_Transcriptome
CLEAN_ILLUMINA_READS=FALSE

# 5) PACbio Mappings # Requires the conda environment TAMA: conda activate TAMA
PACBIO_MAP=FALSE

# 6) STAR Mappings # Requires the conda environment Mixed_Transcriptome: conda activate Mixed_Transcriptome
STAR_MAP=FALSE

# 7) Stringtie Annotation # Requires the conda environment Mixed_Transcriptome: conda activate Mixed_Transcriptome
STRINGTIE_TRANS=FALSE

# 8) Busco evaluation # Requires the conda environment Busco: conda activate Busco
BUSCO_TESTS=FALSE
  BUSCO_GENOME=FALSE # Runs the busco for the entire genome
  BUSCO_RNA=FALSE # Runs the busco for the transcriptome
  BUSCO_PROT_ALL=FALSE # Runs the busco for the predicted proteins

# 9) Find Mitochondria # Requires the basic conda environment: conda activate
FIND_MITOCHONDRIA=FALSE
  FIND_MITOCHONDRIA_RM_TM=FALSE # Removes the blast reference

# 10) Map reads to Mitochondria PACBIO # Requires the conda environment TAMA: conda activate TAMA
MAP_MITOCHONDRIA_PACBIO=FALSE
MAP_MITOCHONDRIA_PACBIO_NEW_ORIGIN=FALSE # Confirm that the change in the origin didn't affected the read mapping counts

# 11) Map reads to Mitochondria ILLUMINA #  Requires the basic conda environment: conda activate
MAP_MITOCHONDRIA_ILLUMINA=FALSE
MAP_MITOCHONDRIA_ILLUMINA_NEW_ORIGIN=FALSE # Confirm that the change in the origin didn't affected the read mapping counts

# 12) Get ORFs  #  Requires the basic conda environment: conda activate
ORF_PROT_SEQ_1=FALSE # runs the necessary predictions and tests
ORF_PROT_SEQ_2=FALSE # Predicts the Proteins

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

