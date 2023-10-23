#!/usr/bin/env bash

################################################################################
################################################################################
################################################################################
# Run tasks!!
CLEAN_RUN=FALSE
TEA_EARLGREY_HOT=FALSE
TEA_EARLGREY_SUMMARY=FALSE
CLEAN_PACBIO_READS=FALSE
CLEAN_PACBIO_READS_RM_TM=FALSE
CLEAN_ILLUMINA_READS=FALSE
PACBIO_MAP=FALSE
STAR_MAP=FALSE
STRINGTIE_TRANS=FALSE
FIND_MITOCHONDRIA=FALSE
FIND_MITOCHONDRIA_RM_TM=FALSE
MAP_MITOCHONDRIA_PACBIO=FALSE
MAP_MITOCHONDRIA_PACBIO_NEW_ORIGIN=FALSE
MAP_MITOCHONDRIA_ILLUMINA=FALSE
MAP_MITOCHONDRIA_ILLUMINA_NEW_ORIGIN=FALSE
ORF_PROT_SEQ_1=FALSE
ORF_PROT_SEQ_2=FALSE
BUSCO_TESTS=FALSE
BUSCO_GENOME=FALSE
BUSCO_RNA=FALSE
BUSCO_PROT_ALL=FALSE
RUN_INTERPROT=FALSE
RUN_EGGNOG=FALSE
RUN_TRNA=FALSE
RUN_RIBOSOME=FALSE
SEARCH_SIGNALS=FALSE
SEARCH_SECRETOMEP=FALSE
DeepTMHMM=FALSE
TABLE_TMHMM=FALSE
SAM_IGSF_SEARCH=FALSE
SEARCH_CREP_FREP_GREP=FALSE
CREP_FREP_PHYLOGENY=FALSE
SUP_TAB_PAPER_CREP_FREP_GREP=FALSE
PREPARE_ORTHOFINDER=FALSE
ORTHOFINDER=FALSE
ORTHOFINDER_EGGNOG=FALSE
ORTHOFINDER_EGGNOG_FOR_PHO=FALSE
CAFE_EXPANSION_PREPARE=FALSE
CAFE_EXPANSION_TEST_RUNS=FALSE
CAFE_EXPANSION_DEF_RUN=FALSE
CAFE_EXPANSION_ANOT_EGGNOG=FALSE
DATA_CROSS_INTERPROT_SIGNALP=FALSE
DATA_CROSS_POLYA_REPEATS=FALSE
SUMMARY=FALSE
LOCATION_TABLE=FALSE
GFF=FALSE
GFF_SORT=FALSE

################################################################################
################################################################################
################################################################################

# General variables
species_name=Bsudanica
threads=30
work_dir=/home/amanda/PROYECTS/Project_Bsudanica/Bsudanica_Runs/Final_Genome_Annotation
genome_file=/home/amanda/PROYECTS/Project_Bsudanica/Bsudanica_Runs/Genome_Info/Bsud111_ccs_assembly.fasta

# Input data (data from outside this script like genome or reads)
raw_pacbio_reads=/home/amanda/PROYECTS/Project_Bsudanica/Bsudanica_Runs/PACBIO_Reads/More_Files/ccs_Q20/m64047_210901_223130.ccs.bam
pacbio_primers=/home/amanda/PROYECTS/Project_Bsudanica/Bsudanica_Runs/PACBIO_Reads/Primer.fasta
orthofinder_folder=/home/amanda/PROYECTS/Project_Bsudanica/Bsudanica_Runs/Final_Genome_Annotation/Mini_Orthofinder_Run
other_species_orthofinder=/home/amanda/PROYECTS/Project_Bsudanica/Mini_Orthofinder_cDNA
other_species_mitochondria=/home/amanda/PROYECTS/Project_Bsudanica/Mini_Orthofinder_Mitochondrias/All_mito.fasta
exclude_seq_manual=/home/amanda/PROYECTS/Project_Bsudanica/Bsudanica_Runs/Final_Genome_Annotation/FREP_like_Search/Best_Candidates_ReviewNeded/Manual_check_Missidentifications.txt # File with FREP and CREP candidates to be excluded because of anomalies in the predicted
raw_illumina_reads_dir=/home/amanda/PROYECTS/Project_Bsudanica/Bsudanica_Runs/Illumina_reads
raw_illumina_1=_1.fq.gz
raw_illumina_2=_2.fq.gz
Illumina_adapters=/home/amanda/miniconda3/pkgs/matic-0.39-hdfd78af_2/share/matic/adapters/TruSeq3-PE-2.fa
mito_genes=/home/amanda/PROYECTS/Project_Bsudanica/Bglabrata_Runs/Genome_Info/Mitochondria_Coding.fa
mito_gff=/home/amanda/PROYECTS/Project_Bsudanica/Bsudanica_Runs/Mitocondrial_tests/Final_Anotation_Files/Mitocondrial_htseq.gff # custom made annotation for counting
mito_htseq_illumina_gff=/home/amanda/PROYECTS/Project_Bsudanica/Bsudanica_Runs/Mitocondrial_tests/Final_Anotation_Files/Mitocondrial_htseq_Illumina.gff # custom made annotation for counting
mito_genome=/home/amanda/PROYECTS/Project_Bsudanica/Bsudanica_Runs/Mitocondrial_tests/Mitocondrial_chr.fa
new_ori_genome=/home/amanda/PROYECTS/Project_Bsudanica/Bsudanica_Runs/Mitocondrial_tests/Final_Anotation_Files/Mitocondrial_New_Origin_chr.fasta
new_ori_gff=/home/amanda/PROYECTS/Project_Bsudanica/Bsudanica_Runs/Mitocondrial_tests/Final_Anotation_Files/Mitocondrial_Mito_New_Origin_Htseqcount.gff
new_ori_htseq_illumina_gff=/home/amanda/PROYECTS/Project_Bsudanica/Bsudanica_Runs/Mitocondrial_tests/Final_Anotation_Files/Mitocondrial_Mito_New_Origin_Htseqcount.gff # custom made annotation for counting
uniprot=/home/amanda/PROYECTS/Project_Bsudanica/Uniprot/uniprot_Swiss-Prot.fasta
pfam=/home/amanda/PROYECTS/Project_Bsudanica/PFAM/Pfam-A.hmm
eggnog_data=/home/amanda/PROYECTS/eggnog_22-7-2022
Lu_FREP_prot=/home/amanda/PROYECTS/Project_Bsudanica/References_FREP_CREP_GREP/Lu_etal2020_FREP.fa
Dheilly_FREP_prot=/home/amanda/PROYECTS/Project_Bsudanica/References_FREP_CREP_GREP/Dehilly_etal_2015_FREP.fa
Dheilly_CREP_prot=/home/amanda/PROYECTS/Project_Bsudanica/References_FREP_CREP_GREP/Dehilly_etal_2015_CREP.fa
Dheilly_GREP_prot=/home/amanda/PROYECTS/Project_Bsudanica/References_FREP_CREP_GREP/Dehilly_etal_2015_GREP.fa
igsf_hmmer_profile=/home/amanda/PROYECTS/Project_Bsudanica/Sam_Team_IgSF/Bg_IgSF.hmm

# Software installation
targetp_file=/home/amanda/PROYECTS/Signalp/targetp-2.0/bin
secretomep=/home/amanda/PROYECTS/secretomep-1.0
earlgrey_instalation=/home/amanda/programas/EarlGrey
interprot_path=/home/amanda/PROYECTS/interproscan-5.56-89.0

# Other variables
dash_ids=$(echo "Biomphalaria_glabrata_VectorBaseBB02 Biomphalaria_pfeifferi_SamTeam Biomphalaria_straminea_gigascience")
star_index_num=13
transdecoder_evalue=1e-5
FBD_signatures=$(echo "IPR002181 IPR014716 IPR020837 IPR036056")
Inmunoglobulin_signatures=$(echo "IPR003598 IPR003599 IPR007110 IPR013098 IPR013783 IPR036179")
C_lectin_signatures=$(echo "IPR001304 IPR016186 IPR016187 IPR018378")
Galectin_signatures=$(echo "IPR001079 IPR015533 IPR044156 IPR000922 IPR043159")
EGF_signatures=$(echo "IPR000152 IPR000742 IPR001881 IPR002049 IPR013032 IPR013111 IPR018097 IPR024731 IPR026823 IPR041161")
crep_frep_grep_blast_ident=50
mini_ortho_min_Nspe=3
end_kvalue=5
max_cafe_runs_perK=10
max_number_iterations=10000
ultrametic_tree_root_age=20
Def_Kvalue=1
exclude_mitochondria=contig_7934

kill_fai () {
  if [ -f $genome_file".fai" ]
  then
    rm $genome_file".fai"
  fi
}

run_frep_crep_phylogeny () {
  mafft --thread $threads --quiet --auto $storage_folder"/Run_Analysis.fasta" > $storage_folder"/Run_Analysis.msa"
  trimal -gt 0.8 -st 0 -in $storage_folder"/Run_Analysis.msa" -out $storage_folder"/Run_Analysis.trimmed" -colnumbering >> $storage_folder"/Run_Analysis.cols"
  iqtree --seqtype "AA" -T AUTO -s $storage_folder"/Run_Analysis.trimmed" --prefix $storage_folder"/"$grupo_phy"_Long_Iso_phylogeny" -m MFP -B 1000
}

interprot_basic_mining () {
echo "Prot_ID;Interpro;Int_Description;"$special_header | tr ";" "\t" > $storage_place
for gene in $genes_interes
do
  echo $gene" "$storage_place
  extra_data=$(grep -w -F $gene $extra_data_file | tr "\t" "\n" | sed 1d | sed "s/$/_ooo_/g" | tr -d "\n" | sed "s/_ooo_$/\n/" | sed "s/ /_oo_/g")
  search_gene=$(echo $gene | awk -F "." '{print $1"."$2"."}')

  grep -F $search_gene $work_dir"/Interprot_Results/Results.tsv" | awk -F "\t" '{if ($12!="-") print $12"\t"$13}' | sort -u | awk -F "\t" -v gen=$gene -v extra=$extra_data '{print gen"\t"$1"\t"$2"\t"extra}' | sed "s/_ooo_/\t/g"  | sed "s/_oo_/ /g" >> $storage_place
done
}

cafe_def_run (){
  rm -r $explore_orthofinder"/CAFE_Expansion/DEF_Run"
  mkdir $explore_orthofinder"/CAFE_Expansion/DEF_Run"
  mkdir $explore_orthofinder"/CAFE_Expansion/DEF_Run/Sequences"
  mkdir $explore_orthofinder"/CAFE_Expansion/DEF_Run/Sequences/Per_Family"
  mkdir $explore_orthofinder"/CAFE_Expansion/DEF_Run/Node_Categories"
  mkdir $explore_orthofinder"/CAFE_Expansion/DEF_Run/Temporal"

  cafe5 -o $explore_orthofinder"/CAFE_Expansion/DEF_Run/CAFE_Kvalue_"$Def_Kvalue -k $Def_Kvalue -I $max_number_iterations -c $threads -i $explore_orthofinder"/CAFE_Expansion/Families_CAFE.in" -t $explore_orthofinder"/CAFE_Expansion/SpeciesTree_ultrametric.tre"

  result_model=$(grep "Final Likelihood (-lnL):" $explore_orthofinder"/CAFE_Expansion/DEF_Run/CAFE_Kvalue_"$Def_Kvalue"/"*".txt" | awk -F ":" '{print $1}' | sed "s/_results.txt//" | sed "s/.*\///")
  interest_families=$(grep -w "y" $explore_orthofinder"/CAFE_Expansion/DEF_Run/CAFE_Kvalue_"$Def_Kvalue"/"$result_model"_family_results.txt" | awk -F "\t" '{print $1}')

  species=$(head -n 1 $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" | tr "\t" "\n" | sed 1,3d)
  head -n 1 $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" | tr "\t" "\n" | grep -n .  | sed 1,3d | tr ":" "\t" >> $explore_orthofinder"/CAFE_Expansion/DEF_Run/Temporal/species_cols"
  head -n 1 $explore_orthofinder"/CAFE_Expansion/DEF_Run/CAFE_Kvalue_"$Def_Kvalue"/"$result_model"_change.tab"  | tr "\t" "\n" | sed "s/.*</</" | grep -n . | tr ":" "\t" >> $explore_orthofinder"/CAFE_Expansion/DEF_Run/Temporal/nodes_cols"

  # Generate header based on the given names for the nodes
  header=FamilyID
  count=2
  stop=$(grep -c . $explore_orthofinder"/CAFE_Expansion/DEF_Run/Temporal/nodes_cols")

  while [ $count -le $stop ]
  do
    node=$(sed -n $count"p" $explore_orthofinder"/CAFE_Expansion/DEF_Run/Temporal/nodes_cols" | awk -F "\t" '{print $2}')
    name=$(grep $node $interest_nodes_file | awk -F "\t" '{print $1}')
    header=$(echo $header" "$name)

    count=$(($count + 1))
  done

  echo $header" N_Vanila_Orthogroups Main_Changes" | tr " " "\t" >> $explore_orthofinder"/CAFE_Expansion/DEF_Run/Changing_Prot_Families.tab"

  for fam in $interest_families
  do
    echo $fam
    echo "Checking Orthofinder splitting"
    vanila_orthogroup=$(grep -F -w $fam $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" | awk -F "\t" '{print $2}' | grep -c .)

    echo "Getting Sequences"
    for spe in $species
    do
      extract=$(grep -w $spe $explore_orthofinder"/CAFE_Expansion/DEF_Run/Temporal/species_cols" | awk -F "\t" '{print $1}')
      check_data=$(grep -F -w $fam $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" | tr "\t" "\n" | sed -n $extract"p" | grep -c .)
      if [ $check_data -gt 0 ]
      then
        grep -F -w $fam $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" | tr "\t" "\n" | sed -n $extract"p" | tr -d " " | tr "," "\n" > $explore_orthofinder"/CAFE_Expansion/DEF_Run/Temporal/Extract_genes.tmp"
        seqkit grep -j $threads -f $explore_orthofinder"/CAFE_Expansion/DEF_Run/Temporal/Extract_genes.tmp" $explore_orthofinder"/"$spe".pep" | sed "s/\*$//" >> $explore_orthofinder"/CAFE_Expansion/DEF_Run/Sequences/"$spe"_interest.pep"
        seqkit grep -j $threads -f $explore_orthofinder"/CAFE_Expansion/DEF_Run/Temporal/Extract_genes.tmp" $explore_orthofinder"/"$spe".pep" | sed "s/\*$//" >> $explore_orthofinder"/CAFE_Expansion/DEF_Run/Sequences/Per_Family/"$fam".pep"
      fi
    done

    echo "Classification"

    info_changes=$(grep -w $fam $explore_orthofinder"/CAFE_Expansion/DEF_Run/CAFE_Kvalue_"$Def_Kvalue"/"$result_model"_change.tab")
    significant_nodes=$(grep $fam $explore_orthofinder"/CAFE_Expansion/DEF_Run/CAFE_Kvalue_"$Def_Kvalue"/"$result_model"_asr.tre" | tr "," "\n" | tr ")" "\n" | grep \* | sed "s/.*</</" | sed "s/>.*/>/")
    check_sig_nodes=$(echo $significant_nodes | grep -c .)

    if [ $check_sig_nodes -gt 0 ]
    then
      print_status=
      for nod in $significant_nodes
      do
        not_col=$(grep $nod $explore_orthofinder"/CAFE_Expansion/DEF_Run/Temporal/nodes_cols" | awk -F "\t" '{print $1}')
        node_name=$(grep $nod $interest_nodes_file | awk -F "\t" '{print $1}')
        status=$(grep -w $fam $explore_orthofinder"/CAFE_Expansion/DEF_Run/CAFE_Kvalue_"$Def_Kvalue"/"$result_model"_change.tab" | tr "\t" "\n" | sed -n $not_col"p" | tr -d "+" | awk -F "\t" '{if ($1>0) print "Expansion"; else if ($1<0) print "Reduction"; else print "Stable"}')

        echo $fam >> $explore_orthofinder"/CAFE_Expansion/DEF_Run/Node_Categories/"$node_name"_"$status".txt"
        echo $nod"||"$not_col"||"$node_name"||"$status # example of mistakes <8>||10||Internal_Bg_plus_African_Biom Root||Expansion

        print_status=$(echo $print_status";"$node_name"("$status")" | sed "s/^;//")
      done
    else
      print_status=Unclear
      echo $fam >> $explore_orthofinder"/CAFE_Expansion/DEF_Run/Node_Categories/Unclear.txt"
    fi

    echo $info_changes" "$vanila_orthogroup" "$print_status | tr " " "\t" >> $explore_orthofinder"/CAFE_Expansion/DEF_Run/Changing_Prot_Families.tab"
  done
}

sorting_gff () {
  for chr in $chromosomes
  do
    grep -w -F $chr $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/features.tmp" | sort -n -k 4 >> $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/"$chr".tmp"
    grep -w -F $chr $work_dir"/Final_Anotation_Files/Full_GFF/Intermediary/Concatenated.gff" >> $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/"$chr"_full.tmp"

    count=1
    recorrer=$(grep -c . $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/"$chr".tmp")
    while [ $count -le $recorrer ]
    do
      term_id=$(sed -n $count"p"  $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/"$chr".tmp" | awk -F "\t" '{print $9}' | tr ";" "\n" | grep "ID=" | awk -F "=" '{print $2}' )
      source=$(sed -n $count"p"  $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/"$chr".tmp" | awk -F "\t" '{print $2}')

      echo $chr" "$count"/"$recorrer": "$source

      sed -n $count"p"  $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/"$chr".tmp" >> $work_dir"/Final_Anotation_Files/Full_GFF/Bsudanica_Anotation_Full.gff"

      if [ "$source" == "transdecoder" ]
      then
        echo "We in Transdecoder"
        transcripts=$(grep "Parent="$term_id";" $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/"$chr"_full.tmp" | awk -F "\t" '{if ($3=="mRNA") print $9}' | tr ";" "\n" | grep ID= | awk -F "=" '{print $2}')
        for trans in $transcripts
        do
          grep "ID="$trans";" $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/"$chr"_full.tmp" | awk -F "\t" '{if ($3=="mRNA") print}' >> $work_dir"/Final_Anotation_Files/Full_GFF/Bsudanica_Anotation_Full.gff"
          grep "Parent="$trans$ $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/"$chr"_full.tmp" >> $work_dir"/Final_Anotation_Files/Full_GFF/Bsudanica_Anotation_Full.gff"
        done
      elif [ "$source" == "StringTie" ]
      then
        echo "We in Stringtie"
        transcripts=$(grep "Parent="$term_id$ $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/"$chr"_full.tmp" | awk -F "\t" '{if ($3=="transcript") print $9}' | tr ";" "\n" | grep ID= | awk -F "=" '{print $2}')
        for trans in $transcripts
        do
          grep "ID="$trans";" $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/"$chr"_full.tmp" | awk -F "\t" '{if ($3=="transcript") print}' >> $work_dir"/Final_Anotation_Files/Full_GFF/Bsudanica_Anotation_Full.gff"
          grep "Parent="$trans$ $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/"$chr"_full.tmp" >> $work_dir"/Final_Anotation_Files/Full_GFF/Bsudanica_Anotation_Full.gff"
        done
      elif [ "$source" == "manual" ]
      then
        grep Parent $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/"$chr"_full.tmp" >> $work_dir"/Final_Anotation_Files/Full_GFF/Bsudanica_Anotation_Full.gff"
      elif [ "$source" == "tRNAscan-SE" ]
      then
        grep "Parent="$term_id$ $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/"$chr"_full.tmp" >> $work_dir"/Final_Anotation_Files/Full_GFF/Bsudanica_Anotation_Full.gff"
      fi
      count=$(($count + 1))
    done
  done
}

retrive_domain_locations () {
  seqkit grep -j $threads -p $trans $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" | seqkit seq -i >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/"$Type"_candidates.aa"

  if [ $integrity == "Full" ]
  then
    seqkit grep -j $threads -p $trans $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" | seqkit seq -i >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/"$Type"_Full_candidates.aa"
    pho=$(grep $trans $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" | awk -F "\t" '{print $1}')
    echo $trans": "$pho >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/"$Type"_PHO_candidates.aa"
  fi

  largo=$(seqkit grep -j $threads -p $trans $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/"$Type"_candidates.aa" | seqkit fx2tab -n -l | awk -F "\t" '{print $2}')
  echo "Protein: "$trans" | Integrity: "$integrity >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/"$Type"_candidates_domain_locations.tab"
  echo "Lenght (AA): "$largo >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/"$Type"_candidates_domain_locations.tab"
  echo "Domain Start End (Prob/Evalue)" | tr " " "\t"  >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/"$Type"_candidates_domain_locations.tab"

  if [ "$res_signalp" != "." ]
  then
    grep -w -F $trans $work_dir"/Location_Signals/SignalP_results.txt" | awk -F "\t" '{print $5}' | sed "s/CS pos: //" | sed "s/. Pr://" | awk '{print "SP "$0}' | tr " " "\t" | tr "-" "\t" >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp.file"
  fi

  if [ "$res_igsf1_hit" != "." ]
  then
    grep -F -w $trans $work_dir"/IGsF_Domain_Search/Per_Domain_Table_work.tab" | grep -w -F "IgSF1"| awk -F "\t" '{print $4" "$20" "$21" "$13}' | tr " " "\t" >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp.file"
  fi

  if [ "$res_igsf2_hit" != "." ]
  then
    grep -F -w $trans $work_dir"/IGsF_Domain_Search/Per_Domain_Table_work.tab" | grep -w -F "IgSF2"| awk -F "\t" '{print $4" "$20" "$21" "$13}' | tr " " "\t" >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp.file"
  fi

  if [ $check_Interprot_InmunoGlobin -gt 0 ]
  then
    echo "Description Immunoglobulin like signatures"  >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp2.file"
    echo "Interprot_ID Analysis Signature_ID Signature_Description Interprot_Description" >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp2.file"
    for sign in $Inmunoglobulin_signatures
    do
      check_hit=$(grep -c -w -F $sign $work_dir"/FREP_like_Search/interport.temp")
      if [ $check_hit -gt 0 ]
      then
        grep -w -F $sign $work_dir"/FREP_like_Search/interport.temp" | awk -F "\t" '{print $12"("$13");;"$7";;"$8";;"$9}' | sed "s/ /_ooo_/g" | sed "s/;;/\t/g" >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp.file"
        grep -w -F $sign $work_dir"/FREP_like_Search/interport.temp" | awk -F "\t" '{print $12";;"$4";;"$5";;"$6";;"$13}' | sed "s/;;/\t/g" | sort -u >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp2.file"
      fi
      echo "" >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp2.file"
    done
  fi

  echo $Type" like signatures"  >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp2.file"
  echo "Interprot_ID Analysis Signature_ID Signature_Description Interprot_Description" >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp2.file"
  for sign in $interest_signatures
  do
    check_hit=$(grep -c -w -F $sign $work_dir"/FREP_like_Search/interport.temp")
    if [ $check_hit -gt 0 ]
    then
      grep -w -F $sign $work_dir"/FREP_like_Search/interport.temp" | awk -F "\t" '{print $12"("$13");;"$7";;"$8";;"$9}' | sed "s/ /_ooo_/g" | sed "s/;;/\t/g" >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp.file"
      grep -w -F $sign $work_dir"/FREP_like_Search/interport.temp" | awk -F "\t" '{print $12";;"$4";;"$5";;"$6";;"$13}' | sed "s/;;/\t/g" | sort -u >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp2.file"
    fi
  done
  echo "" >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp2.file"

  echo $Inmunoglobulin_signatures" "$interest_signatures | tr " " "\n" >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp3.file"
  check_more_interprot_hits=$(grep -v -f $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp3.file" $work_dir"/FREP_like_Search/interport.temp" | awk -F "\t" '{print $12}' | grep -c IPR  )
  if [ $check_more_interprot_hits -gt 0 ]
  then
    echo "Other interprot hits"  >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp2.file"
    echo "Interprot_ID Analysis Signature_ID Signature_Description Interprot_Description" >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp2.file"
    moar_interprot=$(grep -v -f $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp3.file" $work_dir"/FREP_like_Search/interport.temp" | awk -F "\t" '{print $12}' | grep IPR | sort -u)
    for sign in $moar_interprot
    do
      grep -w -F $sign $work_dir"/FREP_like_Search/interport.temp" | awk -F "\t" '{print $12"("$13");;"$7";;"$8";;"$9}' | sed "s/ /_ooo_/g" | sed "s/;;/\t/g" >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp.file"
      grep -w -F $sign $work_dir"/FREP_like_Search/interport.temp" | awk -F "\t" '{print $12";;"$4";;"$5";;"$6";;"$13}' | sed "s/;;/\t/g" | sort -u >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp2.file"
    done
    echo "" >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp2.file"
  fi

  sort -n -k 2 $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp.file" | sed "s/_ooo_/ /g" >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/"$Type"_candidates_domain_locations.tab"
  echo ""  >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/"$Type"_candidates_domain_locations.tab"
  cat $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp2.file" >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/"$Type"_candidates_domain_locations.tab"
  echo ""  >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/"$Type"_candidates_domain_locations.tab"

  echo "###################################################################"  >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/"$Type"_candidates_domain_locations.tab"
  echo ""  >> $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/"$Type"_candidates_domain_locations.tab"

  rm $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp.file"
  rm $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp2.file"
  rm $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/temp3.file"
}

overlap_signatures () {
  if [ $check_hit -gt 0 ]
  then
    grep -F -w $sign $work_dir"/FREP_like_Search/interport.temp" | awk -F "\t" '{print $7"\t"$8}' | sort -n -k 1 >> $work_dir"/FREP_like_Search/Overlap_temp/"$trans"_"$sign"_coords"

    start_first_igsf=$(sed -n "1p" $work_dir"/FREP_like_Search/Overlap_temp/"$trans"_"$sign"_coords" | awk -F "\t" '{print $1}')
    end_first_igsf=$(sed -n "1p" $work_dir"/FREP_like_Search/Overlap_temp/"$trans"_"$sign"_coords" | awk -F "\t" '{print $2}')
    count=2
    overlap_feature=1
    while [ $count -le $check_hit ]
    do
      start_test_igsf=$(sed -n $count"p" $work_dir"/FREP_like_Search/Overlap_temp/"$trans"_"$sign"_coords" | awk -F "\t" '{print $1}')
      end_test_igsf=$(sed -n $count"p" $work_dir"/FREP_like_Search/Overlap_temp/"$trans"_"$sign"_coords" | awk -F "\t" '{print $2}')

      if [ $start_test_igsf -le $end_first_igsf -a $start_first_igsf -le $end_test_igsf ]
      then
        end_first_igsf=$end_test_igsf
        count=$(($count + 1))
      else
        start_first_igsf=$start_test_igsf
        end_first_igsf=$end_test_igsf
        overlap_feature=$(($overlap_feature + 1))
        count=$(($count + 1))
      fi
    done
    sign_res=$(echo "X("$overlap_feature")")
  else
    sign_res=$(echo .)
  fi
}

prepare_cafe () {
  rm -r $explore_orthofinder"/CAFE_Expansion"
  mkdir $explore_orthofinder"/CAFE_Expansion"
  mkdir $explore_orthofinder"/CAFE_Expansion/Temporal"

  make_ultrametric.py -r $ultrametic_tree_root_age $results_orthofinder"/Species_Tree/SpeciesTree_rooted.txt"
  mv $results_orthofinder"/Species_Tree/SpeciesTree_rooted.txt.ultrametric.tre" $explore_orthofinder"/CAFE_Expansion/SpeciesTree_ultrametric.tre"

  species=$(head -n 1 $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" | tr "\t" "\n" | sed 1,3d)
  head -n 1 $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" | tr "\t" "\n" | grep -n . | awk -F ":" '{print $2"\t"$1}' | sed 1,3d >> $explore_orthofinder"/CAFE_Expansion/Temporal/Column_Order.temp"
  count=2
  Hog_terms=$(grep -c . $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")

  echo "Desc Hog "$species | tr " " "\t" >> $explore_orthofinder"/CAFE_Expansion/Families_CAFE.in"
  while [ $count -le $Hog_terms ]
  do
    hog=$(sed -n $count"p" $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" | awk -F "\t" '{print $1}')
    registro_conteos=
    spe_with_data=$(sed -n $count"p" $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" | tr "\t" "\n" | sed 1,3d | grep -c .)

    echo $hog" "$count"/"$Hog_terms
    if [ $spe_with_data -ge $mini_ortho_min_Nspe ]
    then
      for spe in $species
      do
        spe_num=$(grep $spe $explore_orthofinder"/CAFE_Expansion/Temporal/Column_Order.temp" | awk -F "\t" '{print $2}')
        N_genes=$(sed -n $count"p" $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" | tr "\t" "\n" | sed -n $spe_num"p" | tr -d " " | tr "," "\n" | grep -c .)
        registro_conteos=$(echo $registro_conteos" "$N_genes | sed "s/^ //")
      done
      echo "(null) "$hog" "$registro_conteos | tr " " "\t" >> $explore_orthofinder"/CAFE_Expansion/Families_CAFE.in"
    fi
    count=$(($count + 1))
  done
}

run_cafe () {
  rm -r $explore_orthofinder"/CAFE_Expansion/Runs"
  mkdir $explore_orthofinder"/CAFE_Expansion/Runs"

  echo "Kvalue Run Likelyhood Lambda" | tr " " "\t" > $explore_orthofinder"/CAFE_Expansion/Runs/All_runs_results.txt"

  Kvalue=1
  while [ $Kvalue -le $end_kvalue ]
  do
    run=1
    while [ $run -le $max_cafe_runs_perK ]
    do
      echo "Run: "$run" Kvalue: "$Kvalue
      cafe5 -o $explore_orthofinder"/CAFE_Expansion/Runs/CAFE_Kvalue"$Kvalue"_run"$run -k $Kvalue -I $max_number_iterations -c $threads -i $explore_orthofinder"/CAFE_Expansion/Families_CAFE.in" -t $explore_orthofinder"/CAFE_Expansion/SpeciesTree_ultrametric.tre"

      result_file=$(grep "Final Likelihood (-lnL):" $explore_orthofinder"/CAFE_Expansion/Runs/CAFE_Kvalue"$Kvalue"_run"$run"/"*"txt" | awk -F ":" '{print $1}')
      check_failure=$(grep -w -c failures $result_file)

      if [ $check_failure -eq 0 ]
      then
        like=$(grep "Final Likelihood (-lnL):" $result_file | awk -F ":" '{print $2}' | tr -d " ")
        lambda=$(grep "Lambda:" $result_file | awk -F ":" '{print $2}' | tr -d " ")
        echo $Kvalue" "$run" "$like" "$lambda | tr " " "\t" >> $explore_orthofinder"/CAFE_Expansion/Runs/All_runs_results.txt"
        run=$(($run + 1))
      else
        # rm -r $explore_orthofinder"/CAFE_Expansion/Runs/CAFE_Kvalue"$Kvalue"_run"*
        Kvalue=$(($end_kvalue + $end_kvalue))
        run=$(($max_cafe_runs_perK  + $max_cafe_runs_perK))
      fi
    done
    Kvalue=$(($Kvalue + 1))
  done
}

final_gff_interpro () {
  check_interprot=$(grep -F $gene_id $work_dir"/Interprot_Results/Results.tsv" | grep -c "IPR")
  if [ $check_interprot -gt 0 ]
  then
    interprot_anotation=$(grep -F $gene_id $work_dir"/Interprot_Results/Results.tsv" | awk -F "\t" '{print $12" ("$13")"}' | grep "IPR" | sort -u | tr "\n" "%" | sed "s/%$/\n/" | sed "s/%/---/g" | awk '{print "Interpro="$0}' | sed "s/ /_ooo_/g")
  else
    interprot_anotation=$(echo "Interpro=NA")
  fi
}

if [ "$CLEAN_RUN" == "TRUE" ]
then
  echo "You asked for it. Removing previous data"
  rm -r $work_dir
  mkdir $work_dir
fi

if [ "$TEA_EARLGREY_HOT" == "TRUE" ]
then
  # conda activate earlGrey

  echo "Live long. And Prosper" # earlGrey is a Star Trek reference
  echo "Do not stare at the cow"
  rm -r $work_dir"/"$species_name"EarlGrey"
  $earlgrey_instalation"/earlGrey" -g $genome_file -s $species_name -o $work_dir -r eumetazoa -t $threads
fi

if [ "$TEA_EARLGREY_SUMMARY" == "TRUE" ]
then
  families=$(echo "DNA LINE LTR RC Retroposon Satellite SINE srpRNA Unknown")
  rm -r $work_dir"/"$species_name"EarlGrey/Manual_Sumary"
  mkdir $work_dir"/"$species_name"EarlGrey/Manual_Sumary"

  count=2
  stop=$(grep -c . $work_dir"/"$species_name"EarlGrey/"$species_name"_summaryFiles/"$species_name".filteredRepeats.gff")
  current_chr=$(sed -n 1p $work_dir"/"$species_name"EarlGrey/"$species_name"_summaryFiles/"$species_name".filteredRepeats.gff" | awk -F "\t" '{print $1}')
  current_start=$(sed -n 1p $work_dir"/"$species_name"EarlGrey/"$species_name"_summaryFiles/"$species_name".filteredRepeats.gff" | awk -F "\t" '{print $4}')
  current_end=$(sed -n 1p $work_dir"/"$species_name"EarlGrey/"$species_name"_summaryFiles/"$species_name".filteredRepeats.gff" | awk -F "\t" '{print $5}')

  while [ $count -le $stop ]
  do
    echo $count" of "$stop

    next_chr=$(sed -n $count"p" $work_dir"/"$species_name"EarlGrey/"$species_name"_summaryFiles/"$species_name".filteredRepeats.gff" | awk -F "\t" '{print $1}')
    next_start=$(sed -n $count"p" $work_dir"/"$species_name"EarlGrey/"$species_name"_summaryFiles/"$species_name".filteredRepeats.gff" | awk -F "\t" '{print $4}')
    next_end=$(sed -n $count"p" $work_dir"/"$species_name"EarlGrey/"$species_name"_summaryFiles/"$species_name".filteredRepeats.gff" | awk -F "\t" '{print $5}')

    if [ "$next_chr" == "$current_chr" ]
    then
      if [ $next_start -le $current_end -a $current_start -le $next_end ]
      then
        definitive_start=$current_start
        definitive_end=$next_end
        count=$(($count + 1))

        test_next=TRUE
        while [ "$test_next" == "TRUE" ]
        do
          next_chr=$(sed -n $count"p" $work_dir"/"$species_name"EarlGrey/"$species_name"_summaryFiles/"$species_name".filteredRepeats.gff" | awk -F "\t" '{print $1}')
          next_start=$(sed -n $count"p" $work_dir"/"$species_name"EarlGrey/"$species_name"_summaryFiles/"$species_name".filteredRepeats.gff" | awk -F "\t" '{print $4}')
          next_end=$(sed -n $count"p" $work_dir"/"$species_name"EarlGrey/"$species_name"_summaryFiles/"$species_name".filteredRepeats.gff" | awk -F "\t" '{print $5}')

          if [ "$next_chr" == "$current_chr" ]
          then
            echo "$next_start -le $definitive_end"
            if [ $next_start -le $definitive_end ]
            then
              definitive_end=$next_end
              count=$(($count + 1))
            else
              test_next=FALSE
            fi
          else
            test_next=FALSE
          fi
          echo $test_next
        done
        prev=$(($count - 1))
        echo "Combine_entry_"$prev" "$current_chr" "$definitive_start" "$definitive_end | tr " " "\t" >> $work_dir"/"$species_name"EarlGrey/Manual_Sumary/Temp.tab"
      else
        prev=$(($count - 1))
        sed -n $prev"p" $work_dir"/"$species_name"EarlGrey/"$species_name"_summaryFiles/"$species_name".filteredRepeats.gff" | awk -F "\t" -v entry=$prev '{print "Entry_"entry"\t"$1"\t"$4"\t"$5}' >> $work_dir"/"$species_name"EarlGrey/Manual_Sumary/Temp.tab"
        current_chr=$next_chr
        current_start=$next_start
        current_end=$next_end
      fi
    else
      prev=$(($count - 1))
      sed -n $prev"p" $work_dir"/"$species_name"EarlGrey/"$species_name"_summaryFiles/"$species_name".filteredRepeats.gff" | awk -F "\t" -v entry=$prev '{print "Entry_"entry"\t"$1"\t"$4"\t"$5}' >> $work_dir"/"$species_name"EarlGrey/Manual_Sumary/Temp.tab"
      current_chr=$next_chr
      current_start=$next_start
      current_end=$next_end
    fi
    count=$(($count + 1))
  done
  repeat_element_len_bases=$(awk -F "\t" '{print $4-$3+1}' $work_dir"/"$species_name"EarlGrey/Manual_Sumary/Temp.tab" | awk -F "\t" '{ sum += $2 } END {print sum}')
  genome_len_bases=$(seqkit fx2tab -n -l $genome_file | awk -F "\t" '{ sum += $2 } END {print sum}')
  percentage=$(echo $repeat_element_len_bases" "$genome_len_bases | awk -F " " '{printf "%.2f\n", $1/$2 }')

  echo "Total Repeated Bases: "$repeat_element_len_bases"/"$genome_len_bases" ("$percentage"%)" >> $work_dir"/"$species_name"EarlGrey/Manual_Sumary/My_summary.tab"
  echo "" >> $work_dir"/"$species_name"EarlGrey/Manual_Sumary/My_summary.tab"
  echo "" >> $work_dir"/"$species_name"EarlGrey/Manual_Sumary/My_summary.tab"

  total_bases_repeats=0
  for fam in $families
  do
    general_lenght=$(awk -F "\t" '{print $3"\t"$5-$4+1}' $work_dir"/"$species_name"EarlGrey/"$species_name"_summaryFiles/"$species_name".filteredRepeats.gff" | grep -w $fam | awk -F "\t" '{ sum += $2 } END {print sum}')
    total_bases_repeats=$(($total_bases_repeats + $general_lenght))

    general_total=$(awk -F "\t" '{print $3}' $work_dir"/"$species_name"EarlGrey/"$species_name"_summaryFiles/"$species_name".filteredRepeats.gff" | grep -w -c $fam)
    sub_class=$(awk -F "\t" '{print $3}' $work_dir"/"$species_name"EarlGrey/"$species_name"_summaryFiles/"$species_name".filteredRepeats.gff" | grep -w $fam | sort -u)

    echo $fam" General "$general_total" "$general_lenght
    echo $fam" General "$general_total" "$general_lenght | tr " " "\t" >> $work_dir"/"$species_name"EarlGrey/Manual_Sumary/My_summary.tab"

    for sub in $sub_class
    do
      echo $fam" "$sub

      sub_num=$(grep -P -c $sub"\t" $work_dir"/"$species_name"EarlGrey/"$species_name"_summaryFiles/"$species_name".filteredRepeats.gff")
      sub_len=$(awk -F "\t" '{print $3"\t"$5-$4+1}' $work_dir"/"$species_name"EarlGrey/"$species_name"_summaryFiles/"$species_name".filteredRepeats.gff" | grep -P $sub"\t" | awk -F "\t" '{ sum += $2 } END {print sum}')
      echo $fam" "$sub" "$sub_num" "$sub_len | tr " " "\t" >> $work_dir"/"$species_name"EarlGrey/Manual_Sumary/My_summary.tab"
    done

    echo "" >> $work_dir"/"$species_name"EarlGrey/Manual_Sumary/My_summary.tab"
    echo "" >> $work_dir"/"$species_name"EarlGrey/Manual_Sumary/My_summary.tab"
  done
fi

if [ "$CLEAN_PACBIO_READS" == "TRUE" ]
then
  # conda activate TAMA
  rm -r $work_dir"/PACbio_Clean_Reads"
  mkdir $work_dir"/PACbio_Clean_Reads"
  mkdir $work_dir"/PACbio_Clean_Reads/Temp"

  echo "Run Lima"
  lima --ccs --dump-clips $raw_pacbio_reads $pacbio_primers $work_dir"/PACbio_Clean_Reads/Temp/Reads_Pacbio_lima.bam"
  echo "Refine"
  isoseq3 refine --require-polya $work_dir"/PACbio_Clean_Reads/Temp/Reads_Pacbio_lima.bam" $pacbio_primers $work_dir"/PACbio_Clean_Reads/Temp/Reads_Pacbio_full_len_trans.bam"
  echo "make as fastq"
  bam2fastq -o $work_dir"/PACbio_Clean_Reads/Reads_work" $work_dir"/PACbio_Clean_Reads/Temp/Reads_Pacbio_full_len_trans.bam"

  # Confirmar el archivo de salida
  seqkit stats -a -T $work_dir"/PACbio_Clean_Reads/Reads_work.fastq.gz" >> $work_dir"/PACbio_Clean_Reads/Basic_metrics.tab"

  if [ "$CLEAN_PACBIO_READS_RM_TM" == "TURE" ]
  then
    rm -r $work_dir"/PACbio_Clean_Reads/Temp"
  fi
fi

if [ "$CLEAN_ILLUMINA_READS" == "TRUE" ]
then
  # conda activate Mixed_Transcriptome
  rm -r $work_dir"/Illumina_Clean_Reads"
  mkdir $work_dir"/Illumina_Clean_Reads"

  echo "Trimming Reads"
  read_files_illumina=$(ls $raw_illumina_reads_dir"/"*$raw_illumina_1 | sed "s/.*\///" | sed "s/$raw_illumina_1//")

  for read in $read_files_illumina
  do
    echo "Procesando: "$read
    trimmomatic PE $raw_illumina_reads_dir"/"$read$raw_illumina_1 $raw_illumina_reads_dir"/"$read$raw_illumina_2 -baseout $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read.gz' -threads $threads ILLUMINACLIP:$Illumina_adapters":2:30:10" SLIDINGWINDOW:5:20 MINLEN:25
  done
  seqkit stats -a -T $work_dir"/Illumina_Clean_Reads/"*".gz" >> $work_dir"/Illumina_Clean_Reads/Basic_metrics.tab"
fi

if [ "$PACBIO_MAP" == "TRUE" ]
then
  # conda activate TAMA
  # Note: Wait for better fine toonment of the read mappings
  rm -r $work_dir"/PACBIO_MAP"
  mkdir $work_dir"/PACBIO_MAP"

  echo "Minimap2"
  minimap2 --secondary=no -t $threads -ax splice:hq -uf -a $genome_file $work_dir"/PACbio_Clean_Reads/Reads_work.fastq.gz" >> $work_dir"/PACBIO_MAP/Minimap2_Pacbio_aln.sam"
  echo "Bamification and Sorting"
  samtools sort -@ $threads -O BAM -o $work_dir"/PACBIO_MAP/Minimap2_Pacbio_aln.bam" $work_dir"/PACBIO_MAP/Minimap2_Pacbio_aln.sam"
  rm $work_dir"/PACBIO_MAP/Minimap2_Pacbio_aln.sam"
fi

if [ "$STAR_MAP" == "TRUE" ]
then
  # conda activate Mixed_Transcriptome
  read_files_illumina=$(ls $raw_illumina_reads_dir"/"*$raw_illumina_1 | sed "s/.*\///" | sed "s/$raw_illumina_1//")
  echo $read_files_illumina
  rm -r $work_dir"/STAR_MAP"
  mkdir $work_dir"/STAR_MAP"
  mkdir $work_dir"/STAR_MAP/Reference"

  echo "STAR Reference"
  STAR --runMode genomeGenerate --runThreadN $threads --genomeDir $work_dir"/STAR_MAP/Reference" --genomeFastaFiles $genome_file --genomeChrBinNbits $star_index_num --genomeSAindexNbases $star_index_num
  for read in $read_files_illumina
  do
    echo $read
    echo "Unzip 1P"
    gzip -d -k $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_1P.gz'
    echo "Unzip 2P"
    gzip -d -k $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_2P.gz'

    echo "Mapping"
    STAR --runMode alignReads --runThreadN $threads --genomeSAindexNbases $star_index_num --outFileNamePrefix $work_dir"/STAR_MAP/"$read"_" --outSAMtype BAM SortedByCoordinate --genomeDir $work_dir"/STAR_MAP/Reference" --readFilesIn $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_1P' $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_2P'

    rm $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_1P'
    rm $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_2P'
  done
fi

if [ "$STRINGTIE_TRANS" == "TRUE" ]
then
  # conda activate Mixed_Transcriptome
  rm -r $work_dir"/Mixed_Anot_Stringtie"
  rm $work_dir"/STAR_MAP/All_work_illumina.bam"

  mkdir $work_dir"/Mixed_Anot_Stringtie"
  mkdir $work_dir"/Mixed_Anot_Stringtie/FASTA"

  echo "Merging illumina files..."
  samtools merge -@ $threads $work_dir"/STAR_MAP/All_merged.bam" $work_dir"/STAR_MAP/"*"_Aligned.sortedByCoord.out.bam"
  samtools sort -O BAM -@ $threads -o $work_dir"/STAR_MAP/All_work_illumina.bam" $work_dir"/STAR_MAP/All_merged.bam"
  rm $work_dir"/STAR_MAP/All_merged.bam"

  echo "Stringtie time"
  stringtie --mix --conservative -l "BSUD" -v -p $threads -o $work_dir"/Mixed_Anot_Stringtie/PolyA_Transcripts_anot.gff" $work_dir"/STAR_MAP/All_work_illumina.bam" $work_dir"/PACBIO_MAP/Minimap2_Pacbio_aln.bam"

  kill_fai
  gffread -W -w $work_dir"/Mixed_Anot_Stringtie/FASTA/PolyA_Transcripts.fasta" -F -g $genome_file $work_dir"/Mixed_Anot_Stringtie/PolyA_Transcripts_anot.gff"
  seqkit stats -a -T $work_dir"/Mixed_Anot_Stringtie/FASTA/PolyA_Transcripts.fasta" >> $work_dir"/Mixed_Anot_Stringtie/FASTA/PolyA_Transcripts.stats"
fi

# Find Mitochondrial Genes
if [ "$FIND_MITOCHONDRIA" == "TRUE" ]
then
  # conda activate
  rm -r $work_dir"/Mitochondria"
  mkdir $work_dir"/Mitochondria"
  mkdir $work_dir"/Mitochondria/Find_it"
  mkdir $work_dir"/Mitochondria/Find_it/Reference"

  echo "Finding the mitochondria!... again"
  makeblastdb -in $genome_file -dbtype nucl -parse_seqids -input_type fasta -out $work_dir"/Mitochondria/Find_it/Reference/Genome_ref"
  blastn -query $mito_genes -db $work_dir"/Mitochondria/Find_it/Reference/Genome_ref" -outfmt "6 std" -out $work_dir"/Mitochondria/Find_it/Mitochondria_Blast.raw" -perc_identity 80 -num_threads $threads

  if [ "$FIND_MITOCHONDRIA_RM_TM" == "TRUE" ]
  then
    rm -r $work_dir"/Mitochondria/Find_it/Reference"
  fi
fi

if [ "$MAP_MITOCHONDRIA_PACBIO" == "TRUE" ]
then
  # conda activate TAMA
  rm -r $work_dir"/Mitochondria/PACbio_Map"
  rm -r $work_dir"/Mitochondria/Htseq_Counts"

  mkdir $work_dir"/Mitochondria/PACbio_Map"
  mkdir $work_dir"/Mitochondria/Htseq_Counts"

  minimap2 --secondary=no -t $threads -ax splice:hq -uf -a $mito_genome $work_dir"/PACbio_Clean_Reads/Reads_work.fastq.gz" >> $work_dir"/Mitochondria/PACbio_Map/mapped_raw_reads.sam"
  samtools view -@ $threads -F 4 -O BAM -o $work_dir"/Mitochondria/PACbio_Map/mapped_unsort_reads.bam" $work_dir"/Mitochondria/PACbio_Map/mapped_raw_reads.sam"
  samtools sort -@ $threads -O BAM -o $work_dir"/Mitochondria/PACbio_Map/Mito_pacbio_reads.bam" $work_dir"/Mitochondria/PACbio_Map/mapped_unsort_reads.bam"
  rm $work_dir"/Mitochondria/PACbio_Map/mapped_raw_reads.sam"
  rm $work_dir"/Mitochondria/PACbio_Map/mapped_unsort_reads.bam"

  samtools index $work_dir"/Mitochondria/PACbio_Map/Mito_pacbio_reads.bam"

  htseq-count -r pos --nonunique=all --idattr=Name --format=bam --type=gene --samout $work_dir"/Mitochondria/Htseq_Counts/Mito_pacbio_reads_HseqCount.sam" -q $work_dir"/Mitochondria/PACbio_Map/Mito_pacbio_reads.bam" $mito_gff >> $work_dir"/Mitochondria/Htseq_Counts/Pacbio_reads.count"
  sed "s/.*XF:Z://" $work_dir"/Mitochondria/Htseq_Counts/Mito_pacbio_reads_HseqCount.sam" | sort | uniq -c | sed "s/\t+/+/g" | sed "s/\t]/]/" >> $work_dir"/Mitochondria/Htseq_Counts/PACbio_reads_overlap.count"
fi

################################################################################
################################################################################
################################################################################

if [ "$MAP_MITOCHONDRIA_PACBIO_NEW_ORIGIN" == "TRUE" ]
then
  # conda activate TAMA
  rm -r $work_dir"/Mitochondria/New_Origin"
  mkdir $work_dir"/Mitochondria/New_Origin"
  mkdir $work_dir"/Mitochondria/New_Origin/PACbio_Map"
  mkdir $work_dir"/Mitochondria/New_Origin/Htseq_Counts"

  minimap2 --secondary=no -t $threads -ax splice:hq -uf -a $new_ori_genome $work_dir"/PACbio_Clean_Reads/Reads_work.fastq.gz" >> $work_dir"/Mitochondria/New_Origin/PACbio_Map/mapped_raw_reads.sam"
  samtools view -@ $threads -F 4 -O BAM -o $work_dir"/Mitochondria/New_Origin/PACbio_Map/mapped_unsort_reads.bam" $work_dir"/Mitochondria/New_Origin/PACbio_Map/mapped_raw_reads.sam"
  samtools sort -@ $threads -O BAM -o $work_dir"/Mitochondria/New_Origin/PACbio_Map/Mito_pacbio_reads.bam" $work_dir"/Mitochondria/New_Origin/PACbio_Map/mapped_unsort_reads.bam"
  rm $work_dir"/Mitochondria/New_Origin/PACbio_Map/mapped_raw_reads.sam"
  rm $work_dir"/Mitochondria/New_Origin/PACbio_Map/mapped_unsort_reads.bam"

  samtools index $work_dir"/Mitochondria/New_Origin/PACbio_Map/Mito_pacbio_reads.bam"

  htseq-count -r pos --nonunique=all --idattr=Name --format=bam --type=gene --samout $work_dir"/Mitochondria/New_Origin/Htseq_Counts/Mito_pacbio_reads_HseqCount.sam" -q $work_dir"/Mitochondria/New_Origin/PACbio_Map/Mito_pacbio_reads.bam" $new_ori_gff >> $work_dir"/Mitochondria/New_Origin/Htseq_Counts/Pacbio_reads.count"
  sed "s/.*XF:Z://" $work_dir"/Mitochondria/New_Origin/Htseq_Counts/Mito_pacbio_reads_HseqCount.sam" | sort | uniq -c | sed "s/\t+/+/g" | sed "s/\t]/]/" >> $work_dir"/Mitochondria/New_Origin/Htseq_Counts/PACbio_reads_overlap.count"
fi

################################################################################
################################################################################
################################################################################

if [ "$MAP_MITOCHONDRIA_ILLUMINA" == "TRUE" ]
then
  # conda activate
  read_files_illumina=$(ls $raw_illumina_reads_dir"/"*$raw_illumina_1 | sed "s/.*\///" | sed "s/$raw_illumina_1//")

  rm -r $work_dir"/Mitochondria/Illumina_Map"
  mkdir $work_dir"/Mitochondria/Illumina_Map"
  mkdir $work_dir"/Mitochondria/Illumina_Map/Reference"
  echo ""
  echo ""

  for read in $read_files_illumina
  do
    echo  "Align: "$read

    STAR --runMode genomeGenerate --runThreadN $threads --genomeDir $work_dir"/Mitochondria/Illumina_Map_Test/Reference" --genomeFastaFiles $mito_genome --genomeChrBinNbits 5 --genomeSAindexNbases 5

    for read in $read_files_illumina
    do
      echo  "Align: "$read
      echo "Sam managing"

      echo $read
      echo "Unzip 1P"
      gzip -d -k $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_1P.gz'
      echo "Unzip 2P"
      gzip -d -k $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_2P.gz'

      echo "Mapping"
      STAR --runMode alignReads --runThreadN $threads --genomeSAindexNbases 5 --limitBAMsortRAM 1025931299 --outFileNamePrefix $work_dir"/Mitochondria/Illumina_Map/"$read"_illumina_" --outSAMtype BAM SortedByCoordinate --genomeDir $work_dir"/Mitochondria/Illumina_Map_Test/Reference" --readFilesIn $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_1P' $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_2P'

      samtools index $work_dir"/Mitochondria/Illumina_Map/"$read"_illumina_Aligned.sortedByCoord.out.bam"
      htseq-count -r pos --idattr=Name --format=bam --type=gene -q $work_dir"/Mitochondria/Illumina_Map/"$read"_illumina_Aligned.sortedByCoord.out.bam" $mito_htseq_illumina_gff > $work_dir"/Mitochondria/Htseq_Counts/"$read"_Illumina_reads.count"

      rm $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_1P'
      rm $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_2P'
    done
  done
fi

################################################################################
################################################################################
################################################################################

if [ "$MAP_MITOCHONDRIA_ILLUMINA_NEW_ORIGIN" == "TRUE" ]
then
  # conda activate
  read_files_illumina=$(ls $raw_illumina_reads_dir"/"*$raw_illumina_1 | sed "s/.*\///" | sed "s/$raw_illumina_1//")

  rm -r $work_dir"/Mitochondria/New_Origin/Illumina_Map_Test"
  mkdir $work_dir"/Mitochondria/New_Origin/Illumina_Map_Test"
  mkdir $work_dir"/Mitochondria/New_Origin/Illumina_Map_Test/Reference"

  STAR --runMode genomeGenerate --runThreadN $threads --genomeDir $work_dir"/Mitochondria/New_Origin/Illumina_Map_Test/Reference" --genomeFastaFiles $new_ori_genome --genomeChrBinNbits 5 --genomeSAindexNbases 5
  echo ""
  echo ""

  for read in $read_files_illumina
  do
    echo  "Align: "$read
    echo "Sam managing"

    echo $read
    echo "Unzip 1P"
    gzip -d -k $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_1P.gz'
    echo "Unzip 2P"
    gzip -d -k $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_2P.gz'

    echo "Mapping"
    STAR --runMode alignReads --runThreadN $threads --genomeSAindexNbases 5 --limitBAMsortRAM 1025931299 --outFileNamePrefix $work_dir"/Mitochondria/New_Origin/Illumina_Map/"$read"_illumina_" --outSAMtype BAM SortedByCoordinate --genomeDir $work_dir"/Mitochondria/New_Origin/Illumina_Map_Test/Reference" --readFilesIn $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_1P' $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_2P'

    samtools index $work_dir"/Mitochondria/New_Origin/Illumina_Map/"$read"_illumina_Aligned.sortedByCoord.out.bam"
    htseq-count -r pos --idattr=Name --format=bam --type=gene -q $work_dir"/Mitochondria/New_Origin/Illumina_Map/"$read"_illumina_Aligned.sortedByCoord.out.bam" $new_ori_htseq_illumina_gff > $work_dir"/Mitochondria/New_Origin/Htseq_Counts/"$read"_Illumina_reads.count"

    rm $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_1P'
    rm $work_dir"/Illumina_Clean_Reads/"$read'_trimmed_read_2P'
  done
fi

################################################################################
################################################################################
################################################################################

if [ "$ORF_PROT_SEQ_1" == "TRUE" ]
then
  rm -r $work_dir"/TransDecoder"
  rm -r $work_dir"/TransDecoder/"*"__checkpoints"*
  mkdir $work_dir"/TransDecoder"
  mkdir $work_dir"/TransDecoder/WorkFiles"
  mkdir $work_dir"/TransDecoder/LongOrfs"
  mkdir $work_dir"/TransDecoder/BLASTP_Ref"

  # Exclude Mitochondria
  exclude_things=$(echo "loc:"$exclude_mitochondria"\|" )
  seqkit grep -j $threads -v -n -r -p $exclude_things $work_dir"/Mixed_Anot_Stringtie/FASTA/PolyA_Transcripts.fasta" >> $work_dir"/TransDecoder/WorkFiles/PolyA_Transcripts_work.fasta"
  grep -v -w $exclude_mitochondria $work_dir"/Mixed_Anot_Stringtie/PolyA_Transcripts_anot.gff" >> $work_dir"/TransDecoder/WorkFiles/Stringtie_Transcripts_work.gff"

  makeblastdb -in $uniprot -dbtype prot -input_type fasta -out $work_dir"/TransDecoder/BLASTP_Ref/Uniprot_ref"

  TransDecoder.LongOrfs -m 50 -S -t $work_dir"/TransDecoder/WorkFiles/PolyA_Transcripts_work.fasta" -O $work_dir"/TransDecoder/LongOrfs"
  blastp -query $work_dir"/TransDecoder/LongOrfs/longest_orfs.pep" -db $work_dir"/TransDecoder/BLASTP_Ref/Uniprot_ref" -max_target_seqs 1 -outfmt 6 -evalue $transdecoder_evalue -num_threads $threads > $work_dir"/TransDecoder/WorkFiles/blastp.outfmt6"
  hmmscan --cpu $threads --domtblout $work_dir"/TransDecoder/WorkFiles/pfam.domtblout" $pfam $work_dir"/TransDecoder/LongOrfs/longest_orfs.pep"
fi

################################################################################
################################################################################
################################################################################

if [ "$ORF_PROT_SEQ_2" == "TRUE" ]
then
  rm -r $work_dir"/TransDecoder/Predict"
  rm -r $work_dir"/TransDecoder/"*"__checkpoints"*

  rm $work_dir"/TransDecoder/Bsudanica_annotation_polyA.gff3"
  rm $work_dir"/TransDecoder/Bsudanica_annotation_mapped_prot.fasta"
  mkdir $work_dir"/TransDecoder/Predict"

  TransDecoder.Predict --single_best_only -t $work_dir"/TransDecoder/WorkFiles/PolyA_Transcripts_work.fasta" --retain_pfam_hits $work_dir"/TransDecoder/WorkFiles/pfam.domtblout" --retain_blastp_hits $work_dir"/TransDecoder/WorkFiles/blastp.outfmt6" -O $work_dir"/TransDecoder/LongOrfs"
  rm pipeliner.*
  mv PolyA_Transcripts_work* $work_dir"/TransDecoder/Predict"

  gtf_to_alignment_gff3.pl $work_dir"/TransDecoder/WorkFiles/Stringtie_Transcripts_work.gff" >> $work_dir"/TransDecoder/transcripts.gff3"
  cdna_alignment_orf_to_genome_orf.pl $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.gff3" $work_dir"/TransDecoder/transcripts.gff3" $work_dir"/TransDecoder/WorkFiles/PolyA_Transcripts_work.fasta" >> $work_dir"/TransDecoder/Bsudanica_annotation_polyA.gff3"
  rm $work_dir"/TransDecoder/transcripts.gff3"

  kill_fai
  gffread -W -C -w $work_dir"/TransDecoder/Bsudanica_annotation_mapped_prot.fasta" -g $genome_file $work_dir"/TransDecoder/Bsudanica_annotation_polyA.gff3"
fi

################################################################################
################################################################################
################################################################################

if [ "$BUSCO_TESTS" == "TRUE" ]
then
  # conda activate Busco
  rm -r $work_dir"/Mixed_Anot_Stringtie/BUSCO"
  mkdir $work_dir"/Mixed_Anot_Stringtie/BUSCO"

  # Genome
  if [ "$BUSCO_GENOME" == "TRUE" ]
  then
    busco -i $genome_file -o "Genome_Busco_molusca" -m genome -l mollusca
    mv "Genome_Busco_molusca" $work_dir"/Mixed_Anot_Stringtie/BUSCO"
  fi

  # RNA
  if [ "$BUSCO_RNA" == "TRUE" ]
  then
    busco -i $work_dir"/Mixed_Anot_Stringtie/FASTA/PolyA_Transcripts.fasta" -o "PolyA_Transcripts_busco_molusca" -m tran -l mollusca
    mv "PolyA_Transcripts_busco_molusca" $work_dir"/Mixed_Anot_Stringtie/BUSCO"
  fi

  # Protein
  if [ "$BUSCO_PROT_ALL" == "TRUE" ]
  then
    busco -i $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" -o "Predicted_Proteins_busco_molusca" -m prot -l mollusca
    mv "Predicted_Proteins_busco_molusca" $work_dir"/Mixed_Anot_Stringtie/BUSCO"
  fi
fi

################################################################################
################################################################################
################################################################################

if [ "$RUN_INTERPROT" == "TRUE" ]
then
  # conda activate Interprot
  rm -r $work_dir"/Interprot_Results"
  mkdir $work_dir"/Interprot_Results"
  mkdir $work_dir"/Interprot_Results/Sequence"

  seqkit seq -i $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" | sed "s/\*$//" >> $work_dir"/Interprot_Results/Sequence/work.pep"
  $interprot_path"/interproscan.sh" -f TSV --goterms -pa -b $work_dir"/Interprot_Results/Results" -dra -i $work_dir"/Interprot_Results/Sequence/work.pep"
fi

################################################################################
################################################################################
################################################################################

if [ "$RUN_EGGNOG" == "TRUE" ]
then
  rm -r $work_dir"/Eggnog"
  mkdir $work_dir"/Eggnog"

  emapper.py --override --cpu $threads -i $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" -o $work_dir"/Eggnog" --data_dir $eggnog_data
  mv $work_dir"/Eggnog."* $work_dir"/Eggnog"
  awk -F "\t" '{print $1"\t"$10}' $work_dir"/Eggnog/Eggnog.emapper.annotations" >> $work_dir"/Eggnog/Eggnog.emapper.go"
fi

################################################################################
################################################################################
################################################################################

if [ "$RUN_TRNA" == "TRUE" ]
then
  rm -r $work_dir"/tRNA_data"
  mkdir $work_dir"/tRNA_data"
  tRNAscan-SE -E -o $work_dir"/tRNA_data/trna_main_results.txt" -f $work_dir"/tRNA_data/trna_structure.txt" -j $work_dir"/tRNA_data/trna_anot.gff" $genome_file
fi

################################################################################
################################################################################
################################################################################

if [ "$RUN_RIBOSOME" == "TRUE" ]
then
  rm -r $work_dir"/Ribosome"
  mkdir $work_dir"/Ribosome"

  barrnap --kingdom euk --threads $threads --outseq $work_dir"/Ribosome/Ribosome_Seqs.fasta" $genome_file  > $work_dir"/Ribosome/Ribosomes.gff"
fi

################################################################################
################################################################################
################################################################################

if [ "$SEARCH_SIGNALS" == "TRUE" ]
then
  rm -r $work_dir"/Location_Signals"
  mkdir $work_dir"/Location_Signals"
  mkdir $work_dir"/Location_Signals/SignalP6"
  mkdir $work_dir"/Location_Signals/TargetP"

  seqkit seq -i $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" >> $work_dir"/Location_Signals/temp.fasta"

  signalp6 --organism euk --fastafile $work_dir"/Location_Signals/temp.fasta" --output_dir $work_dir"/Location_Signals/SignalP6"
  $targetp_file"/targetp" -fasta $work_dir"/Location_Signals/temp.fasta" -format "long" -prefix $work_dir"/Location_Signals/TargetP/Results"

  rm $work_dir"/Location_Signals/temp.fasta"
fi

################################################################################
################################################################################
################################################################################

if [ "$SEARCH_SECRETOMEP" == "TRUE" ]
then
  rm -r $work_dir"/Location_Signals/SECRETOMEP"
  mkdir $work_dir"/Location_Signals/SECRETOMEP"

  grep -v "Prediction" $work_dir"/Location_Signals/SignalP_results.txt" | grep -v "#" | awk -F "\t" '{print $1}' >> $work_dir"/Location_Signals/SECRETOMEP/temp.txt"
  grep -v "Prediction" $work_dir"/Location_Signals/Mitocondria_results.txt" | awk -F "\t" '{print $1}' >> $work_dir"/Location_Signals/SECRETOMEP/temp.txt"

  seqkit seq -n -i $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" | grep -v -f $work_dir"/Location_Signals/SECRETOMEP/temp.txt" >> $work_dir"/Location_Signals/SECRETOMEP/temp2.txt"
  seqkit grep -j $threads -f $work_dir"/Location_Signals/SECRETOMEP/temp2.txt" $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" | seqkit seq -i | tr -d "*" >> $work_dir"/Location_Signals/SECRETOMEP/Search_Secuences.fa"
  seqkit seq -n -i $work_dir"/Location_Signals/SECRETOMEP/Search_Secuences.fa" | grep -n . | awk -F ":" '{print "Gen-"$1"\t"$2}' >> $work_dir"/Location_Signals/SECRETOMEP/Convertion.id"

  count=1
  stop=$(grep -c . $work_dir"/Location_Signals/SECRETOMEP/Convertion.id")

  while [ $count -le $stop ]
  do
    echo $count" of "$stop
    original=$(sed -n $count"p" $work_dir"/Location_Signals/SECRETOMEP/Convertion.id" | awk -F "\t" '{print $2"$"}')
    replace=$(sed -n $count"p" $work_dir"/Location_Signals/SECRETOMEP/Convertion.id" | awk -F "\t" '{print $1}')
    sed -i "s/$original/$replace/" $work_dir"/Location_Signals/SECRETOMEP/Search_Secuences.fa"
    count=$(($count + 1))
  done

  $secretomep"/"secretomep -v -s $work_dir"/Location_Signals/SECRETOMEP/Search_Secuences.fa" >& $work_dir"/Location_Signals/SECRETOMEP/Secretome.out"

  rm $work_dir"/Location_Signals/SECRETOMEP/temp.txt" $work_dir"/Location_Signals/SECRETOMEP/temp2.txt"
fi
################################################################################
################################################################################
################################################################################

if [ "$DeepTMHMM" == "TRUE" ]
then
  rm -r biolib_results
  rm -r $work_dir"/DeepTMHMM"

  biolib run DTU/DeepTMHMM --fasta $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep"
  mv biolib_results $work_dir"/DeepTMHMM"
fi

################################################################################
################################################################################
################################################################################

if [ "$TABLE_TMHMM" == "TRUE" ]
then
  rm -r $work_dir"/DeepTMHMM/Membrane_candidates"
  mkdir $work_dir"/DeepTMHMM/Membrane_candidates"

  echo "ProtID Status N_Transmembrane SignalP SP_Prob SecretomeP NN-score TargetP TP_Prob" | tr " " "\t" >> $work_dir"/DeepTMHMM/Membrane_candidates/Summary_candidates.tab"
  check_prot=$(grep Length: $work_dir"/DeepTMHMM/TMRs.gff3" | sed "s/Length:.*//" | sed "s/# //")

  for prot in $check_prot
  do
    res_signalP=$(grep -F -w $prot $work_dir"/Location_Signals/SignalP6/prediction_results.txt" | awk -F "\t" '{print $2,$4}')
    res_targetP=$(grep -F -w $prot $work_dir"/Location_Signals/TargetP/Results_summary.targetp2" | awk -F "\t" '{print $2,$5}')
    check_secretomeP=$(grep -F -c -w $prot $work_dir"/Location_Signals/SECRETOMEP/Final/SecretomeP_work_output.tab")
    if [ $check_secretomeP -gt 0 ]
    then
      res_secretomeP=$(grep -F -w $prot $work_dir"/Location_Signals/SECRETOMEP/Final/SecretomeP_work_output.tab" | awk -F "\t" '{if ($2>=0.95) print "Secreted "$2; else print "NA "$2}')
    else
      res_secretomeP=$(echo "NA NA")
    fi

    number_TMHMM=$(grep $prot" Number of predicted TMRs:" $work_dir"/DeepTMHMM/TMRs.gff3" | awk -F ":" '{print $2}' | tr -d " ")

    # Status!!!
    check_signalp=$(echo $res_signalP | grep -c -w SP)
    check_targetp=$(echo $res_targetP | grep -c -w mTP)
    check_secretomeP=$(echo $res_secretomeP | grep -c -w Secreted)

    echo $prot
    echo $res_signalP": "$check_signalp
    echo $res_targetP": "$check_targetp
    echo $res_secretomeP": "$check_secretomeP
    echo ""

    if [ $number_TMHMM -eq 0 ]
    then
      if [ $check_signalp -gt 0 ] || [ $check_secretomeP -gt 0 ]
      then
        status=Secreted
        seqkit grep -j $threads -p $prot $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" >> $work_dir"/DeepTMHMM/Membrane_candidates/Secreted_prot.aa"
      elif [ $check_targetp -gt 0 ]
      then
        status=Mitochondia_interior
        seqkit grep -j $threads -p $prot $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" >> $work_dir"/DeepTMHMM/Membrane_candidates/Mitochondia_interior_prot.aa"
      else
        status=Other
      fi
    else
      if [ $check_signalp -gt 0 ] || [ $check_secretomeP -gt 0 ]
      then
        status=Plasmatic_Membrane
        seqkit grep -j $threads -p $prot $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" >> $work_dir"/DeepTMHMM/Membrane_candidates/Plasmatic_Membrane_prot.aa"
      elif [ $check_targetp -gt 0 ]
      then
        status=Mitochondia_membrane
        seqkit grep -j $threads -p $prot $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" >> $work_dir"/DeepTMHMM/Membrane_candidates/Mitochondia_membrane_prot.aa"
      else
        status=Potentially_membrane
        seqkit grep -j $threads -p $prot $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" >> $work_dir"/DeepTMHMM/Membrane_candidates/Potentially_membrane_prot.aa"
      fi
    fi

    echo $prot" "$status
    echo $prot" "$status" "$number_TMHMM" "$res_signalP" "$res_secretomeP" "$res_targetP | tr " " "\t" >> $work_dir"/DeepTMHMM/Membrane_candidates/Summary_candidates.tab"
  done
  seqkit stats -T $work_dir"/DeepTMHMM/Membrane_candidates/"*"_prot.aa" >> $work_dir"/DeepTMHMM/Membrane_candidates/Stats.tab"
fi

################################################################################
################################################################################
################################################################################

if [ "$SAM_IGSF_SEARCH" == "TRUE" ]
then
  rm -r $work_dir"/IGsF_Domain_Search"
  mkdir $work_dir"/IGsF_Domain_Search"
  seqkit seq -i $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" >> $work_dir"/IGsF_Domain_Search/Work.aa"
  hmmsearch -E 1e-10 --cpu $threads -o $work_dir"/IGsF_Domain_Search/General_Output.txt" --domtblout $work_dir"/IGsF_Domain_Search/Per_Domain_Table.tab" $igsf_hmmer_profile $work_dir"/IGsF_Domain_Search/Work.aa"

  echo "Prot_ID Prot_Accession Prot_len DomainID Prot_Accession Domain_Len Full_Seq_Evalue Full_Seq_Score Full_Seq_Bias N_Domain Total_Domains cEvalue i-Evalue Domain_Score Domain_Bias Hmm_Start Hmm_End Ali_Start Ali_End Env_Start Env_End acc Target_Description" | tr " " "\t" >> $work_dir"/IGsF_Domain_Search/Per_Domain_Table_work.tab"
  grep -v "#" $work_dir"/IGsF_Domain_Search/Per_Domain_Table.tab" | tr -s ' ' | tr " " "\t" >> $work_dir"/IGsF_Domain_Search/Per_Domain_Table_work.tab"
  rm $work_dir"/IGsF_Domain_Search/Work.aa"
fi


################################################################################
################################################################################
################################################################################


if [ "$SEARCH_CREP_FREP_GREP" == "TRUE" ]
then
  explore_orthofinder=$orthofinder_folder
  results_orthofinder=$(ls -d $explore_orthofinder"/OrthoFinder/Results_"*)

  rm -r $work_dir"/FREP_like_Search"
  mkdir $work_dir"/FREP_like_Search"
  mkdir $work_dir"/FREP_like_Search/BLAST_Reference"
  mkdir $work_dir"/FREP_like_Search/Overlap_temp"
  mkdir $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded"

  cat $Lu_FREP_prot $Dheilly_FREP_prot $Dheilly_CREP_prot $Dheilly_GREP_prot >> $work_dir"/FREP_like_Search/Reference_proteins.fasta"
  makeblastdb -in $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" -dbtype prot -parse_seqids -input_type fasta -out $work_dir"/FREP_like_Search/BLAST_Reference/Bsud_ref"
  blastp -query $work_dir"/FREP_like_Search/Reference_proteins.fasta" -db $work_dir"/FREP_like_Search/BLAST_Reference/Bsud_ref" -outfmt "6 std qcovs" -out $work_dir"/FREP_like_Search/FREP_Candidates_full.blastp" -qcov_hsp_perc 30 -num_threads $threads
  awk -F "\t" -v iden=$crep_frep_grep_blast_ident '{if ($3>=iden) print}' $work_dir"/FREP_like_Search/FREP_Candidates_full.blastp" >> $work_dir"/FREP_like_Search/FREP_Candidates.blastp"

  print_FBD_signatures=$(echo $FBD_signatures | tr " " "\n" | sed "s/^/FBD_Like:/" | tr "\n" " " | sed "s/ $//")
  print_C_lectin_signatures=$(echo $C_lectin_signatures | tr " " "\n" | sed "s/^/C-lectin_Like:/" | tr "\n" " " | sed "s/ $//")
  print_Galectin_signatures=$(echo $Galectin_signatures | tr " " "\n" | sed "s/^/Galectin_Like:/" | tr "\n" " " | sed "s/ $//")
  print_EGF_signatures=$(echo $EGF_signatures | tr " " "\n" | sed "s/^/EGF_Like:/" | tr "\n" " " | sed "s/ $//")
  print_Inmunoglobulin_signatures=$(echo $Inmunoglobulin_signatures | tr " " "\n" | sed "s/^/Inmunoglobulin_Like:/" | tr "\n" " " | sed "s/ $//")

  echo "GeneID ProtID Type N_Blast_Hits Best_Hit Identity Qcoverage Signalp Secretomep Targetp N_TMHMM Total_IGsF Non_Overlaping_IGsF IGsF1 IGsF2 "$print_FBD_signatures" "$print_C_lectin_signatures" "$print_Galectin_signatures" "$print_EGF_signatures" "$print_Inmunoglobulin_signatures | tr " " "\t" >> $work_dir"/FREP_like_Search/Candidates.tab"

  candidate_genes=$(awk -F "\t" '{print $2}' $work_dir"/FREP_like_Search/FREP_Candidates.blastp" | awk -F "." '{print $1"."$2}')

  # Genes with Fibrogen like Interprot signatures
  for sign in $FBD_signatures
  do
    add_genes=$(grep -w -F $sign $work_dir"/Interprot_Results/Results.tsv" | awk -F "\t" '{print $1}' | awk -F "." '{print $1"."$2}')
    candidate_genes=$(echo $candidate_genes" "$add_genes)
  done

  # Genes with c-lectin like Interprot signatures
  for sign in $C_lectin_signatures
  do
    add_genes=$(grep -w -F $sign $work_dir"/Interprot_Results/Results.tsv" | awk -F "\t" '{print $1}' | awk -F "." '{print $1"."$2}')
    candidate_genes=$(echo $candidate_genes" "$add_genes)
  done

  # Genes with Galactin like Interprot signatures
  for sign in $Galectin_signatures
  do
    add_genes=$(grep -w -F $sign $work_dir"/Interprot_Results/Results.tsv" | awk -F "\t" '{print $1}' | awk -F "." '{print $1"."$2}')
    candidate_genes=$(echo $candidate_genes" "$add_genes)
  done

  # Doing Clasification
  candidate_genes=$(echo $candidate_genes | tr " " "\n" | sort -u)
  for gen in $candidate_genes
  do
    transcripts=$(seqkit seq -n -i $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" | grep -F $gen"." | sort -u)
    for trans in $transcripts
    do
      grep -F -w $trans $work_dir"/Interprot_Results/Results.tsv" >> $work_dir"/FREP_like_Search/interport.temp"
      check_blast_hit=$(grep -c -w -F $trans $work_dir"/FREP_like_Search/FREP_Candidates.blastp")
      if [ $check_blast_hit -gt 0 ]
      then
        best_hit_info=$(grep -w -F $trans $work_dir"/FREP_like_Search/FREP_Candidates.blastp" | awk -F "\t" '{print $1"\t"$3"\t"$13"\t"$12*100}' | sort -n -k 4 | tail -n 1 | awk -F "\t" '{print $1"\t"$2"\t"$3}')
      else
        best_hit_info=$(echo "NA NA NA")
      fi

      res_signalp=$(grep -c -w -F $trans $work_dir"/Location_Signals/SignalP_results.txt" | awk '{if ($1>0) print "X"; else print "."}')
      res_secretomep=$(grep -c -w -F $trans $work_dir"/Location_Signals/SecretomeP_results.txt" | awk '{if ($1>0) print "X"; else print "."}')
      res_targetp=$(grep -c -w -F $trans $work_dir"/Location_Signals/Mitocondria_results.txt" | awk '{if ($1>0) print "X"; else print "."}')

      registro_Inmunoglobulin_signs=
      for sign in $Inmunoglobulin_signatures
      do
        check_hit=$(grep -c -w -F $sign $work_dir"/FREP_like_Search/interport.temp")
        overlap_signatures
        registro_Inmunoglobulin_signs=$(echo $registro_Inmunoglobulin_signs" "$sign_res)
      done

      registro_fibrinogen_signs=
      for sign in $FBD_signatures
      do
        check_hit=$(grep -c -w -F $sign $work_dir"/FREP_like_Search/interport.temp")
        overlap_signatures
        registro_fibrinogen_signs=$(echo $registro_fibrinogen_signs" "$sign_res)
      done

      registro_C_lectin_signs=
      for sign in $C_lectin_signatures
      do
        check_hit=$(grep -c -w -F $sign $work_dir"/FREP_like_Search/interport.temp")
        overlap_signatures
        registro_C_lectin_signs=$(echo $registro_C_lectin_signs" "$sign_res)
      done

      registro_Galectin_signs=
      for sign in $Galectin_signatures
      do
        check_hit=$(grep -c -w -F $sign $work_dir"/FREP_like_Search/interport.temp")
        overlap_signatures
        registro_Galectin_signs=$(echo $registro_Galectin_signs" "$sign_res)
      done

      registro_EGF_signatures=
      for sign in $EGF_signatures
      do
        check_hit=$(grep -c -w -F $sign $work_dir"/FREP_like_Search/interport.temp")
        overlap_signatures
        registro_EGF_signatures=$(echo $registro_EGF_signatures" "$sign_res)
      done

      # Check Sam's domains
      res_igsf1_hit=$(grep -F -w $trans $work_dir"/IGsF_Domain_Search/Per_Domain_Table_work.tab" | grep -c -w -F "IgSF1"| awk '{if ($1>0) print "X"; else print "."}')
      res_igsf2_hit=$(grep -F -w $trans $work_dir"/IGsF_Domain_Search/Per_Domain_Table_work.tab" | grep -c -w -F "IgSF2"| awk '{if ($1>0) print "X"; else print "."}')
      total_igsf_detections=$(grep -F -w $trans $work_dir"/IGsF_Domain_Search/Per_Domain_Table_work.tab" | grep -c .)
      grep -F -w $trans $work_dir"/IGsF_Domain_Search/Per_Domain_Table_work.tab" | awk -F "\t" '{print $20"\t"$21}' | sort -n -k 1 >> $work_dir"/FREP_like_Search/Overlap_temp/"$trans"_igsf_coords"

      if [ $total_igsf_detections -gt 1 ]
      then
        start_first_igsf=$(sed -n "1p" $work_dir"/FREP_like_Search/Overlap_temp/"$trans"_igsf_coords" | awk -F "\t" '{print $1}')
        end_first_igsf=$(sed -n "1p" $work_dir"/FREP_like_Search/Overlap_temp/"$trans"_igsf_coords" | awk -F "\t" '{print $2}')
        count=2
        non_overlaping_igsf=1
        while [ $count -le $total_igsf_detections ]
        do
          start_test_igsf=$(sed -n $count"p" $work_dir"/FREP_like_Search/Overlap_temp/"$trans"_igsf_coords" | awk -F "\t" '{print $1}')
          end_test_igsf=$(sed -n $count"p" $work_dir"/FREP_like_Search/Overlap_temp/"$trans"_igsf_coords" | awk -F "\t" '{print $2}')

          if [ $start_test_igsf -le $end_first_igsf -a $start_first_igsf -le $end_test_igsf ]
          then
            end_first_igsf=$end_test_igsf
            count=$(($count + 1))
          else
            start_first_igsf=$start_test_igsf
            end_first_igsf=$end_test_igsf
            non_overlaping_igsf=$(($non_overlaping_igsf + 1))
            count=$(($count + 1))
          fi
        done
      elif [ $total_igsf_detections -eq 1 ]
      then
        non_overlaping_igsf=1
      else
        non_overlaping_igsf=0
      fi

      # Check Trans-membrane domain
      Hit_DeepTMHMM=$(grep -F -w $trans" Number of predicted TMRs:" $work_dir"/DeepTMHMM/TMRs.gff3" | awk -F ":" '{print $2}' | tr -d " ")

      # Check components
      check_secretion=$(echo "$res_signalp" "$res_secretomep" | grep -c "X" )
      check_mito=$(echo $res_targetp | grep -c "X" )
      check_IGSF_1=$(echo $res_igsf1_hit | grep -c "X" )
      check_IGSF_2=$(echo $res_igsf2_hit | grep -c "X" )
      check_Interprot_InmunoGlobin=$(echo "$registro_Inmunoglobulin_signs" | grep -c "X" )
      check_Fibrogen=$(echo "$registro_fibrinogen_signs" | grep -c "X" )
      check_Clectin=$(echo "$registro_C_lectin_signs" | grep -c "X" )
      check_Galectin=$(echo "$registro_Galectin_signs" | grep -c "X" )
      check_EGF=$(echo "$registro_EGF_signatures" | grep -c "X" )

      check_functional_clusterfuck=$(($check_Clectin + $check_Fibrogen + $check_Galectin))
      integrity=ERROR
      if [ $Hit_DeepTMHMM -gt 0 ]
      then
        integrity=Membrane_Gene
      elif [ $check_secretion -gt 0 ]
      then
        if [ $check_IGSF_1 -gt 0 ] || [ $check_IGSF_2 -gt 0 ]
        then
          integrity=Full
        else
          if [ $check_Interprot_InmunoGlobin -gt 0 ]
          then
            integrity=Divergent_IGSF
          else
            integrity=No_IGSF
          fi
        fi
      elif [ $check_secretion -eq 0 ]
      then
        if [ $check_IGSF_1 -gt 0 ] || [ $check_IGSF_2 -gt 0 ]
        then
          integrity=NO_SP
        else
          if [ $check_Interprot_InmunoGlobin -gt 0 ]
          then
            integrity=No_SP_Divergent_IGSF
          else
            integrity=NO_SP_NO_IGSF
          fi
        fi
      fi

      if [ "$integrity" == "Membrane_Gene" ]
      then
        Gene_putative_Type=$(echo "Membrane_Gene")
      elif [ "$integrity" == "NO_SP_NO_IGSF" ]
      then
        if [ $check_EGF -eq 1 ]
        then
          if [ $check_Fibrogen -eq 1 ]
          then
            Gene_putative_Type=$(echo "FREM")
          else
            Gene_putative_Type=$(echo "Protein_with_EGF_Domain")
          fi
        else
          Gene_putative_Type=$(echo "Other:"$integrity)
        fi
      else
        if [ $check_functional_clusterfuck -gt 1 ]
        then
          Gene_putative_Type=Multiple_signature_Domains
        else
          if [ $check_EGF -eq 1 ]
          then
            if [ $check_Fibrogen -eq 1 ]
            then
              Gene_putative_Type=$(echo "FREM")
            else
              Gene_putative_Type=$(echo "Protein_with_EGF_Domain")
            fi
          elif [ $check_Fibrogen -eq 1 ]
          then
            interest_signatures=$FBD_signatures
            Type=FREP
            retrive_domain_locations
            Gene_putative_Type=$(echo "FREP:"$integrity)
          elif [ $check_Clectin -eq 1 ]
          then
            interest_signatures=$C_lectin_signatures
            Type=CREP
            retrive_domain_locations
            Gene_putative_Type=$(echo "CREP:"$integrity)
          elif [ $check_Galectin -eq 1 ]
          then
            interest_signatures=$Galectin_signatures
            Type=GREP
            retrive_domain_locations
            Gene_putative_Type=$(echo "GREP:"$integrity)
          else
            Gene_putative_Type=$(echo "Other:No_key_domain")
          fi
        fi
      fi
      echo $gen" "$trans" "$Gene_putative_Type
      echo $gen" "$trans" "$Gene_putative_Type" "$check_blast_hit" "$best_hit_info" "$res_signalp" "$res_secretomep" "$res_targetp" "$Hit_DeepTMHMM" "$total_igsf_detections" "$non_overlaping_igsf" "$res_igsf1_hit" "$res_igsf2_hit" "$registro_fibrinogen_signs" "$registro_C_lectin_signs" "$registro_Galectin_signs" "$registro_EGF_signatures" "$registro_Inmunoglobulin_signs | tr " " "\t" >> $work_dir"/FREP_like_Search/Candidates.tab"
      rm $work_dir"/FREP_like_Search/interport.temp"
    done
  done
fi

################################################################################
################################################################################
################################################################################

if [ "$CREP_FREP_PHYLOGENY" == "TRUE" ]
then
  explore_orthofinder=$orthofinder_folder
  results_orthofinder=$(ls -d $explore_orthofinder"/OrthoFinder/Results_"*)

  rm -r $work_dir"/FREP_like_Search/Filogeny"
  mkdir $work_dir"/FREP_like_Search/Filogeny"
  mkdir $work_dir"/FREP_like_Search/Filogeny/FREP"
  mkdir $work_dir"/FREP_like_Search/Filogeny/FREP/Luetal_ref"
  mkdir $work_dir"/FREP_like_Search/Filogeny/FREP/Dheilly_ref"
  mkdir $work_dir"/FREP_like_Search/Filogeny/FREP/ALL_ref"
  mkdir $work_dir"/FREP_like_Search/Filogeny/CREP"
  mkdir $work_dir"/FREP_like_Search/Filogeny/CREP/Dheilly_ref"

  if [ -f $exclude_seq_manual]
  then
    echo "Excluding sequences on: "$exclude_seq_manual
    genes_Bsud=$(seqkit seq -n $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/FREP_Full_candidates.aa" | awk -F "." '{print $1"."$2}' | sort -u | grep -w -F -v -f $exclude_seq_manual)
  else
    echo "No file with manually excluded sequences at: "$exclude_seq_manual
    genes_Bsud=$(seqkit seq -n $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/FREP_Full_candidates.aa" | awk -F "." '{print $1"."$2}' | sort -u)
  fi

  for gen in $genes_Bsud
  do
    echo $gen" FREP"
    selected=$(seqkit grep -r -p $gen"\." $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/FREP_Full_candidates.aa" | seqkit fx2tab -n -l | sort -r -n -k 2 | head -n 1 | awk -F "\t" '{print $1}')
    seqkit grep -p $selected $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/FREP_Full_candidates.aa" >> $work_dir"/FREP_like_Search/Filogeny/FREP/FREP_Bsud.fasta"

    pho=$(grep -F $gen"." $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" | awk -F "\t" '{print $1}' | sort -u)
    echo $gen" "$pho | tr " " "\t" >> $work_dir"/FREP_like_Search/Filogeny/FREP/gene_pho.txt"
  done

  # 1) FREP Lu Reference
  storage_folder=$work_dir"/FREP_like_Search/Filogeny/FREP/Luetal_ref"
  grupo_phy=FREP_Lu
  genes_ref=$(seqkit seq -n $Lu_FREP_prot | sed "s/_--_/\t/" | awk -F "\t" '{print $1}' | awk -F "-" '{print $1}' | sort -u)
  for gen in $genes_ref
  do
    echo $gen" Lu"
    selected=$(seqkit grep -r -p $gen"-" $Lu_FREP_prot | seqkit fx2tab -n -l | sort -r -n -k 2 | head -n 1 | awk -F "\t" '{print $1}')
    seqkit grep -p $selected $Lu_FREP_prot >> $storage_folder"/Reference.fasta"
  done
  cat $work_dir"/FREP_like_Search/Filogeny/FREP/FREP_Bsud.fasta" $storage_folder"/Reference.fasta" >> $storage_folder"/Run_Analysis.fasta"
  run_frep_crep_phylogeny

  # 2) FREP Dheilly Reference
  storage_folder=$work_dir"/FREP_like_Search/Filogeny/FREP/Dheilly_ref"
  grupo_phy=FREP_Dheilly
  cat $Dheilly_FREP_prot $work_dir"/FREP_like_Search/Filogeny/FREP/FREP_Bsud.fasta" >> $storage_folder"/Run_Analysis.fasta"
  run_frep_crep_phylogeny

  # 3) FREP All
  storage_folder=$work_dir"/FREP_like_Search/Filogeny/FREP/ALL_ref"
  grupo_phy=FREP_All
  cat $work_dir"/FREP_like_Search/Filogeny/FREP/FREP_Bsud.fasta" $work_dir"/FREP_like_Search/Filogeny/FREP/Luetal_ref/Reference.fasta" $Dheilly_FREP_prot >> $storage_folder"/Run_Analysis.fasta"
  run_frep_crep_phylogeny

  # 4) CREP
  genes_Bsud=$(seqkit seq -n $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/CREP_Full_candidates.aa" | awk -F "." '{print $1"."$2}' | sort -u | grep -w -F -v -f $exclude_seq_manual)
  for gen in $genes_Bsud
  do
    echo $gen" CREP"
    selected=$(seqkit grep -r -p $gen"\." $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/CREP_Full_candidates.aa" | seqkit fx2tab -n -l | sort -r -n -k 2 | head -n 1 | awk -F "\t" '{print $1}')
    seqkit grep -p $selected $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/CREP_Full_candidates.aa" >> $work_dir"/FREP_like_Search/Filogeny/CREP/CREP_Bsud.fasta"

    pho=$(grep -F $gen"." $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" | awk -F "\t" '{print $1}' | sort -u)
    echo $gen" "$pho | tr " " "\t" >> $work_dir"/FREP_like_Search/Filogeny/FREP/gene_pho.txt"
  done

  storage_folder=$work_dir"/FREP_like_Search/Filogeny/CREP/Dheilly_ref"
  grupo_phy=CREP
  cat $work_dir"/FREP_like_Search/Filogeny/CREP/CREP_Bsud.fasta" $Dheilly_CREP_prot >> $storage_folder"/Run_Analysis.fasta"
  run_frep_crep_phylogeny
fi

################################################################################
################################################################################
################################################################################

if [ "$SUP_TAB_PAPER_CREP_FREP_GREP" == "TRUE" ]
then
  rm $work_dir"/Final_Anotation_Files/Paper_Special_Tables/prot.temp"
  rm $work_dir"/Final_Anotation_Files/Paper_Special_Tables/CREP_FREP_GREP_Sup_Table.txt"
  rm $work_dir"/Final_Anotation_Files/Paper_Special_Tables/interpro.temp"

  families=$(echo "CREP FREP GREP")
  echo "Protein_ID Family Final_Selection Analysis Signature_accession Score Start End InterPro_annotation Description Signature_Type Final_Selection" | tr " " "\t"  >> $work_dir"/Final_Anotation_Files/Paper_Special_Tables/CREP_FREP_GREP_Sup_Table.txt"
  for fam in $families
  do
    protein_candidates=$(grep "Protein:" $work_dir"/FREP_like_Search/Best_Candidates_ReviewNeded/"$fam"_candidates_domain_locations.tab" | awk -F " " '{print $2}')

    for prot in $protein_candidates
    do
      echo $fam": "$prot

      if [ "$fam" != "GREP" ]
      then
        is_final_set=$(grep -w -F -c $prot $work_dir"/FREP_like_Search/Filogeny/"$fam"/"$fam"_Bsud.fasta" | awk '{if ($1>0) print "X"; else print "."}')
      else
        is_final_set=.
      fi

      check_signal_peptide=$(grep -w -F -c $prot $work_dir"/Location_Signals/SignalP_results.txt")
      if [ $check_signal_peptide -gt 0 ]
      then
        sp_coords=$(grep -w -F $prot $work_dir"/Location_Signals/SignalP_results.txt"  | awk -F "\t" '{print $5}' | awk -F " " '{print $3}' | sed "s/.$//" | tr "-" ";")
        probability=$(grep -w -F $prot $work_dir"/Location_Signals/SignalP_results.txt"  | awk -F "\t" '{print $4}')
        echo $prot";"$fam";"$is_final_set";SignalP;SP;"$probability";"$sp_coords";NA;NA;Secreted Protein" | tr ";" "\t" >> $work_dir"/Final_Anotation_Files/Paper_Special_Tables/prot.temp"
      fi

      check_igsf=$(grep -w -F -c $prot $work_dir"/IGsF_Domain_Search/Per_Domain_Table_work.tab" )
      if [ $check_igsf -gt 0 ]
      then
        grep -w -F $prot $work_dir"/IGsF_Domain_Search/Per_Domain_Table_work.tab" > $work_dir"/Final_Anotation_Files/Paper_Special_Tables/igsf.temp"
        count_igsf=1
        recorrer_igsf=$(grep -c . $work_dir"/Final_Anotation_Files/Paper_Special_Tables/igsf.temp")

        while [ $count_igsf -le $recorrer_igsf ]
        do
          domain=$(sed -n $count_igsf"p" $work_dir"/Final_Anotation_Files/Paper_Special_Tables/igsf.temp" | awk -F "\t" '{print $4}')
          probability=$(sed -n $count_igsf"p" $work_dir"/Final_Anotation_Files/Paper_Special_Tables/igsf.temp" | awk -F "\t" '{print $12}')
          igsf_coords=$(sed -n $count_igsf"p" $work_dir"/Final_Anotation_Files/Paper_Special_Tables/igsf.temp" | awk -F "\t" '{print $18";"$19}')
          echo $prot";"$fam";"$is_final_set";Hmmsearch;"$domain";"$probability";"$igsf_coords";NA;NA;IgSF like Domain" | tr ";" "\t" >> $work_dir"/Final_Anotation_Files/Paper_Special_Tables/prot.temp"

          count_igsf=$(($count_igsf + 1))
        done
      fi

      grep $prot $work_dir"/Interprot_Results/Results.tsv" | awk -F "\t" -v family=$fam -v is_final=$is_final_set '{print $1";"family";"is_final";"$4";"$5";"$9";"$7";"$8";"$12";"$13}' | tr ";" "\t"  >> $work_dir"/Final_Anotation_Files/Paper_Special_Tables/interpro.temp"

      count=1
      recorrer=$(grep -c . $work_dir"/Final_Anotation_Files/Paper_Special_Tables/interpro.temp")

      while [ $count -le $recorrer ]
      do
        current_interpro=$(sed -n $count"p" $work_dir"/Final_Anotation_Files/Paper_Special_Tables/interpro.temp" | awk -F "\t" '{print $9}')

        if [ $current_interpro == "-" ]
        then
          signature=Other
        else
          check_FBD=$(echo $FBD_signatures | grep -w -c $current_interpro)
          check_immuno=$(echo $Inmunoglobulin_signatures | grep -w -c $current_interpro)
          check_clectin=$(echo $C_lectin_signatures | grep -w -c $current_interpro)
          check_galectin=$(echo $Galectin_signatures | grep -w -c $current_interpro)

          if [ $check_FBD -gt 0 ]
          then
            signature=$(echo "FBD like Domain")
          elif [ $check_immuno -gt 0 ]
          then
            signature=$(echo "Immunoglobulin like Domain")
          elif [ $check_clectin -gt 0 ]
          then
            signature=$(echo "C-Lectin like Domain")
          elif [ $check_galectin -gt 0 ]
          then
            signature=$(echo "Galectin like Domain")
          else
            signature=Other
          fi

          sed -n $count"p" $work_dir"/Final_Anotation_Files/Paper_Special_Tables/interpro.temp" | awk -v add="$signature" '{print $0"\t"add}' >> $work_dir"/Final_Anotation_Files/Paper_Special_Tables/prot.temp"
        fi
        count=$(($count + 1))
      done

      sort -n -k 7 $work_dir"/Final_Anotation_Files/Paper_Special_Tables/prot.temp" >> $work_dir"/Final_Anotation_Files/Paper_Special_Tables/CREP_FREP_GREP_Sup_Table.txt"

      rm $work_dir"/Final_Anotation_Files/Paper_Special_Tables/prot.temp"

      if [ -f $work_dir"/Final_Anotation_Files/Paper_Special_Tables/interpro.temp" ]
      then
        rm $work_dir"/Final_Anotation_Files/Paper_Special_Tables/interpro.temp"
      fi

      if [ -f $work_dir"/Final_Anotation_Files/Paper_Special_Tables/igsf.temp" ]
      then
        rm $work_dir"/Final_Anotation_Files/Paper_Special_Tables/igsf.temp"
      fi
    done
  done
fi

################################################################################
################################################################################
################################################################################

if [ "$DATA_CROSS_INTERPROT_SIGNALP" == "TRUE" ]
then
  storage_place=$work_dir"/Data_Crossroads/Secretion_Interpro.txt"
  genes_interes=$(grep -v -w "Prediction" $work_dir"/Location_Signals/SignalP_results.txt" | awk -F "\t" '{print $1}' | sort -u)
  special_header=$(echo "Prediction;OTHER;SP(Sec/SPI);CS Position")
  extra_data_file=$work_dir"/Location_Signals/SignalP_results.txt"
  interprot_basic_mining

  storage_place=$work_dir"/Data_Crossroads/Mitochondria_Interpro.txt"
  genes_interes=$(grep -v -w "Prediction" $work_dir"/Location_Signals/Mitocondria_results.txt" | awk -F "\t" '{print $1}' | sort -u)
  special_header=$(echo "ID;Prediction;noTP;SP;mTP;CS Position")
  extra_data_file=$work_dir"/Location_Signals/Mitocondria_results.txt"
  interprot_basic_mining
fi

################################################################################
################################################################################
################################################################################

if [ "$SUMMARY" == "TRUE" ]
then
  rm -r $work_dir"/Final_Anotation_Files"
  mkdir $work_dir"/Final_Anotation_Files"
  mkdir $work_dir"/Final_Anotation_Files/PolyA_Annotation"
  mkdir $work_dir"/Final_Anotation_Files/PolyA_Annotation/Transdecoter_intermediary_files"
  mkdir $work_dir"/Final_Anotation_Files/Repeats"
  mkdir $work_dir"/Final_Anotation_Files/BUSCO"

  # Annotation of vanilla genes (protein coding)
  cp $work_dir"/TransDecoder/Bsudanica_annotation_polyA.gff3" $work_dir"/Final_Anotation_Files/PolyA_Annotation/Bsudanica_annotation_polyA.gff3"
  cp $work_dir"/TransDecoder/Bsudanica_annotation_mapped_prot.fasta" $work_dir"/Final_Anotation_Files/PolyA_Annotation/Bsudanica_annotation_mapped_prot.fas"
  cp $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.cds" $work_dir"/Final_Anotation_Files/PolyA_Annotation/Transdecoter_intermediary_files/Bsudanica_CDS.fas"
  cp $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" $work_dir"/Final_Anotation_Files/PolyA_Annotation/Transdecoter_intermediary_files/Bsudanica_prot.faa"

  # Repeated elements
  cp -r $work_dir"/"$species_name"EarlGrey/Bsudanica_summaryFiles" $work_dir"/Final_Anotation_Files/Repeats/Summary_Files"
  cp -r $work_dir"/"$species_name"EarlGrey/Bsudanica_RepeatMasker_Against_Custom_Library" $work_dir"/Final_Anotation_Files/Repeats/RepeatMasker_Custom_Library"
  cp -r $work_dir"/"$species_name"EarlGrey/Bsudanica_clusTErs/Bsudanica.detailedClusTErs.bed" $work_dir"/Final_Anotation_Files/Repeats/RepeatMasker_Custom_Library/Reapeated_Elements_Clusters.bed"

  # Busco Results
  extract_things=$(ls -d $work_dir"/Mixed_Anot_Stringtie/BUSCO/"*)
  for ext in $extract_things
  do
    cp $ext"/"*txt $work_dir"/Final_Anotation_Files/BUSCO"
  done
fi

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

if [ "$PREPARE_ORTHOFINDER" == "TRUE" ]
then
  rm -r $orthofinder_folder
  mkdir $orthofinder_folder
  mkdir $orthofinder_folder"/Get_Orfs"
  mkdir $orthofinder_folder"/Mito_Exclusion"
  mkdir $orthofinder_folder"/Work_cDNA"

  work_species=$(ls $other_species_orthofinder"/"*".cds" | sed "s/.cds//" | sed "s/.*\///")

  makeblastdb -in $other_species_mitochondria -dbtype nucl -parse_seqids -input_type fasta -out $orthofinder_folder"/Mito_Exclusion/Reference_mit"

  for spe in $work_species
  do
    echo "Peparing: "$spe
    mkdir $orthofinder_folder"/Get_Orfs/"$spe
    mkdir $orthofinder_folder"/Get_Orfs/"$spe"/Mitocondria"

    check_dash=$(echo $dash_ids | grep -c $spe)
    if [ $check_dash -eq 0 ]
    then
      transcripts=$(seqkit seq -n -i $other_species_orthofinder"/"$spe".cds" | awk -F "." '{print $1}' | sort -u)
    else
      transcripts=$(seqkit seq -n -i $other_species_orthofinder"/"$spe".cds" | awk -F "-" '{print $1}' | sort -u)
    fi

    echo "Mitochondra BLAST: "$spe
    blastn -query $other_species_orthofinder"/"$spe".cds" -db $orthofinder_folder"/Mito_Exclusion/Reference_mit" -outfmt "6 std" -out $orthofinder_folder"/Get_Orfs/"$spe"/Mitocondria/"$spe_mito".blastn" -perc_identity 70 -num_threads $threads -qcov_hsp_perc 50

    check_hits=$(grep -c . $orthofinder_folder"/Get_Orfs/"$spe"/Mitocondria/"$spe_mito".blastn")
    if [ $check_hits -gt 0 ]
    then
      awk -F "\t" '{print $1}' $orthofinder_folder"/Get_Orfs/"$spe"/Mitocondria/"$spe_mito".blastn" >> $orthofinder_folder"/Get_Orfs/"$spe"/Exclude.temp"
      seqkit grep -j $threads -v -f $orthofinder_folder"/Get_Orfs/"$spe"/Exclude.temp" $other_species_orthofinder"/"$spe".cds" | seqkit seq -i >> $orthofinder_folder"/Get_Orfs/"$spe"/"$spe"_temp.cds"
    else
      seqkit seq -i $other_species_orthofinder"/"$spe".cds" >> $orthofinder_folder"/Get_Orfs/"$spe"/"$spe"_temp.cds"
    fi

    # New issue! Transdecoder doesn't works if you already removed the UTRs, as is the for most of these sequeces.
    # I will keep it simple since in theory the work is already done
    getorf -sequence $orthofinder_folder"/Get_Orfs/"$spe"/"$spe"_temp.cds" -outseq $orthofinder_folder"/Get_Orfs/"$spe"/"$spe"_temp.orf" -noreverse -minsize 150

    for trans in $transcripts
    do
      echo "Procesing: "$trans" --- "$spe
      if [ $check_dash -eq 0 ]
      then
        check_orf=$(grep -F -c $trans"." $orthofinder_folder"/Get_Orfs/"$spe"/"$spe"_temp.orf" )
      else
        check_orf=$(grep -F -c $trans"-" $orthofinder_folder"/Get_Orfs/"$spe"/"$spe"_temp.orf" )
      fi

      if [ $check_orf -gt 0 ]
      then
        if [ $check_dash -eq 0 ]
        then
          seqkit grep -j $threads -r -p $trans\\. $orthofinder_folder"/Get_Orfs/"$spe"/"$spe"_temp.orf" | seqkit fx2tab -n -l | tr -d " " | sort -n -k 2 > $orthofinder_folder"/Get_Orfs/orf_len.tmp"
        else
          seqkit grep -j $threads -r -p $trans"-" $orthofinder_folder"/Get_Orfs/"$spe"/"$spe"_temp.orf" | seqkit fx2tab -n -l | tr -d " " | sort -n -k 2 > $orthofinder_folder"/Get_Orfs/orf_len.tmp"
        fi

        best_len=$(tail -n 1 $orthofinder_folder"/Get_Orfs/orf_len.tmp" | awk -F "\t" '{print $2}')
        check_num_orf=$(awk -F "\t" -v len=$best_len '{if ($2==len) print}' $orthofinder_folder"/Get_Orfs/orf_len.tmp" | grep -c .)

        if [ $check_num_orf -eq 1 ]
        then
          best_orf=$(awk -F "\t" -v len=$best_len '{if ($2==len) print $1}' $orthofinder_folder"/Get_Orfs/orf_len.tmp" | tr -d " " | awk -F "\t" '{print $1}')
        else
          echo $trans >> $orthofinder_folder"/Get_Orfs/"$spe"/Trouble_Makers.txt"
          best_coord=$(awk -F "\t" -v len=$best_len '{if ($2==len) print $1}' $orthofinder_folder"/Get_Orfs/orf_len.tmp" | tr -d " " | awk -F "\t" '{print $1}' | sed "s/.*\[//" | tr -d "]" | tr "-" "\t" | sort -n -k 1 | head -n 1 | tr "\t" "-")
          best_orf=$(awk -F "\t" -v len=$best_len '{if ($2==len) print $1}' $orthofinder_folder"/Get_Orfs/orf_len.tmp" | tr -d " " | awk -F "\t" '{print $1}' | grep -w $best_coord | head -n 1)
        fi

        best_orf_id=$(echo $best_orf | sed "s/\[.*//")
        remove_from_id=$(echo $best_orf_id | tr "_" "\n" | tail -n 1 | awk '{print "_"$1"END"}')
        print_id=$(echo $best_orf | sed "s/\[.*//" | awk '{print $1"END"}' | sed "s/$remove_from_id//")

        coords=$(echo $best_orf | sed "s/.*\[//" | sed "s/\]//" | tr "-" ":")

        seqkit grep -j $threads -p $best_orf_id $orthofinder_folder"/Get_Orfs/"$spe"/"$spe"_temp.orf" | sed "s/>.*/>$print_id/" >> $orthofinder_folder"/"$spe".pep"
        seqkit grep -j $threads -p $print_id $orthofinder_folder"/Get_Orfs/"$spe"/"$spe"_temp.cds" | seqkit subseq -r $coords | seqkit seq -i >> $orthofinder_folder"/Work_cDNA/"$spe".cds"
      else
        echo $trans >> $orthofinder_folder"/Get_Orfs/"$spe"/Genes_no_Good_orf.txt"
      fi
    done
  done

  seqkit fx2tab -n -i -l $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" >> $orthofinder_folder"/lenght.temp"
  extract_peps=$(awk -F "\t" '{print $1}' $orthofinder_folder"/lenght.temp" | awk -F "." '{print $1"."$2}' | sort -u)
  for ext in $extract_peps
  do
    take=$(grep -w $ext $orthofinder_folder"/lenght.temp" | sort -n -k 2 | tail -n 1 | awk -F "\t" '{print $1}')
    echo $ext" "$take
    seqkit grep -j $threads -p $take $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.cds" | seqkit seq -i | seqkit seq -w 0 | sed "s/TAA$//"  | sed "s/TGA$//"  | sed "s/TAG$//" | seqkit seq -w 60 >> $orthofinder_folder"/Work_cDNA/Biomphalaria_sudanica.cds"
    seqkit grep -j $threads -p $take $work_dir"/TransDecoder/Predict/PolyA_Transcripts_work.fasta.transdecoder.pep" | seqkit seq -i | sed "s/\*$//" >> $orthofinder_folder"/Biomphalaria_sudanica.pep"
  done
  rm $orthofinder_folder"/lenght.temp"
fi

################################################################################
################################################################################
################################################################################

if [ "$ORTHOFINDER" == "TRUE" ]
then
  rm -r $orthofinder_folder"/OrthoFinder"
  orthofinder -f $orthofinder_folder -t $threads -a $threads
fi

################################################################################
################################################################################
################################################################################

if [ "$ORTHOFINDER_EGGNOG" == "TRUE" ]
then
  working_folder=$orthofinder_folder

  rm -r $working_folder"/eggnog_Orthofinder"
  mkdir $working_folder"/eggnog_Orthofinder"

  seqkit seq -i $working_folder"/"*".pep" >> $working_folder"/eggnog_Orthofinder/work_seqs.fa"
  emapper.py --override --cpu $threads -i  $working_folder"/eggnog_Orthofinder/work_seqs.fa" -o $work_dir"/Eggnog" --data_dir $eggnog_data

  mv $work_dir"/Eggnog."* $working_folder"/eggnog_Orthofinder"
  awk -F "\t" '{print $1"\t"$10}' $working_folder"/eggnog_Orthofinder/Eggnog.emapper.annotations" >> $working_folder"/eggnog_Orthofinder/Eggnog.emapper.go"
fi

################################################################################
################################################################################
################################################################################

if [ "$ORTHOFINDER_EGGNOG_FOR_PHO" == "TRUE" ]
then
  working_folder=$orthofinder_folder
  results_orthofinder=$(ls -d $orthofinder_folder"/OrthoFinder/Results_"*)

  rm -r $working_folder"/eggnog_Orthofinder/For_Enrichment"
  mkdir $working_folder"/eggnog_Orthofinder/For_Enrichment"

  count=2
  recorrer=$(grep -c . $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")

  while [ $count -le $recorrer ]
  do
    pho=$(sed -n $count"p" $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" | awk -F "\t" '{print $1}')
    sed -n $count"p" $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"  | tr "\t" "\n" | sed 1,3d | tr -d " " | tr "," "\n" | grep . > $working_folder"/eggnog_Orthofinder/For_Enrichment/Temp_file"
    go_terms=$(grep -w -F -f $working_folder"/eggnog_Orthofinder/For_Enrichment/Temp_file" $working_folder"/eggnog_Orthofinder/Eggnog.emapper.go" | awk -F "\t" '{print $2}' | tr "," "\n" | sort -u | grep GO | tr "\n" "," | sed "s/,$//")

    echo $pho" "$count"/"$recorrer
    echo $pho" "$go_terms | tr " " "\t" >> $working_folder"/eggnog_Orthofinder/For_Enrichment/Background_PHO_All.txt"
    count=$(($count + 1))
  done
fi

################################################################################
################################################################################
################################################################################

if [ "$CAFE_EXPANSION_PREPARE" == "TRUE" ]
then
  explore_orthofinder=$orthofinder_folder
  results_orthofinder=$(ls -d $explore_orthofinder"/OrthoFinder/Results_"*)
  prepare_cafe
fi

################################################################################
################################################################################
################################################################################

if [ "$CAFE_EXPANSION_TEST_RUNS" == "TRUE" ]
then
  explore_orthofinder=$orthofinder_folder
  results_orthofinder=$(ls -d $explore_orthofinder"/OrthoFinder/Results_"*)
  run_cafe
fi

################################################################################
################################################################################
################################################################################

if [ "$CAFE_EXPANSION_DEF_RUN" == "TRUE" ]
then
  explore_orthofinder=$orthofinder_folder
  results_orthofinder=$(ls -d $explore_orthofinder"/OrthoFinder/Results_"*)
  interest_nodes_file=$explore_orthofinder"/Interest_nodes_CAFE.txt"
  cafe_def_run
fi

################################################################################
################################################################################
################################################################################

if [ "$CAFE_EXPANSION_ANOT_EGGNOG" == "TRUE" ]
then
  working_folder=$orthofinder_folder"/CAFE_Expansion/DEF_Run"

  rm -r $working_folder"/Eggnog_Species"
  mkdir $working_folder"/Eggnog_Species"

  seq_files=$(ls $working_folder"/Sequences/"*".pep" | sed "s/_interest.pep//" | sed "s/.*\///")

  for seq in $seq_files
  do
    mkdir $working_folder"/Eggnog_Species/"$seq

    emapper.py --override --cpu $threads -i $working_folder"/Sequences/"$seq"_interest.pep" -o $working_folder"/Eggnog_Species/"$seq"/"$seq --data_dir $eggnog_data
    awk -F "\t" '{print $1"\t"$10}' $working_folder"/Eggnog_Species/"$seq"/"$seq".emapper.annotations"  >> $working_folder"/Eggnog_Species/"$seq"/"$seq".emapper.go"
  done
fi

################################################################################
################################################################################
################################################################################

if [ "$LOCATION_TABLE" == "TRUE" ]
then
  rm -r $work_dir"/Final_Anotation_Files/Paper_Special_Tables"
  mkdir $work_dir"/Final_Anotation_Files/Paper_Special_Tables"

  signalp_prot=$(grep -w -v "Prediction" $work_dir"/Location_Signals/SignalP_results.txt"| grep -v "#" | awk -F "\t" '{print $1}')
  targetp_prot=$(grep -w -v "Prediction" $work_dir"/Location_Signals/Mitocondria_results.txt" | grep -v "#" | awk -F "\t" '{print $1}')
  secretomep_prot=$(grep -w -v "Weighted_by_prior" $work_dir"/Location_Signals/SecretomeP_results.txt" | awk -F "\t" '{print $1}')

  work_genes=$(echo $signalp_prot" "$targetp_prot" "$secretomep_prot | tr " " "\n" | sort -u)
  echo "ProtID Type Program Signal Prob/Score Details N_TMHMM Interprot_Data" | tr " " "\t" >> $work_dir"/Final_Anotation_Files/Paper_Special_Tables/Location_Signals_paper.tab"
  echo $work_genes | tr " " "\n" >> $work_dir"/Final_Anotation_Files/Paper_Special_Tables/genes.temp"
  grep -w -F -f $work_dir"/Final_Anotation_Files/Paper_Special_Tables/genes.temp" $work_dir"/Interprot_Results/Results.tsv" | awk -F "\t" '{if ($12!="-") print $1"\t"$12"("$13")"}' | sort -u >> $work_dir"/Final_Anotation_Files/Paper_Special_Tables/interprot.temp"

  for prot in $work_genes
  do
    echo $prot
    check_signalP=$(grep -F -c -w $prot $work_dir"/Location_Signals/SignalP_results.txt")
    check_targetP=$(grep -F -c -w $prot $work_dir"/Location_Signals/Mitocondria_results.txt")
    check_secretomeP=$(grep -F -c -w $prot $work_dir"/Location_Signals/SecretomeP_results.txt")

    number_TMHMM=$(grep $prot" Number of predicted TMRs:" $work_dir"/DeepTMHMM/TMRs.gff3" | awk -F ":" '{print $2}' | tr -d " ")
    interprot_data=$(grep -w -F $prot $work_dir"/Final_Anotation_Files/Paper_Special_Tables/interprot.temp" | awk -F "\t" '{print $2}' | sort -u | sed "s/ /_ooo_/g" | tr "\n" ";" | sed "s/;$//" )

    echo $prot": "$check_signalP"||"$check_targetP"||"$check_secretomeP

    if [ $check_signalP -gt 0 ]
    then
      program=SignalP
      res_signal=$(grep -F -w $prot $work_dir"/Location_Signals/SignalP6/prediction_results.txt" | awk -F "\t" '{print $2"_x_"$4"_x_"$5}' | sed "s/ /_ooo_/g")

      if [ $number_TMHMM -eq 0 ]
      then
        type=Secreted
      else
        type=Membrane
      fi
    elif [ $check_targetP -gt 0 ]
    then
      program=TargetP
      res_signal=$(grep -F -w $prot $work_dir"/Location_Signals/TargetP/Results_summary.targetp2" | awk -F "\t" '{print $2"_x_"$5"_x_"$6}' | sed "s/ /_ooo_/g")

      if [ $number_TMHMM -eq 0 ]
      then
        type=Mitochondrial_Cytoplasm
      else
        type=Mitochondrial_Membrane
      fi
    elif [ $check_secretomeP -gt 0 ]
    then
      program=SecretomeP
      res_signal=$(grep -F -w $prot $work_dir"/Location_Signals/SECRETOMEP/Final/SecretomeP_work_output.tab" | awk -F "\t" '{print "NA_x_"$2"_x_NA"}')

      if [ $number_TMHMM -eq 0 ]
      then
        type=Secreted
      else
        type=Membrane
      fi
    fi
    echo $prot"_x_"$type"_x_"$program"_x_"$res_signal"_x_"$number_TMHMM"_x_"$interprot_data | sed "s/_x_/\t/g" | sed "s/_ooo_/ /g" | sed "s/;/; /g" >> $work_dir"/Final_Anotation_Files/Paper_Special_Tables/Location_Signals_paper.tab"
  done
fi

################################################################################
################################################################################
################################################################################

if [ "$DATA_CROSS_POLYA_REPEATS" == "TRUE" ]
then
  awk -F "\t" '{if ($3=="gene") print $1"\t"$4"\t"$5"\t"$9}' $work_dir"/TransDecoder/Bsudanica_annotation_polyA.gff3" | sed "s/;Name=.*//" | sed "s/ID=//" > $work_dir"/Data_Crossroads/Genes_coords.tab"
  awk -F "\t" '{if ($3=="five_prime_UTR") print $1"\t"$4"\t"$5"\t"$9}' $work_dir"/TransDecoder/Bsudanica_annotation_polyA.gff3" | sed "s/;Name=.*//" | sed "s/ID=//" > $work_dir"/Data_Crossroads/five_prime_UTR_coords.tab"
  awk -F "\t" '{if ($3=="CDS") print $1"\t"$4"\t"$5"\t"$9}' $work_dir"/TransDecoder/Bsudanica_annotation_polyA.gff3" | sed "s/;Name=.*//" | sed "s/ID=//" > $work_dir"/Data_Crossroads/CDS_coords.tab"
  awk -F "\t" '{if ($3=="three_prime_UTR") print $1"\t"$4"\t"$5"\t"$9}' $work_dir"/TransDecoder/Bsudanica_annotation_polyA.gff3" | sed "s/;Name=.*//" | sed "s/ID=//" > $work_dir"/Data_Crossroads/three_prime_UTR_coords.tab"

  echo "Gen Gen_Chr Gen_Start Gen_End N_Element Repeat_element Rep_Chr Rep_Start Rep_End Features_Overlap" | tr " " "\t" > $work_dir"/Data_Crossroads/Genes_overlapping_repeats.tab"
  count=1
  total_genes=$(grep -c . $work_dir"/Data_Crossroads/Genes_coords.tab")
  while [ $count -le $total_genes ]
  do
    chr=$(sed -n $count"p" $work_dir"/Data_Crossroads/Genes_coords.tab" | awk -F "\t" '{print $1}')
    start=$(sed -n $count"p" $work_dir"/Data_Crossroads/Genes_coords.tab" | awk -F "\t" '{print $2}')
    end=$(sed -n $count"p" $work_dir"/Data_Crossroads/Genes_coords.tab" | awk -F "\t" '{print $3}')
    genid=$(sed -n $count"p" $work_dir"/Data_Crossroads/Genes_coords.tab" | awk -F "\t" '{print $4}')

    echo $genid" "$chr": "$start".."$end
    awk -F "\t" -v chr=$chr -v start=$start -v end=$end '{if ($1==chr && $4<=end &&  start<=$5) print}' $work_dir"/BsudanicaEarlGrey/Bsudanica_summaryFiles/Bsudanica.filteredRepeats.gff" > $work_dir"/Data_Crossroads/Repeat_info.tab"
    check_data=$(grep -c . $work_dir"/Data_Crossroads/Repeat_info.tab")

    if [ $check_data -gt 0 ]
    then
      awk -F "\t" -v genID=$genid '{print $3" "$1" "$4" "$5}' $work_dir"/Data_Crossroads/Repeat_info.tab" > $work_dir"/Data_Crossroads/Genes_overlapping_repeats.temp"
      extract=1

      while [ $extract -le $check_data ]
      do
        current_element=$(echo $extract"_of_"$check_data)
        Rep_element=$(sed -n $extract"p" $work_dir"/Data_Crossroads/Genes_overlapping_repeats.temp" | awk -F " " '{print $1}')
        Rep_chr=$(sed -n $extract"p" $work_dir"/Data_Crossroads/Genes_overlapping_repeats.temp" | awk -F " " '{print $2}')
        Rep_start=$(sed -n $extract"p" $work_dir"/Data_Crossroads/Genes_overlapping_repeats.temp" | awk -F " " '{print $3}')
        Rep_end=$(sed -n $extract"p" $work_dir"/Data_Crossroads/Genes_overlapping_repeats.temp" | awk -F " " '{print $4}')

        echo "Debug:"$Rep_chr"||"$Rep_start"||"$Rep_end

        check_5UTR=$(grep -F $genid"." $work_dir"/Data_Crossroads/five_prime_UTR_coords.tab" | awk -F "\t" -v chr=$Rep_chr -v start=$Rep_start -v end=$Rep_end '{if ($1==chr && $2<=end &&  start<=$3) print}' | grep -c . | awk '{if ($1>=1) print "5UTR"; else print "NA"}')
        check_CDS=$(grep -F $genid"." $work_dir"/Data_Crossroads/CDS_coords.tab" | awk -F "\t" -v chr=$Rep_chr -v start=$Rep_start -v end=$Rep_end '{if ($1==chr && $2<=end &&  start<=$3) print}' | grep -c . | awk '{if ($1>=1) print "CDS"; else print "NA"}')
        check_3UTR=$(grep -F $genid"." $work_dir"/Data_Crossroads/three_prime_UTR_coords.tab" | awk -F "\t" -v chr=$Rep_chr -v start=$Rep_start -v end=$Rep_end '{if ($1==chr && $2<=end &&  start<=$3) print}' | grep -c . | awk '{if ($1>=1) print "3UTR"; else print "NA"}')

        echo "5UTR:"$check_5UTR" CDS:"$check_CDS" 3UTR:"$check_3UTR

        if [ "$check_5UTR" == "NA" ] && [ "$check_CDS" == "NA" ] && [ "$check_3UTR" == "NA" ]
        then
          Rep_Overlap=$(echo "Intron")
        else
          Rep_Overlap=$(echo $check_5UTR" "$check_CDS" "$check_3UTR | tr " " "\n" | grep -w -v "NA" | tr "\n" ";")
        fi

        echo $genid" "$chr" "$start" "$end" "$current_element" "$Rep_element" "$Rep_chr" "$Rep_start" "$Rep_end" "$Rep_Overlap | tr " " "\t" >> $work_dir"/Data_Crossroads/Genes_overlapping_repeats.tab"
        extract=$(($extract + 1))
      done
    fi
    count=$(($count + 1))
  done
  rm $work_dir"/Data_Crossroads/Genes_coords.tab"
  rm $work_dir"/Data_Crossroads/Repeat_info.tab"
  rm $work_dir"/Data_Crossroads/Genes_overlapping_repeats.temp"
  rm $work_dir"/Data_Crossroads/five_prime_UTR_coords.tab"
  rm $work_dir"/Data_Crossroads/CDS_coords.tab"
fi

################################################################################
################################################################################
################################################################################

if [ "$GFF" == "TRUE" ]
then
  # All the GFFs are on different formats. This section puts them all toguether and homogenises some key fields

  # All RNA genes:$work_dir/Mixed_Anot_Stringtie/PolyA_Transcripts_anot.gff
  # stringtie --mix --conservative -l BSUD -v -p 10 -o /home/amanda/PROYECTS/Project_Bsudanica/Bsudanica_Runs/Final_Genome_Annotation/Mixed_Anot_Stringtie/PolyA_Transcripts_anot.gff /home/amanda/PROYECTS/Project_Bsudanica/Bsudanica_Runs/Final_Genome_Annotation/STAR_MAP/All_work_illumina.bam /home/amanda/PROYECTS/Project_Bsudanica/Bsudanica_Runs/Final_Genome_Annotation/PACBIO_MAP/Minimap2_Pacbio_aln.bam
  # StringTie version 2.2.1
  # contig_1000     StringTie       transcript      10340   12894   1000    -       .       gene_id "BSUD.1"; transcript_id "BSUD.1.1"; cov "16.095701"; FPKM "0.998120"; TPM "2.803999";
  # contig_1000     StringTie       exon    10340   12327   1000    -       .       gene_id "BSUD.1"; transcript_id "BSUD.1.1"; exon_number "1"; cov "16.260857";
  # contig_1000     StringTie       exon    12798   12894   1000    -       .       gene_id "BSUD.1"; transcript_id "BSUD.1.1"; exon_number "2"; cov "10.405499";

  # Annotation Info:
    # cov -> Exclude
    # exon_number -> Change for ID ID=BSUD.9999.5.p1.exon4
    # FPKM -> Exclude
    # gene_id -> Change for ID
    # TPM -> Exclude
    # transcript_id -> Change for Parent

  # Other changes: 1) Duplicate the transcriot line with the "gene" feature and without the Parent annotation

  # Coding genes: $work_dir"/TransDecoder/Bsudanica_annotation_polyA.gff3"
  # contig_1000     transdecoder    gene    10340   12894   .       -       .       ID=BSUD.1;Name=ORF%20type%3A5prime_partial%20len%3A146%20%28%2B%29%2Cscore%3D16.72
  # contig_1000     transdecoder    mRNA    10340   12894   .       -       .       ID=BSUD.1.1.p1;Parent=BSUD.1;Name=ORF%20type%3A5prime_partial%20len%3A146%20%28%2B%29%2Cscore%3D16.72
  # contig_1000     transdecoder    exon    12798   12894   .       -       .       ID=BSUD.1.1.p1.exon1;Parent=BSUD.1.1.p1

  # Annotation Info:
    # ID -> Keep
    # Name -> Replace for Interpro hits.
    # Parent-> Keep

  # Ribosome: $work_dir"/Ribosome/Ribosomes.gff"
  ##gff-version 3
  # contig_3817     barrnap:0.9     rRNA    23316   23468   1.9e-24 -       .       Name=5_8S_rRNA;product=5.8S ribosomal RNA
  # contig_3818     barrnap:0.9     rRNA    36968   37120   1.9e-24 +       .       Name=5_8S_rRNA;product=5.8S ribosomal RNA
  # contig_3819     barrnap:0.9     rRNA    2384    2535    4.4e-24 -       .       Name=5_8S_rRNA;product=5.8S ribosomal RNA

  # Annotation Info:
    # Name -> Changed to ID
    # note  -> Keep
    # product -> Change to Name

  # Annotation Features:
    # CDS
    # exon
    # five_prime_UTR
    # gene -> Annotation info
    # mRNA -> Annotation info
    # three_prime_UTR


  # tRNA: $work_dir"/tRNA_data/trna_anot.gff"
  ##gff-version 3
  # contig_1003     tRNAscan-SE     pseudogene      40272   40343   51.5    +       .       ID=contig_1003.trna1;Name=contig_1003.tRNA1-GlyTCC;isotype=Gly;anticodon=TCC;gene_biotype=pseudogene;
  # contig_1003     tRNAscan-SE     exon    40272   40343   .       +       .       ID=contig_1003.trna1.exon1;Parent=contig_1003.trna1;
  # contig_1003     tRNAscan-SE     pseudogene      71607   71675   25.7    -       .       ID=contig_1003.trna5;Name=contig_1003.tRNA5-G

  # Annotation Info:
    # anticodon -> Keep
    # gene_biotype -> Keep
    # ID -> Keep
    # isotype -> Keep
    # Name -> Keep
    # Parent -> Keep

  rm -r $work_dir"/Final_Anotation_Files/Full_GFF"
  mkdir $work_dir"/Final_Anotation_Files/Full_GFF"
  mkdir $work_dir"/Final_Anotation_Files/Full_GFF/Intermediary"
  mkdir $work_dir"/Final_Anotation_Files/Full_GFF/Temp"

  echo Wololooo >> $work_dir"/Final_Anotation_Files/Full_GFF/Temp/Already_done.temp"

  sed "s/\tpseudogene\t/\tpseudogene_tRNA\t/" $work_dir"/tRNA_data/trna_anot.gff" | sed "s/;$//" >> $work_dir"/Final_Anotation_Files/Full_GFF/Intermediary/trna_anot.gff"
  sed "s/Name=/ID=/g" $work_dir"/Ribosome/Ribosomes.gff" | sed "s/product=/Name=/g" >> $work_dir"/Final_Anotation_Files/Full_GFF/Intermediary/Ribosomes.gff"

  chromosomes=$(seqkit seq -n -i $genome_file)

  for chr in $chromosomes
  do
    echo $chr
    if [ "$chr" == "$exclude_mitochondria" ]
    then
      lenght=$(seqkit grep -p $chr $genome_file | seqkit fx2tab -n -i -l | awk -F "\t" '{print $2}')
      echo $chr" manual gene 1 "$lenght" . + . ID=Mitochondrial_Genome" | tr " " "\t" >> $work_dir"/Final_Anotation_Files/Full_GFF/Intermediary/mitochondria.gff"
      echo $chr" manual mRNA 1 "$lenght" . + . ID=Mitochondrial_Genome_Polycistron;Parent=Mitochondrial_Genome" | tr " " "\t" >> $work_dir"/Final_Anotation_Files/Full_GFF/Intermediary/mitochondria.gff"
    else
      check_protein=$(grep -w -F -c $chr $work_dir"/TransDecoder/Bsudanica_annotation_polyA.gff3")
      check_rna_only=$(grep -w -F -c $chr $work_dir"/Mixed_Anot_Stringtie/PolyA_Transcripts_anot.gff")

      gene_registry=AVOID
      if [ $check_protein -gt 0 ]
      then
        grep -w -F $chr $work_dir"/TransDecoder/Bsudanica_annotation_polyA.gff3" > $work_dir"/Final_Anotation_Files/Full_GFF/Temp/current_chr_prot.gff"

        count=1
        recorrer=$(grep -c . $work_dir"/Final_Anotation_Files/Full_GFF/Temp/current_chr_prot.gff")

        while [ $count -le $recorrer ]
        do
          echo "Debug "$chr" Transdecoder: "$count"/"$recorrer
          feature=$(sed -n $count"p" $work_dir"/Final_Anotation_Files/Full_GFF/Temp/current_chr_prot.gff" | awk -F "\t" '{print $3}')
          if [ $feature == "gene" ]
          then
            information=$(sed -n $count"p" $work_dir"/Final_Anotation_Files/Full_GFF/Temp/current_chr_prot.gff" | awk -F "\t" '{print $1" "$2" "$3" "$4" "$5" "$6" "$7" "$8}')
            ID_Part=$(sed -n $count"p" $work_dir"/Final_Anotation_Files/Full_GFF/Temp/current_chr_prot.gff" | awk -F "\t" '{print $9}' | tr ";" "\n" | grep "ID=")
            gene_id=$(echo $ID_Part | awk -F "=" '{print $2"."}')
            gene_registry=$(echo $gene_registry" "$gene_id | sed "s/\.$//")

            final_gff_interpro

            echo $information" "$ID_Part";"$interprot_anotation | tr " " "\t" | sed "s/_ooo_/ /g" >> $work_dir"/Final_Anotation_Files/Full_GFF/Intermediary/Bsudanica_Prot.gff"
          elif [ $feature == "mRNA" ]
          then
            information=$(sed -n $count"p" $work_dir"/Final_Anotation_Files/Full_GFF/Temp/current_chr_prot.gff" | awk -F "\t" '{print $1" "$2" "$3" "$4" "$5" "$6" "$7" "$8}' | sed "s/\tmRNA\t/\ttranscript\t/")
            ID_Part=$(sed -n $count"p" $work_dir"/Final_Anotation_Files/Full_GFF/Temp/current_chr_prot.gff" | awk -F "\t" '{print $9}' | tr ";" "\n" | grep "ID=")
            Parent_Part=$(sed -n $count"p" $work_dir"/Final_Anotation_Files/Full_GFF/Temp/current_chr_prot.gff" | awk -F "\t" '{print $9}' | tr ";" "\n" | grep "Parent=")
            gene_id=$(echo $Parent_Part | awk -F "=" '{print $2"."}')

            final_gff_interpro
            echo $information" "$ID_Part";"$Parent_Part";"$interprot_anotation | tr " " "\t" | sed "s/_ooo_/ /g" >> $work_dir"/Final_Anotation_Files/Full_GFF/Intermediary/Bsudanica_Prot.gff"
          else
            sed -n $count"p" $work_dir"/Final_Anotation_Files/Full_GFF/Temp/current_chr_prot.gff" >> $work_dir"/Final_Anotation_Files/Full_GFF/Intermediary/Bsudanica_Prot.gff"
          fi
          count=$(($count + 1))
        done
      fi

      if [ $check_rna_only -gt 0 ]
      then
        grep -w -F $chr $work_dir"/Mixed_Anot_Stringtie/PolyA_Transcripts_anot.gff" > $work_dir"/Final_Anotation_Files/Full_GFF/Temp/current_chr_prot.gff"

        count=1
        recorrer=$(grep -c . $work_dir"/Final_Anotation_Files/Full_GFF/Temp/current_chr_prot.gff")

        while [ $count -le $recorrer ]
        do
          echo "Debug "$chr" RNA Only:"$count"/"$recorrer
          confirm_gene=$(sed -n $count"p" $work_dir"/Final_Anotation_Files/Full_GFF/Temp/current_chr_prot.gff" | awk -F "\t" '{print $9}' | tr ";" "\n" | sed "s/^ //" | grep gene_id | awk -F " " '{print $2}' | tr -d '"')
          check_skipped_gene=$(echo $gene_registry | tr " " "\n" | grep -w -F -c $confirm_gene )

          if [ $check_skipped_gene -eq 0 ]
          then
            information1=$(sed -n $count"p" $work_dir"/Final_Anotation_Files/Full_GFF/Temp/current_chr_prot.gff" | awk -F "\t" '{print $1" "$2}')
            feature=$(sed -n $count"p" $work_dir"/Final_Anotation_Files/Full_GFF/Temp/current_chr_prot.gff" | awk -F "\t" '{print $3}')
            information2=$(sed -n $count"p" $work_dir"/Final_Anotation_Files/Full_GFF/Temp/current_chr_prot.gff" | awk -F "\t" '{print $4" "$5" "$6" "$7" "$8}')
            pre_annotation=$(sed -n $count"p" $work_dir"/Final_Anotation_Files/Full_GFF/Temp/current_chr_prot.gff" | awk -F "\t" '{print $9}')

            GenID=$(echo $pre_annotation | tr ";" "\n" | sed "s/^ //" | grep gene_id | awk -F " " '{print $2}' | tr -d '"')
            TranscriptID=$(echo $pre_annotation | tr ";" "\n" | sed "s/^ //" | grep transcript_id | awk -F " " '{print $2}' | tr -d '"')

            if [ $feature == "transcript" ]
            then

              check_complete=$(grep -w -F -c $GenID $work_dir"/Final_Anotation_Files/Full_GFF/Temp/Already_done.temp")
              if [ $check_complete -eq 0 ]
              then
                echo $information1" gene "$information2" ID="$GenID | tr " " "\t" >> $work_dir"/Final_Anotation_Files/Full_GFF/Intermediary/Bsudanica_RNA_only.gff"
                echo $GenID >> $work_dir"/Final_Anotation_Files/Full_GFF/Temp/Already_done.temp"
              fi
              echo $information1" "$feature" "$information2" ID="$TranscriptID";Parent="$GenID | tr " " "\t" >> $work_dir"/Final_Anotation_Files/Full_GFF/Intermediary/Bsudanica_RNA_only.gff"
            elif [ $feature == "exon" ]
            then
              exon_number=$(echo $pre_annotation | tr ";" "\n" | sed "s/^ //" | grep exon_number | awk -F " " '{print $2}' | tr -d '"')
              echo $information1" "$feature" "$information2" ID="$TranscriptID".exon"$exon_number";Parent="$TranscriptID  | tr " " "\t" >> $work_dir"/Final_Anotation_Files/Full_GFF/Intermediary/Bsudanica_RNA_only.gff"
            fi

            # Other changes: 1) Duplicate the transcriot line with the "gene" feature and without the Parent annotation
            # Annotation Info:
              # cov -> Exclude
              # exon_number -> Change for ID ID=BSUD.9999.5.p1.exon4
              # FPKM -> Exclude
              # gene_id -> Change for ID
              # TPM -> Exclude
              # transcript_id -> Change for Parent
          fi

          count=$(($count + 1))
        done
      fi
    fi
  done
  cat $work_dir"/Final_Anotation_Files/Full_GFF/Intermediary/"*".gff" | grep -v "#" >> $work_dir"/Final_Anotation_Files/Full_GFF/Intermediary/Concatenated.gff"
fi

################################################################################
################################################################################
################################################################################

if [ "$GFF_SORT" == "TRUE" ]
then
  rm $work_dir"/Final_Anotation_Files/Full_GFF/Bsudanica_Anotation_Full.gff"
  rm -r $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort"
  mkdir $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort"

  grep -v "Parent=" $work_dir"/Final_Anotation_Files/Full_GFF/Intermediary/Concatenated.gff" >> $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/features.tmp"

  chromosomes=$(awk -F "\t" '{print $1}' $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/features.tmp" | grep contig | sort -u | tr "_" "\t" | sort -n -k 2 | tr "\t" "_")
  sorting_gff

  chromosomes=$(awk -F "\t" '{print $1}' $work_dir"/Final_Anotation_Files/Full_GFF/Temp_sort/features.tmp" | grep scaffold | sort -u | tr "_" "\t" | sort -n -k 2 | tr "\t" "_")
  sorting_gff

fi
