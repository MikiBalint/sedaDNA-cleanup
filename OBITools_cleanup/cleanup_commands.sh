# Trial lines some more

# General details
  #All done on the malloy
  #Data in: /phylodata/mbalint/workdir/Stechlin_preexp/data
  #Workdir: /phylodata/mbalint/workdir/Stechlin_preexp
  #OBITools v1.2.7, downloaded 160717

# 0: create a small dataset to check the script
#parallel citation: Tange 2011

    mkdir data
    cd data
    ln -s /phylodata/mbalint/datasets/160329_Stechlin-preexp/data/*.fastq .
    mkdir complete_files
    mv *.fastq complete_files/
    cd complete_files
    ls *.fastq | parallel -j 60 'head -n 400 {} > trial_{.}.fastq'
    mv trial*.fastq ../
    cd ../

# 1: Quality check:

    mkdir 01_fastqc
    ls *.fastq | parallel -j 60 'fastqc -o 01_fastqc'
    mv 01_fastqc /phylodata/mbalint/workdir/Stechlin_preexp
    cd /phylodata/mbalint/workdir/Stechlin_preexp/01_fastqc
    rm *.zip

# 2: Trimming
  #Check the trimming length with the primers.
  #Primers (from Taberlet): i18S_V9_F (TCACAGACCTGTTATTGC ), i18S_V9_R (TYTGTCTGSTTRATTSCG)
  #Fragment length: most species ~110 bp, some ~150 bp + primers 2 x 20 bp
  #All reads trimmed to 100 bp: this should be shorter then the fragment length, and allows the assembly of the ~150 bp reads (fragment + primer)

    cd /phylodata/mbalint/workdir/Stechlin_preexp
    mkdir 02_trimmed
    cd 02_trimmed
    ln -s ../data/*.fastq . # link because the parallel input is simpler for me this way
    ls *.fastq | parallel -j 60 'fastx_trimmer -t 51 -i {} -o {.}_trimmed.fastq' # -t: remove 51 nucleotides from the ends
    rm *001.fastq # remove the links. the -j 60 option of the parallel uses 60 processors

# 3: Paired-end assembly:

    cd /phylodata/mbalint/workdir/Stechlin_preexp
    mkdir 03_paired-end
    cd 03_paired-end

  #link all files into present folder, because it is less complicated to use the parallel inputting for me

    ln -s ../02_trimmed/*trimmed.fastq .

    ls *trimmed.fastq | grep -v R2_001_trimmed.fastq | sed 's/R1_001_trimmed.fastq//' | parallel -j 60 'illuminapairedend --score-min=40  -r {.}R2_001_trimmed.fastq {.}R1_001_trimmed.fastq > {.}paired.fastq'

  #1. list the filenames.
  #2. keep the file name stems by removing the forward-reverse specific file information with 'sed'.
  #3. give the filename stems to parallel as '{.}'.

  #remove the links

    rm *001_trimmed.fastq

  #Paired quality check

    mkdir fastqc_paired # create a folder for the fastqc output
    ls *_paired.fastq | parallel -j 60 'fastqc -o fastqc_paired {}'
    rm fastqc_paired/*.zip # keep only the htmls

# 4. Remove unaligned

    cd /phylodata/mbalint/workdir/Stechlin_preexp
    mkdir 04_unaligned
    cd 04_unaligned

  #the steps are different from the conventional OBITools pipeline from now on: the samples are already demultiplexed, but the sample names are not yet in the fastq header lines.

  #The sample names should be inserted as 'sample=xxx; ' before the 'ali_length=xxx;'
  #The rename.pl is a small script modified from **Bálint Ecol Evol**
  #Original in /phylodata/mbalint/scripts

    cp /phylodata/mbalint/scripts/rename.pl . # copy the rename script
    ln -s ../03_paired-end/*.fastq .
    perl rename.pl
    rm *paired.fastq

  #Now combine the files

    cat *renamed.fastq > combined.fastq
    obigrep -p 'mode!="joined"' combined.fastq > combined_ali.fastq
    grep -c "@M01271" combined_ali.fastq

# 5 De-replicate into unique sequences

    cd /phylodata/mbalint/workdir/Stechlin_preexp
    mkdir 05_derep
    cd 05_derep
    obiuniq -m sample ../04_unaligned/combined_ali.fastq > stechlin_derep.fasta
    grep -c "^>" stechlin_derep.fasta

# 6 Denoise the dataset

  #Get the counts of the 20 rarest sequences

    cd /phylodata/mbalint/workdir/Stechlin_preexp
    mkdir 06_denoise
    cd 06_denoise
    obistat -c count ../05_derep/stechlin_derep.fasta | sort -nk1 | head -20 > rare_counts.txt

#Keep only the sequence variants having a count greater or equal to 10:

    obigrep -p 'count>=10' ../05_derep/stechlin_derep.fasta > stechlin_c10.fasta
    grep -c "^>" stechlin_c10.fasta

# 7 Clean the sequences

    cd /phylodata/mbalint/workdir/Stechlin_preexp
    mkdir 07_clean
    cd 07_clean

    **#remove the primers from the ends of the sequences**

  #keep only head sequences ( -H option) if these are sequences with no variants with a count greater than 5% of their own count ( -r 0.05 option). This also annotates sequences as head, internal or singleton in a sample.

    obiclean -s merged_sample -r 0.05 -H ../06_denoise/stechlin_c10.fasta > stechlin_clean.fasta
    grep -c "^>" stechlin_clean.fasta

# 8 ecoPCR databases

  ## Download and format the databases. These steps are commented out if the databases exist.
  #Actual database EMBL v128, get it

    cd /phylodata/mbalint/databases/
    mkdir EMBL_128
    cd EMBL_128
    wget ftp://ftp.ebi.ac.uk/pub/databases/embl/release/std/rel_est*.*

  #Get the actual GenBank taxonomy

    cd /phylodata/mbalint/databases/
    mkdir Taxonomy_160717
    cd Taxonomy_160717
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -zxvf taxdump.tar.gz

  #Create the EMBL v128 ecoPCR database. Use '--skip-on-error' to get through some problematic EMBL entries.
  #It may not be necessary to unpack 'gz' beforehand, try using '*.dat.gz'

    cd /phylodata/mbalint/databases
    mkdir ecoPCR_embl_128
    cd ecoPCR_embl_128/
    obiconvert --skip-on-error --embl -t ../Taxonomy_160717 --ecopcrdb-output=ecopcr_embl_128 ../EMBL_128/*.dat

  #Generate an ecoPCR assignment database with the i18S_V9_F, i18S_V9_R primers

    cd /phylodata/mbalint/workdir/Stechlin_preexp
    mkdir 08_ecopcr
    cd 08_ecopcr
    ecoPCR -d /phylodata/mbalint/databases/ecoPCR_embl_128/ -e 3 -l 60 -L 500 TCACAGACCTGTTATTG#C# TYTGTCTGSTTRATTSC#G# > i18S_V9.ecoPCR


  grep -cv '#' */*.ecopcr
  16S/16S.ecopcr:189360

  Clean the databases

  - filter sequences so that they have a good taxonomic description at the species, genus, and family levels

  obigrep -d ../../../../../databases/ecoPCR_embl_125/ecoPCR_embl_125 --require-rank=species --require-rank=genus --require-rank=family 12S.ecopcr > 12S_clean.fasta
  obigrep -d ../../../../../databases/ecoPCR_embl_125/ecoPCR_embl_125 --require-rank=species --require-rank=genus --require-rank=family 16S.ecopcr > 16S_clean.fasta

  - remove redundant sequences (obiuniq command below).

  obiuniq -d ../../../../../databases/ecoPCR_embl_125/ecoPCR_embl_125 12S_clean.fasta > 12S_clean_uniq.fasta
  obiuniq -d ../../../../../databases/ecoPCR_embl_125/ecoPCR_embl_125 16S_clean.fasta > 16S_clean_uniq.fasta

  - ensure that the dereplicated sequences have a taxid at the family level (obigrep command below).

  obigrep -d ../../../../../databases/ecoPCR_embl_125/ecoPCR_embl_125 --require-rank=family 12S_clean_uniq.fasta > 12S_clean_uniq_clean.fasta
  obigrep -d ../../../../../databases/ecoPCR_embl_125/ecoPCR_embl_125 --require-rank=family 16S_clean_uniq.fasta > 16S_clean_uniq_clean.fasta

  - ensure that sequences each have a unique identification (obiannotate command below)

  obiannotate --uniq-id 12S_clean_uniq_clean.fasta > db_12S.fasta
  obiannotate --uniq-id 16S_clean_uniq_clean.fasta > db_16S.fasta




Assign so far not assigned sequences to EMBL :
10_assign
Own references
ecotag -d 16S_bolivian_frogs.db -R 16S_fasta_nonredundant_noprimers.fasta --sort=count -m 0.98 -r frogs_clean.fasta > frogs_own16S_assigned.fasta

# grep assigned
obigrep -v -a 'scientific_name:root' frogs_own16S_sort.fasta > frogs_16S_own_assigned.fasta

# grep not assigned
obigrep -a 'scientific_name:root' frogs_own16S_sort.fasta > frogs_16S_not_assigned.fasta

# counts

frogs_16S_assigned.fasta:4805
frogs_16S_not_assigned.fasta:9637

# assignment against the GenBank database
ecotag -d /phylodata/mbalint/databases/ecoPCR_embl_125/ecoPCR_embl_125 -R ../09_ecopcr/NCBI/16S/db_16S.fasta --sort=count -m 0.98 -r frogs_16S_not_assigned.fasta > frogs_16S_EMBL_assigned.fasta

# counts
grep -c "scientific_name=" *.fasta
16S_fasta_nonredundant_noprimers.fasta:86
frogs_16S_EMBL_assigned.fasta:9637
frogs_16S_not_assigned.fasta:9637
frogs_16S_own_assigned.fasta:4805
frogs_clean.fasta:0
frogs_own16S_annot.fasta:14442
frogs_own16S_assigned.fasta:14442
frogs_own16S_sort.fasta:14442

grep -c "scientific_name=root" *.fasta
16S_fasta_nonredundant_noprimers.fasta:0
frogs_16S_EMBL_assigned.fasta:7807
frogs_16S_not_assigned.fasta:9637
frogs_16S_own_assigned.fasta:0
frogs_clean.fasta:0
frogs_own16S_annot.fasta:9637
frogs_own16S_assigned.fasta:9637
frogs_own16S_sort.fasta:9637

11_abundance_tables
# Remove non-informative annotations
obiannotate  --delete-tag=scientific_name_by_db --delete-tag=obiclean_samplecount --delete-tag=obiclean_count --delete-tag=obiclean_singletoncount --delete-tag=obiclean_cluster --delete-tag=obiclean_internalcount --delete-tag=obiclean_head --delete-tag=taxid_by_db --delete-tag=obiclean_headcount --delete-tag=id_status --delete-tag=rank_by_db --delete-tag=order_name --delete-tag=order --delete_tag=sminR --delete_tag=sminL --delete_tag= frogs_16S_own_assigned.fasta > frogs_16S_own_annot.fasta

obiannotate  --delete-tag=scientific_name_by_db --delete-tag=obiclean_samplecount --delete-tag=obiclean_count --delete-tag=obiclean_singletoncount --delete-tag=obiclean_cluster --delete-tag=obiclean_internalcount --delete-tag=obiclean_head --delete-tag=taxid_by_db --delete-tag=obiclean_headcount --delete-tag=id_status --delete-tag=rank_by_db --delete-tag=order_name --delete-tag=order frogs_16S_EMBL_assigned.fasta > frogs_16S_EMBL_annot.fasta

# Sort them
obisort -k count -r frogs_16S_own_annot.fasta > frogs_16S_own_sort.fasta

obisort -k count -r frogs_16S_EMBL_annot.fasta > frogs_16S_EMBL_sort.fasta

# create 16S assigned abundance matrix
obitab -o --output-field-separator=, frogs_16S_own_sort.fasta > frogs_16S_own.tab

obitab -o --output-field-separator=, frogs_16S_EMBL_sort.fasta > frogs_16S_EMBL.tab
