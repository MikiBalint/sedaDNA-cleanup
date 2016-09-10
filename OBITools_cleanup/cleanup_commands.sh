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
  #Primers (from Taberlet): i18S_V9_F (TCACAGACCTGTTATTGC 18bp), i18S_V9_R (TYTGTCTGSTTRATTSCG 18bp)
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

  #Remove primers from the read ends
  #Primers are 18 bp
  #Insert the sequence length (seq_length) into the headers
  ls *.fastq | parallel -j 60 'obiannotate --length --fasta-output {} > {.}_wlength.fasta'

  #Trim the reverse primers from the read ends. Reverse primer bases calculated with the seq_length
  ls *.fasta | parallel -j 60 'obicut -b 19 -e seq_length-18 --uppercase {} > {.}_noprimer.fasta'

# 4. Remove unaligned

    cd /phylodata/mbalint/workdir/Stechlin_preexp
    mkdir 04_unaligned
    cd 04_unaligned

  #the steps are different from the conventional OBITools pipeline from now on: the samples are already demultiplexed, but the sample names are not yet in the fastq header lines.

  #The sample names should be inserted as 'sample=xxx; ' before the 'ali_length=xxx;'
  #The rename.pl is a small script modified from **Bálint Ecol Evol**
  #Original in /phylodata/mbalint/scripts

    cp /phylodata/mbalint/scripts/rename_fasta.pl . # copy the rename script
    ln -s ../03_paired-end/*_noprimer.fasta .
    perl rename_fasta.pl
    rm *noprimer.fasta

AATGATACGGCGACCACCGAGATCTACACAAGTCGGATCGTCGGCAGCGTC
CAAGCAGAAGACGGCATACGAGATTTGGAGTGGTCTCGTGGGCTCGG0
  #Now combine the files

    cat *renamed.fasta > combined.fasta
    obigrep -p 'mode!="joined"' combined.fasta > combined_ali.fasta

  #Successfully assembled sequences
    grep -c "^>" combined_ali.fasta

  #Remove sequences shorter, that 50 bp
    obigrep -l 50 combined_ali.fasta > combined_noshort.fasta

    grep -c "^>" combined_noshort.fasta

# 5 chimera checking
    cd /phylodata/mbalint/workdir/Stechlin_preexp
    mkdir 05_chimera
    cd 05_chimera
    ln -s ../04_unaligned/combined_noshort.fasta .

#vsearch for chimera checking
##Documentation here:
##https://wiki.gacrc.uga.edu/wiki/Vsearch#Documentation

  vsearch --threads 60 --uchime_denovo combined_noshort.fasta --uchimeout chimera_results.tab --nonchimeras stechlin_nonchimera.fasta --chimeras stechlin_chimera.fasta

# 6 De-replicate into unique sequences

    cd /phylodata/mbalint/workdir/Stechlin_preexp
    mkdir 06_derep
    cd 06_derep
  #Since all were non-chimeric and the VSEARCH chimera checking screws the header, let's use the chimieracheck input file
    obiuniq -m sample ../05_chimera/combined_noshort.fasta > stechlin_derep.fasta
    grep -c "^>" stechlin_derep.fasta

# 7 Denoise the dataset

  #Get the counts of the 20 rarest sequences

    cd /phylodata/mbalint/workdir/Stechlin_preexp
    mkdir 07_denoise
    cd 07_denoise
    obistat -c count ../06_derep/stechlin_derep.fasta | sort -nk1 | head -20 > rare_counts.txt

#Keep only the sequence variants having a count greater or equal to 10:

    obigrep -p 'count>=10' ../06_derep/stechlin_derep.fasta > stechlin_c10.fasta
    grep -c "^>" stechlin_c10.fasta

# 8 Clean the sequences

    cd /phylodata/mbalint/workdir/Stechlin_preexp
    mkdir 08_clean
    cd 08_clean

  #keep only head sequences ( -H option) if these are sequences with no variants with a count greater than 5% of their own count ( -r 0.05 option). This also annotates sequences as head, internal or singleton in a sample.

    obiclean -s merged_sample -r 0.05 -H ../07_denoise/stechlin_c10.fasta > stechlin_clean.fasta
    grep -c "^>" stechlin_clean.fasta

# 9 ecoPCR databases

  ## Download and format the databases. These steps are commented out if the databases exist.
  #Actual database EMBL v128, get it

    cd /phylodata/mbalint/databases/
    mkdir EMBL_128
    cd EMBL_128
    wget ftp://ftp.ebi.ac.uk/pub/databases/embl/release/std/rel_std*.*

  #Get the actual GenBank taxonomy

    cd /phylodata/mbalint/databases/
    mkdir Taxonomy_160908
    cd Taxonomy_160908
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -zxvf taxdump.tar.gz

  #Create the EMBL v128 ecoPCR database. Use '--skip-on-error' to get through some problematic EMBL entries.
  #It may not be necessary to unpack 'gz' beforehand, try using '*.dat.gz'

    cd /phylodata/mbalint/databases
    mkdir ecoPCR_embl_128
    cd ecoPCR_embl_128/
    obiconvert --skip-on-error --embl -t ../Taxonomy_160908 --ecopcrdb-output=ecopcr_embl_128 ../EMBL_128/*.dat.gz

  #Generate an ecoPCR assignment database with the i18S_V9_F, i18S_V9_R primers

    cd /phylodata/mbalint/workdir/Stechlin_preexp
    mkdir 09_ecopcr
    cd 09_ecopcr
    ecoPCR -d /phylodata/mbalint/databases/ecoPCR_embl_128/ecopcr_embl_128 -e 3 -l 50 -L 500 TCACAGACCTGTTATTG#C# TYTGTCTGSTTRATTSC#G# > i18S_V9.ecoPCR

    grep -cv '#' *.ecoPCR

  #Clean the databases

  #- filter sequences so that they have a good taxonomic description at the species, genus, and family levels
  obigrep -d /phylodata/mbalint/databases/ecoPCR_embl_128/ecopcr_embl_128 --require-rank=species --require-rank=genus --require-rank=family i18S_V9.ecoPCR > i18S_V9_clean.fasta

  #remove redundant sequences (obiuniq command below).
  obiuniq -d /phylodata/mbalint/databases/ecoPCR_embl_128/ecopcr_embl_128 i18S_V9_clean.fasta > i18S_V9_clean_uniq.fasta

  #- ensure that the dereplicated sequences have a taxid at the family level (obigrep command below).
  obigrep -d /phylodata/mbalint/databases/ecoPCR_embl_128/ecopcr_embl_128 --require-rank=family i18S_V9_clean_uniq.fasta > i18S_V9_clean_uniq_clean.fasta

  #- ensure that sequences each have a unique identification (obiannotate command below)
  obiannotate --uniq-id i18S_V9_clean_uniq_clean.fasta > db_i18S_v9.fasta

# 10 Taxonomic assignment

  cd /phylodata/mbalint/workdir/Stechlin_preexp
  mkdir 10_assign
  cd 10_assign
  ecotag -d /phylodata/mbalint/databases/ecoPCR_embl_128/ecopcr_embl_128 -R ../09_ecopcr/db_i18S_v9.fasta --sort=count -m 0.98 -r ../08_clean/stechlin_clean.fasta > stechlin_assigned.fasta

  grep -c ">" stechlin_assigned.fasta

# 11_abundance_tables

  cd /phylodata/mbalint/workdir/Stechlin_preexp
  mkdir 11_abundance_tables
  cd 11_abundance_tables

  # create assigned abundance matrix
  obitab -o --output-field-separator=, ../10_assign/stechlin_assigned.fasta > stechlin_assigned.csv
