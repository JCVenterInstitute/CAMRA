#!/bin/bash


if [ "$1" == "--help" ]; then
  echo"
  PURPOSE: Pulls the reads from a GoogleBucket and pushesthem to the BVBRC to be assembled.  
  TO DO BEFORE RUNNING SCRIPT: 
    1. Create a BV-BRC account.  https://www.bv-brc.org/
    2. Install the BV-BRC command line interface called PATRIC. https://www.bv-brc.org/docs/cli_tutorial/index.html
    3. Open the application, login with your BV-BRC credencials using the p3-login script.
    4. Make sure you have your raw-read directory and output-assembly directory made in BV-BRC. See "INPUT" section for more information.
    5. Make your working directory. 
    6. Know the location of the raw-reads in the google bucket.
    7. Assure exsisting location to deposit the assemblies in the google bucket.
    8. Run script Pull_many_Gfiles_push_to_BVBRC-assemb.sh
    9. After all assemblies are complete, run From_BVBRC-JSON_getFiles.py
  SCRIPT STEPS:
    - Makes a sample{datetime}.txt that contains all the locations of the files in the G-bucket
    - Finds the read pairs (SAMPLE123_r1.fq.gz and SAMPLE123_r2.fq.gz)
    - Pushes these pairs onto BRVRC for assembly. 
  USE: From_BVBRC-JSON_getFiles.py <BVBRC-ASSEMBLY LOCATION> <GBUCKET-ASEMBLY-LOCATION>
  INPUTS EXAMPLES:
    Google Bucket Raw Read Location: gs://fc-1bde9971-12a8-416d-a2ff-e314c7ade234/CAMRA-IFAIN/Retrospective/BatchA_121022/Raw_Read/
    BVBRC raw reads storage location: /DanyMB@bvbrc/home/IFAIN_Retro_BatchA/Raw_Read/
    BVBRC Assembly storage Location: /DanyMB@bvbrc/home/IFAIN_Retro_BatchA/Assemblies/
    G Bucket location where assemblies will be deposited. This will only be needed for the second script.
        Example gs://fc-1bde9971-12a8-416d-a2ff-e314c7ade234/CAMRA-IFAIN/WGS/Retrospective/CAMRA_Seq_BatchB_25003/"
  exit 0
fi

# User input, examples underneath
read -p "Google Bucket Raw Read Location: " gbucket 
read -p "BVBRC raw reads storage location: " raw_location
read -p "BVBRC Assembly storage Location: " bvbrc_location



current_datetime=$(date "+%d%h%y_%H%M")

# Check if there are files in the google bucket
if gsutil -q ls "$gbucket" > "samples${current_datetime}.txt"; then
  echo "Google bucket and files exists. GoogleBucket's file locations were stored in samples${current_datetime}.txt"
  
else
  echo "Google bucket and/or files do not exist."
  exit 1
fi



count_this=0
prev_line="" 

# Read the list of file locations from the input file and find pairs
while IFS= read -r current_line; do
  # Extract the file name by removing the path
  current_basename=$(basename "$current_line")

  # Extract the common prefix by removing the last portion of the file name
  # Use this as BVBRC file name
  current_prefix="${current_basename%_R*}"

  # Check if the current line has the same common prefix as the previous line
  if [ "$current_prefix" = "$prev_prefix" ]; then
    file_1=$(basename "$prev_line")
    file_2=$(basename "$current_line")

    echo "$file_1"
    echo "$file_2"

    temp_file1=$(mktemp "/tmp/$file_1")
    temp_file2=$(mktemp "/tmp/$file_2")

    gsutil cp $prev_line $temp_file1
    gsutil cp $current_line $temp_file2

    # submits a job into BVBRC
    p3-submit-genome-assembly --trim-reads --workspace-upload-path $raw_location --recipe unicycler --paired-end-lib $temp_file1 $temp_file2  --min-contig-len 300 --min-contig-cov 30 $bvbrc_location $current_prefix 
    count_this=$((count_this + 1))
  fi

  # Set the current prefix as the previous prefix for the next iteration
  prev_prefix="$current_prefix"
  prev_line="$current_line"
  
done < "samples${current_datetime}.txt"

echo "Your assembly jobs have been submitted to the BVBRC. Check their status on the BVBRC website. Run From_BVBRC-JSON_getFiles.py when all the submitted assemblies have finished."
