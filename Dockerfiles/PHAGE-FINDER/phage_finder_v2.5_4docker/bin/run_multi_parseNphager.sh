#!/bin/sh

# Bash script to parse .gb files from the Downloader and run Phage_Finder in parallel on a local Linux machine

# Usage: run_multi_parseNphager.sh -l <listfile> -b <destination base directory>

NO_ARGS=0
OPTERROR=65

if [ $# -eq "$NO_ARGS" ] # Script invoked with no command-line args?
  then
    echo "Usage: `basename $0` options (-hbl)"
    exit $OPTERROR      # Exit and explain usage, if no argument(s) given.
fi

while getopts ":hb:l:" Option
do
  case $Option in
    l) listfile=$OPTARG;;
    b) base=$OPTARG;;
    h) echo "Usage: `basename $0` options (-b [path of base dir] -l [filename containing GenBank accessions])"
       exit 0;;
  esac
done
shift $(($OPTIND - 1))  # Decrements the argument pointer so it points to next argument.
extentions=( att ent con pep cds )
concatarray=( pep con ent cds RNA )
tmp_dir="tmp_gbff_dir"
src="/scratch/dfouts/bin"

# store pids of subprocesses in a bash array
pids=()

let jobcounter=0
let totalcnter=0
let batchnumber=0
let batchsize=50

task(){
  cntr=$1
  i=$2
  #echo "  ${cntr}: ${i}"
  #export ei=${i}
  if [ ! -d ${base}/${i} ] # if directory NOT present, make it
  then
    echo "  ${base}/${i} is not present, generating ..."
    mkdir ${base}/${i}
    mv ${i}.gb ${base}/${i}/
  else
    echo " ${base}/${i} is already present, skipping directory making step ..."

    # check if tmp_dir is present and delete
    if [ -d ${base}/${i}/${tmp_dir} ]
    then
      echo "  removing existing ${tmp_dir} ..."
      rm -rf ${base}/${i}/${tmp_dir}
    fi

    # check if concatarray files are present and delete
    echo "  perging any existing parsed files ..."
    if [ -s ${base}/${i}/phage_finder_info.txt ]
    then
      rm ${base}/${i}/phage_finder_info.txt
    fi
    for c in ${extentions[@]}  
    do
      find ${base}/${i}/ -name "*${c}" -exec rm {} \;
    done
  fi
  cd ${base}/${i}
  gb_file="${i}.gb"
  num_CDSs=`grep -c "^\s*CDS" ${gb_file}`
  if [ "${num_CDSs}" -gt "0" ] # if CDSs, then proceed, otherwise record which genome has no annotation
  then
    num_asmbls=`grep -c "^LOCUS" ${gb_file}`
    if [ "$num_asmbls" -gt "1" ] # if more than one assembly, then split it, else just clean up the names
    then
      mkdir ${tmp_dir}
      cp ${gb_file} ${tmp_dir}/
      cd ${tmp_dir}

      # Split multi-Genbank file ($gb_file)
      echo "Splitting $num_asmbls assemblies from ${gb_file}..."
      ${src}/split_multi_GenbankFile.pl -f ${gb_file} -b .

      # Parse each Genbank file ($gb_file)
      ${src}/GenbankGenomeParser.pl -accession_file accessions.lst > /dev/null 2>&1

      # Concatenate .pep, .con, .ent and .cds (.seq) files
      echo "Concatenating .pep, .con, .ent, and .cds files..."
      for c in ${concatarray[@]}
      do
        #echo " working on ${c} files ..."
        find ${base}/${i}/${tmp_dir} -name "*${c}" -exec cat {} >> ${base}/${i}/${tmp_dir}/${i}.${c}.tmp \; -exec rm {} \;
      done

     ####### Generate a single .att file to describe the entire genome rather than each contig #######
      echo "Creating a single genome.att file..."
      cat `find -name "*.att"` | tail -1 | perl -ne 'chomp; print "$_";' > ${i}.att.tmp

      # make a unordered pseudomolecule to determine total genome size, gc content and number of protein sequences for use in the combined .att file
      echo ">${i} unordered pseudomolecule" > ${i}_pseudomolecule.tmp 
      cat ${i}.con.tmp | sed '/^>.*$/d' >> ${i}_pseudomolecule.tmp # remove all contig header lines

      ${src}/cleanFasta ${i}_pseudomolecule.tmp > /dev/null 2>&1  # clean up the fasta format
      rm ${i}_pseudomolecule.tmp_orig

      ${src}/getlengths ${i}_pseudomolecule.tmp | perl -ne 'chomp; @a=split(/\s+/); print "\t$a[1]";' >> ${i}.att.tmp # append total size of the unordered pseudomolecule

      ${src}/gc.pl ${i}_pseudomolecule.tmp | perl -ne 'chomp; @a=split(/\s+/); print "\t$a[1]";' >> ${i}.att.tmp # append gc content of the unordered pseudomolecule

      grep -c "^>" ${i}.pep.tmp | perl -ne 'chomp; print "\t$_\t";' >> ${i}.att.tmp # append the total number of proteins to the att temp file
      echo "${i}" >> ${i}.att.tmp # add the parent accession number for genome to end

      # Make a single clean genome.att file
      cat ${i}.att.tmp | perl -ne 'chomp; @a=split(/\t/); $a[10] =~ s/%//; print "$a[14]\t$a[11]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[13]\t$a[7]\t$a[8]\t$a[9]\t$a[10]\n";' > ${base}/${i}/${i}.att  # generate new genome att file     
     #######-----------------------------------------------------------------------------------#######
      # Rename and copy over concatenated .tmp files one directory up
      for c in ${concatarray[@]}
      do
        cp ${i}.${c}.tmp ${base}/${i}/${i}.${c};
      done

      # return to genome base directory
      cd ${base}/${i}

   ############################## F O R  S I N G L E  A S S E M B L Y .gb  F I L E  ####################################
    else
        echo "$i only has a single assembly, so just parsing and generating phage_finder files..."

        # Generate accessions.lst file
        echo "${gb_file%.*}" > accessions.lst

        # Parse Genbank files
        echo "parsing ${gb_file} ..."
        ${src}/GenbankGenomeParser.pl -accession_file accessions.lst -keep > /dev/null 2>&1

        # rename the .con, .pep, .cds and .ent files to match directory name and .gb filename
        for c in ${extentions[@]}
        do
          #echo " working on $c files ..."
          mv ${base}/${i}/*.${c} ${base}/${i}/${i}.${c};
        done

        # process .att file
        cat ${i}.att | perl -ne 'chomp; print "$_\t";' > ${i}.att.tmp # get rid of carriage return and replace with trailing tab
        echo "${i}" >> ${i}.att.tmp # add the parent accession number for genome to end
        cat ${i}.att.tmp | perl -ne 'chomp; @a=split(/\t/); print "$a[11]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[9]\t$a[10]\n";' > ${i}.att  # generate new genome att file
      
      # remove tmp files
        rm *.tmp
    fi
##############################-----------------------------------------------------####################################
  else
    echo ${i} >> ${base}/no_annotation.lst
  fi

  # Generate phage_finder_info.txt file
  echo "Generating the phage_finder_info.txt file..."
  ${src}/convert_new_ent_to_phageFinder_v2.pl ${i}.ent

  # add genome.att file to master.att
  cat ${base}/${i}/${i}.att >> ${base}/master.att

  # remove the temp dir
  rm -rf ${tmp_dir}

  # remove any RNAs (I do my own predictions)
    if [ -e `ls *RNA | head -1` ]
    then
      rm *RNA
    fi

  #run phage_finder wrapper (can comment out if not running phage_finder)
  #echo "initiating Phage_Finder or ${i} ..."
  #${src}/phage_finder_v2.5/bin/parallel_phage_finder.sh ${i}
  #$src/phage_finder_v2.5/bin/phage_finder_v3.sh $i

  cd $base
}

#### M A I N  L O O P ####

for genome in `cat ${listfile}`
do
    let jobcounter=${jobcounter}+1
    let totalcnter=${totalcnter}+1

    ## Spawn multiple jobs
    
    if [ "${jobcounter}" -eq "1" ]
    then
      let tmpbatchnbr=${batchnumber}+1
      echo "${totalcnter}) Launching batch ${tmpbatchnbr}"
    fi
    #echo -n "${jobcounter}: "
    echo "  ${jobcounter}: ${genome}"

    # run the task as a subprocess
    task "${jobcounter}" "${genome}" > /dev/null 2>&1 &

    # store the subprocess id
    pids+=("$!")

    # check if we are at the batchsize
    if [ "${jobcounter}" -eq "${batchsize}" ]
    then
      let batchnumber=${batchnumber}+1
      sleep 0.5s
      echo "Waiting for batch ${batchnumber} jobs to complete"
      
      # obtained this "OPTION 2" code from Gabriel Staples' post on stackoverflow.com
      # source url:  https://stackoverflow.com/questions/356100/how-to-wait-in-bash-for-several-subprocesses-to-finish-and-return-exit-code-0
      # when each process terminates, and print out when each process finishes!
      while true
      do
        for i in "${!pids[@]}"
        do
          pid="${pids[$i]}"
          # echo "pid = $pid"  # debugging

          # See if PID is still running; see my answer here:
          # https://stackoverflow.com/a/71134379/4561887
          ps --pid "$pid" > /dev/null
          if [ "$?" -ne 0 ]
          then
            # PID doesn't exist anymore, meaning it terminated
            # 1st, read its return code

            wait "$pid"
            return_code="$?"

            # 2nd, remove this PID from the `pids` array by `unset`ting the
            # element at this index; NB: due to how bash arrays work, this does
            # NOT actually remove this element from the array. Rather, it
            # removes its index from the `"${!pids[@]}"` list of indices,
            # adjusts the array count(`"${#pids[@]}"`) accordingly, and it sets
            # the value at this index to either a null value of some sort, or
            # an empty string (I'm not exactly sure).

            unset "pids[$i]"
            num_pids="${#pids[@]}"
            echo "  PID $pid is done; return_code = $return_code;" \
                 "$num_pids PIDs remaining."
          fi
        done

        # exit the while loop if the `pids` array is empty
        if [ "${#pids[@]}" -eq 0 ]
        then
          break
        fi

        # Do some small sleep here to keep your polling loop from sucking up
        # 100% of one of your CPUs unnecessarily. Sleeping allows other processes
        # to run during this time.
        sleep 0.1
      done
      echo
      echo "Batch ${batchnumber} jobs are complete"
      # reset jobcounter
      let jobcounter=0
    fi
done
let batchnumber=$batchnumber+1
echo -n "Waiting for batch $batchnumber jobs to complete"
wait


echo
echo "Batch $batchnumber jobs are complete"
echo "........................................"
echo "A total of $totalcnter jobs were completed."
#cd $basedir
echo "Running Post-Processing script for <strict> mode..."
#$src/parse_multi_phph_4complete_2mode.pl -a master.att -m strict
#echo "Finished processing $count datasets."
#echo "Finished searching $count genomes for Phages!" >> phage_phinder_postprocess_strict.out
#echo "Running Post-Processing script for <nonstrict> mode..."
#$src/parse_multi_phph_4complete_2mode.pl -a master.att -m nonstrict
#echo "Finished processing $count datasets."
#echo "Finished searching $count genomes for Phages!" >> phage_phinder_postprocess_nonstrict.out
counter=0
#for i in `cat $listfile`
#do
#    if [ "$counter" -eq "0" ] # printer the whole PFPR_tab.txt file to include the header for the first genome (i.e., we only want one header line in this file)
#    then
#	cat $i/strict_dir/PFPR_tab.txt >> combined_PFPR_strict_tab.txt
#    else
#        cat $i/strict_dir/PFPR_tab.txt | grep -v "^#" >> combined_PFPR_strict_tab.txt
#    fi
#    let counter=counter+1 # increment the genome counter
#done 
