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
let batchsize=10

task(){
  cntr=$1
  i=$2

  cd ${base}/${i}

  #run phage_finder wrapper (can comment out if not running phage_finder)
  #echo "initiating Phage_Finder or ${i} ..."
  ${src}/phage_finder_v2.5/bin/parallel_phage_finder.sh ${i}

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
