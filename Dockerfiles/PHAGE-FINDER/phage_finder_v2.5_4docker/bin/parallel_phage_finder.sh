#!/bin/sh

# Bash script to run the Phage_Finder pipeline in parallel on a single computing node

# Usage: parallel_phage_finder.sh <prefix of .pep/.ffa, .ptt and .con/.fna file>

# NOTE: a phage_finder_info.txt file will be searched before a .ptt file
# .pep is the multifasta protein sequence file
# .ptt is a GenBank .ptt file that has the coordinates and ORF names with annotation
# .con is a file that contains the complete nucleotide sequence of the genome being searched

# find where the bin directory is located
Bin="$(cd -- "$(dirname "$0")" && pwd)"

# set the phage home (phome) directory where all phage_finder scripts are located
phome=${Bin%/*}

# store the prefix of the genome
prefix=$1

# set the base directory
base=`pwd`

# set the path to the source (src) directory where tRNAScan, Aragorn and cmscan are installed
src="/scratch/dfouts/bin"

# set the path to your BLAST+ package
blastbin="${src}/ncbi-blast-2.14.0+/bin"

# store pids of subprocesses in a bash array
pids=()

# identify the name of the file containing protein sequences
  if [ -s ${base}/${prefix}.pep ] # check if .pep file is present
  then
      pepfile="${prefix}.pep"
  elif [ -s ${base}/${prefix}.faa ]
  then
      pepfile="${prefix}.faa"
  else
     echo "Could not find a ${prefix}.pep or ${prefix}.faa file.  Please check to make sure the file is present and contains data"
     exit 1
  fi

# identify the name of the file containing contig fasta sequences
  if [ -s ${base}/${prefix}.con ]
  then
     contigfile="${prefix}.con"
  elif [ -s ${base}/${prefix}.fna ]
  then
     contigfile="${prefix}.fna"
  else
     echo "Could not find a ${prefix}.con or ${prefix}.fna file.  Please make sure one of these files is present and contains data"
     contigfile=""
  fi
 
# Determine if the phage_finder_info.txt file or a GenBank .ptt file were available
  if [ -s ${base}/phage_finder_info.txt ] # check for phage_finder info file and if it has contents
  then
     infofile="phage_finder_info.txt"
  elif [ -s ${base}/${prefix}.ptt ]
  then
     infofile="${prefix}.ptt"
  else
     echo "Could not find a phage_finder_info.txt file or ${prefix}.ptt file.  Please make sure one of these files is present and contains data."
     exit 1
  fi

### HMM searches section ###

# if HMM results not present, then search, otherwise skip
  if [ ! -s ${base}/combined.hmm3 ] # if GLOCAL HMM results not present, search
  then
     ## conduct GLOCAL HMMER3 searches
     echo "  HMMER3 searches ..."
     ${phome}/bin/HMM3_searches.sh ${pepfile} > /dev/null 2>&1 &

     # store the subprocess id
     pids+=("$!")
   
     # let us know what is going on
     echo "  Running HMMER3 searches ..."
  fi

### NCBI searches section ###

# if BLAST results not present, then search, otherwise skip
  if [ ! -s ${base}/ncbi.out ]
  then
     ## do NCBI BLASTP searches
     echo "  BLASTing ${pepfile} against the Phage DB ..."

     # NCBI BLAST+ code provided by Dr. Camille Hankel at Dartmouth #
     ${blastbin}/blastp -db $phome/DB/phage_03_25_19.db -outfmt 6 -evalue 0.001 -query $pepfile -out ncbi.out -max_target_seqs 5 -num_threads 12 &
     
     # store the subprocess id
     pids+=("$!")
   
     # let us know what is going on
     echo "  Running BLASTP searches ..."
  fi

### RNA searches section ###

# if tRNAscan (tRNA) results not present, then search, othersie skip
  if [ ! -e ${base}/tRNAscan.out ] && [ ${base}/${contigfile} ] # if tRNAscan.out file not present, and contig file present, then search
  then
     ## find tRNAs
     echo "  find tRNA sequences ..."
     ${src}/tRNAscan-SE --thread 4 -B -o tRNAscan.out ${base}/${contigfile} > /dev/null 2>&1 &

     # store the subprocess id
     pids+=("$!")

     # let us know what is going on
     echo "  Running tRNAscan searches ..."
  fi

# if Aragorn (tmRNA) results not present, then	search,	othersie skip
  if [ ! -e ${base}/tmRNA_aragorn.out ] && [ -e ${base}/${contigfile} ] # if tRNAscan.out file not present, and contig file present, then search
  then
     ## find tmRNAs
     echo "  find tmRNA sequences ..."
     ${src}/aragorn -m -o tmRNA_aragorn.out ${base}/${contigfile} &

     # store the subprocess id
     pids+=("$!")

     # let us know what is going on
     echo "  Running Aragorn searches ..."
  fi

# if cmscan (ncRNA) results not present, then	search,	othersie skip
  if [ ! -e ${base}/ncRNA_cmscan.out ] && [ -e ${base}/${contigfile} ] # if ncRNA_cmscan.out file not present, and contig file present, then search
  then
     ## find select ncRNAs
     echo "  find ncRNA sequences ..."
     Z=`$src/esl-seqstat ${prefix}.con | perl -ne 'chomp; if (/^Total # residues:\s+(\d+)/) {$n = $1; $Z=($n*2)/1000000; print "$Z\n";}'`
     #echo "Z = $Z"        
     ${src}/cmscan --cpu 4 -Z $Z --cut_ga --rfam --nohmmonly --tblout ${base}/ncRNA_cmscan.out --fmt 2 ${phome}/RfamDB/Rfam_PhageFinder.cm ${prefix}.con > ${prefix}.cmscan &

     # store the subprocess id
     pids+=("$!")
   
     # let us know what is going on
     echo "  Running cmscan searches ..."
  fi

### Waiting for all subprocesses to finish ###

# obtained this "OPTION 2" code from Gabriel Staples' post on stackoverflow.com
# source url:  https://stackoverflow.com/questions/356100/how-to-wait-in-bash-for-several-subprocesses-to-finish-and-return-exit-code-0

# when each process terminates, and print out when each process finishes!
while true; do
    for i in "${!pids[@]}"; do
        pid="${pids[$i]}"
        # echo "pid = $pid"  # debugging

        # See if PID is still running; see my answer here:
        # https://stackoverflow.com/a/71134379/4561887
        ps --pid "$pid" > /dev/null
        if [ "$?" -ne 0 ]; then
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
    if [ "${#pids[@]}" -eq 0 ]; then
        break
    fi

    # Do some small sleep here to keep your polling loop from sucking up
    # 100% of one of your CPUs unnecessarily. Sleeping allows other processes
    # to run during this time.
    sleep 0.1
done

### Phage_finder and visualizations section ###

# Run phage_finder
  echo "  searching for Prophage regions ..."
  ${phome}/bin/Phage_Finder_v2.5.pl -t ncbi.out -i $infofile -r tRNAscan.out -n tmRNA_aragorn.out -c ncRNA_cmscan.out -A $contigfile -S > /dev/null 2>&1

# generate linear figures of each prophage region
  echo "  generating linear illustrations..."
  ${phome}/bin/LinearDisplay_detailed.pl -A ${base}/strict_dir/PFPR.frag -F ${base}/strict_dir/PFPR.att -L -nt > ${base}/strict_dir/PFPR_linear.fig

# check that PFPR_linear.fig exists and has size

  if [ -s ${base}/strict_dir/PFPR_linear.fig ]
  then
    ${phome}/bin/fig2dev -L pdf ${base}/strict_dir/PFPR_linear.fig ${base}/strict_dir/PFPR_linear.pdf
  fi
