#!/bin/sh

# Bash script to run the Phage_Finder pipeline with HMMer3 models

# Usage: phage_finder.sh <prefix of .pep/.ffa, .ptt and .con/.fna file>

# NOTE: a phage_finder_info.txt file will be searched before a .ptt file
# .pep is the multifasta protein sequence file
# .ptt is a GenBank .ptt file that has the coordinates and ORF names with annotation
# .con is a file that contains the complete nucleotide sequence of the genome being searched

#very cool Perl FindBin functional equivalent for BASH by Nicholas Dronen of crashingdaily.wordpress.com
#for Linux
#$Bin="$(dirname $(readlink -f -- "$0"))"
#for OSX (readlink does not work for OSX)
Bin="$(cd -- "$(dirname "$0")" && pwd)"
phome=${Bin%/*}
base=`pwd`
# set the path to your BLAST+ package
src="/opt"
blastbin="${src}/ncbi-blast-2.15.0+/bin"

prefix=$1

  if [ -s $base/$prefix.pep ] # check if .pep file is present
  then
      pepfile="$prefix.pep"
  elif [ -s $base/$prefix.faa ]
  then
      pepfile="$prefix.faa"
  else
     echo "Could not find a $prefix.pep or $prefix.faa file.  Please check to make sure the file is present and contains data"
     exit 1
  fi 
    if [ -s $base/phage_finder_info.txt ] # check for phage_finder info file and if it has contents
    then
        infofile="phage_finder_info.txt"
    elif [ -s $base/$prefix.ptt ]
    then
          infofile="$prefix.ptt"
    else
      echo "Could not find a phage_finder_info.txt file or $prefix.ptt file.  Please make sure one of these files is present and contains data."
      exit 1
    fi
    if [ ! -s $base/combined.hmm3 ] # if GLOCAL HMM results not present, search
    then
        ## conduct GLOCAL HMMER3 searches
        echo "  HMMER3 searches ..."
        $phome/bin/HMM3_searches.sh $pepfile
    fi
    if [ ! -s $base/ncbi.out ] # if BLAST results not present, search
    then
        ## do NCBI BLASTP searches
        echo "  BLASTing $pepfile against the Phage DB ..."
	# NCBI BLAST+ code provided by Dr. Camille Hankel at Dartmouth #
	$blastbin/blastp -db $phome/DB/phage_03_25_19.db -outfmt 6 -evalue 0.001 -query $pepfile -out ncbi.out -max_target_seqs 5 -num_threads 2
        #old#blastall -p blastp -d $phome/DB/phage_01_03_14_internal.db -m 8 -e 0.001 -i $pepfile -o ncbi.out -v 4 -b 4 -a 2 -F F
    fi
    if [ -s $base/$prefix.con ]
    then
        contigfile="$prefix.con"
    elif [ -s $base/$prefix.fna ]
    then
        contigfile="$prefix.fna"
    else
        echo "Could not find a phage_finder_info.txt file or $prefix.ptt file.  Please make sure one of these files is present and contains data.  In the meantime, I will go ahead and run phage_finder.pl without this information, but beware... NO att sites will be found!"
        contigfile=""
    fi
    if [ ! -e $base/tRNAscan.out ] && [ $base/$contigfile ] # if tRNAscan.out file not present, and contig file present, then search
    then
        ## find tRNAs
        echo "  find tRNA sequences ..."
        $src/tRNAscan-SE -B -o tRNAscan.out $base/$contigfile > /dev/null 2>&1
    fi

    if [ ! -e $base/tmRNA_aragorn.out ] && [ -e $base/$contigfile ] # if tRNAscan.out file not present, and contig file present, then search
    then
        ## find tmRNAs
        echo "  find tmRNA sequences ..."
        $src/aragorn -m -o tmRNA_aragorn.out $base/$contigfile
    fi

    if [ ! -e $base/ncRNA_cmscan.out ] && [ -e $base/$contigfile ] # if ncRNA_cmscan.out file not present, and contig file present, then search
    then
        ## find select ncRNAs
        echo "  find ncRNA sequences ..."
        Z=`$src/esl-seqstat $prefix.con | perl -ne 'chomp; if (/^Total # residues:\s+(\d+)/) {$n = $1; $Z=($n*2)/1000000; print "$Z\n";}'`
        #echo "Z = $Z"        
        $src/cmscan -Z $Z --cut_ga --rfam --nohmmonly --tblout $base/ncRNA_cmscan.out --fmt 2 $phome/RfamDB/Rfam_PhageFinder.cm $prefix.con > $prefix.cmscan
    fi

    ## find the phage
    echo "  searching for Prophage regions ..."
    $phome/bin/Phage_Finder_v2.5.pl -t ncbi.out -i $infofile -r tRNAscan.out -n tmRNA_aragorn.out -c ncRNA_cmscan.out -A $contigfile -S

    ## generate linear figures of each prophage region

    echo "  generating linear illustrations..."
    $phome/bin/LinearDisplay_detailed.pl -A strict_dir/PFPR.frag -F strict_dir/PFPR.att -L -nt > strict_dir/PFPR_linear.fig
    $phome/bin/fig2dev -L pdf strict_dir/PFPR_linear.fig strict_dir/PFPR_linear.pdf
