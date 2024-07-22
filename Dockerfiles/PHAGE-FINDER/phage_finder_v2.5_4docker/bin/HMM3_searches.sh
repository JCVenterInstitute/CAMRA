#!/bin/sh

#Bin="$(dirname $(stat -f "$0"))"
Bin="$(cd -- "$(dirname "$0")" && pwd)"
phome=${Bin%/*}
pepfile=$1
src="/scratch/dfouts/bin"
base=`pwd`

if [ ! -d HMM3_searches_dir ] # check if directory is present
then
    mkdir ${base}/HMM3_searches_dir
else
    rm -rf ${base}/HMM3_searches_dir/
    mkdir ${base}/HMM3_searches_dir
fi
cd ${base}/HMM3_searches_dir/
total=`cat ${phome}/hmm3.lst | wc -l | sed 's/^ *//'`
for i in `cat ${phome}/hmm3.lst`
do
  let "count = count + 1"
  Result=`echo $count $total | awk '{printf( "%3.1f\n", (($1/$2) * 100))}'`
  echo -ne "\r    $Result% complete"
  #echo -ne "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b$Result% complete"
  ${src}/hmmsearch ${phome}/PHAGE_HMM3s_dir/${i}.HMM ${base}/${pepfile} >> $i.out
done
echo
cat *.out >> ../combined.hmm3
cd ..
rm -rf HMM3_searches_dir/
