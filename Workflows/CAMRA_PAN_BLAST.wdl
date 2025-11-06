version 1.0

task makeFastaFiles {
  input {
    Array[File] gb_files
    String db_type
  }
  command <<<
     mkdir ./gb_dir/
     mkdir ./fasta_dir/
     mkdir ./pep/
     mkdir ./nucl/
     echo ~{sep = ' ' gb_files}
     #for amr_file in ${file_paths[@]}; do
     for fl in ~{sep = ' ' gb_files}; do
        echo $fl
        cp $fl ./
        #echo "./"$(basename $fl) >> genomes.list
        echo $(pwd)
        echo $(basename $fl) | awk -F'\.' -v dir="$(pwd)" '{ print $1"\t"dir"/"$0; }' >> genomes.list
        echo $(basename $fl) | awk -F'\.' -v dir="$(pwd)" '{ echo $1"\t"dir"/"$0; }'
     done
     echo "Start\n"
     /pangenome/bin/core_hmm_checker_prep.pl -g genomes.list  -o ./
     touch all_sequences.fasta
     cd pep
     for f in *pep; do
        echo $f
        chmod +rw $f
        if [ $f != "genomes.pep" ]; then
            cat "$f" >> ../all_sequences.fasta
            mv $f ../
        fi
     done
     echo "END\n"
     cd ..
     makeblastdb -in all_sequences.fasta -dbtype ~{db_type} -out blast_db
     find "$(pwd -P)" -type f -name "blast_db*" > blast_list.txt
     find "$(pwd -P)" -type f -name "*pep" > fasta_list.txt
  >>>

  output {
        Array[File] blast_files = glob("blast_db*")
        String blast_db_prefix = "blast_db"
		File genomes_list = "genomes.list"
        Array[File] input_fastas = glob("*pep")
  }

  runtime{
        docker: "thclarke/pangenomepipeline:latest"
	mem: "20 GiB"
  }
}



task RunBlast {
  input {
    File query_fasta
    String blast_db_prefix
    String db_type
    Array[File] blast_files 
  }

  command <<<
    for fl in ~{sep = " " blast_files}; do
        cp $fl ./
        echo $fl
    done 
    
     blastp -query ~{query_fasta} -db ~{blast_db_prefix} -out combined.blast -outfmt 6
  >>>

  output {
    File blast_output = "combined.blast"
  }

  runtime {
    docker: "thclarke/pangenomepipeline:latest"
  }
}

task makeCombinedBlast {
  input {
    Array[File]  blast_output
    Int hdd_size = 10
  }

  command <<<
     cat ~{ sep = ' '  blast_output} > combined.blast
  >>>

  output {
    File blast_output = "combined.blast"
  }

  runtime {
    docker: "thclarke/pangenomepipeline:latest"
    disks: "local-disk ~{ hdd_size } HDD"  # Request 4 GB of memory
  }
}

workflow make_blast   {
    meta {
        author: "Thomas Clarke"
        email: "tclarke@jcvi.org"
        description: "Run Pangenome."
    }

    input {
        Array[File] gb_files
        String db_type
	Int gb_req
    }


    call makeFastaFiles{
     input:
        gb_files = gb_files,
        db_type = db_type
  }


  scatter(fasta_file in makeFastaFiles.input_fastas) {
    call RunBlast {
      input:
        query_fasta = fasta_file,
        blast_db_prefix = makeFastaFiles.blast_db_prefix,
        blast_files = makeFastaFiles.blast_files,
	db_type = db_type
      }
    }

    call makeCombinedBlast {
        input:
          blast_output = RunBlast.blast_output

    }

    output {
     File combined_blast_output = makeCombinedBlast.blast_output
   }

}
