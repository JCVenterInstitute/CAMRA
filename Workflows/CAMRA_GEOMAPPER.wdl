version 1.0

task makeGeomapperFile {
  input {
    File Hamonization
    File Card_data
    String id
    String tax
    String geo
    String date
    
  }
  command <<<
    perl /opt/make_geomapper_file.pl -f ~{Hamonization} -o geomapper.tsv -c ~{Card_data} -s ~{id} -t ~{tax}  -d ~{date} -n ~{geo}
  >>>

  output {
       File geomapper_file = "geomapper.tsv"
  }

  runtime{
    docker: "thclarke/geomapper:latest"
	mem: "20 GiB"
  }
}
workflow make_geomapper_file  {
    meta {
        author: "THomas Clarke"
        email: "tclarke@jcvi.org"
        description: "Make GeoMapper File"
    }

    parameter_meta {
        description: "Analysis of genome, AMR focused."
    }
    input {
            File Hamonization
            File Card_data
            String id
            String tax
            String geo
            String date
    }

    # Task to combine genus and species



    call makeGeomapperFile {
        input:
        	Hamonization = Hamonization,
        	Card_data = Card_data,
        	id = id,
        	tax = tax,
        	geo =geo,
        	date = date
    }
    
    output {
        File mob_contig_out = makeGeomapperFile.geomapper_file

    }
}
