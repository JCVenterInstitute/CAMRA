version 1.0

workflow ncbi_pgap {
    input {
        String sample_id
        File assembly_fasta
        File supplemental_data
        String taxon
        String? ignore_errors
        String? report_usage
        Int? cpu
        Int? memory
        Int? disk_size
        String? last_name
        String? first_name
        String? email
        String? organization
        String? department
        String? phone
        String? fax
        String? street
        String? city
        String? state
        String? postal_code
        String? country
        String? middle_initial
        String? pgap_version
    }

    call pgap_annotate {
        input:
            sample_id = sample_id,
            assembly_fasta = assembly_fasta,
            supplemental_data = supplemental_data,
            taxon=taxon,
            ignore_errors = ignore_errors,
            report_usage = report_usage,
            cpu = cpu,
            memory = memory,
            disk_size = disk_size,
            last_name = last_name,
            first_name = first_name,
            email = email,
            organization = organization,
            department = department,
            phone = phone,
            fax = fax,
            street = street,
            city = city,
            state = state,
            postal_code = postal_code,
            country = country,
            middle_initial = middle_initial,
            pgap_version = pgap_version
    }

    output {
        File? pgap_ani_tax_report = pgap_annotate.pgap_ani_tax_report
        File? pgap_ani_tax_report_xml = pgap_annotate.pgap_ani_tax_report_xml
        File? pgap_annotated_proteins_fasta = pgap_annotate.pgap_annotated_proteins_fasta
        File? pgap_genomic_fasta = pgap_annotate.pgap_genomic_fasta
        File? pgap_gbk_annotation = pgap_annotate.pgap_gbk_annotation
        File? pgap_gff3_annotation = pgap_annotate.pgap_gff3_annotation
        File? pgap_annotated_sequence = pgap_annotate.pgap_annotated_sequence
        File? pgap_annotated_cds = pgap_annotate.pgap_annotated_cds
        File? pgap_translated_cds = pgap_annotate.pgap_translated_cds
        File? pgap_gff_with_sequence = pgap_annotate.pgap_gff_with_sequence
        File? pgap_foreign_sequence = pgap_annotate.pgap_foreign_sequence
        File? pgap_completeness_contamination = pgap_annotate.pgap_completeness_contamination
        File? pgap_log = pgap_annotate.pgap_log
    }
}

task pgap_annotate {
    input {
        String sample_id
        File assembly_fasta
        File supplemental_data
        String taxon
        String ignore_errors = 'true'
        String report_usage = 'false'
        Int cpu = 8
        Int memory = 32
        Int disk_size = 100
        String last_name = 'Doe'
        String first_name = 'Jane'
        String email = 'jane_doe@gmail.com'
        String organization = 'NIH'
        String department = 'NCBI'
        String phone = '301-555-0245'
        String fax = '301-555-1234'
        String street = '9000 Rockville Pike'
        String city = 'Bethesda'
        String state = 'MD'
        String postal_code = '20850'
        String country = 'USA'
        String middle_initial = 'A'
        String pgap_version = '2024-07-18.build7555'
    }

    command <<<

        tar -xzvf ~{supplemental_data}

        if [[ "~{assembly_fasta}" == *.fasta.gz || "~{assembly_fasta}" == *.fa.gz ]]; then 
            gunzip -c ~{assembly_fasta} > assembly.fasta 
        elif [[ "~{assembly_fasta}" == *.fasta || "~{assembly_fasta}" == *.fa ]]; then
            mv ~{assembly_fasta} assembly.fasta
        fi

        python3 <<CODE

        import yaml
        def str_to_bool(s):
            return s.lower() in ["true", "t", "1", "yes", "y"]

        submol = {
            'contact_info': {
                'last_name': '~{last_name}',
                'first_name': '~{first_name}',
                'email': '~{email}',
                'organization': '~{organization}',
                'department': '~{department}',
                'phone': '~{phone}',
                'fax': '~{fax}',
                'street': '~{street}',
                'city': '~{city}',
                'state': '~{state}',
                'postal_code': '~{postal_code}',
                'country': '~{country}'
            },
            'authors': [
                {
                    'author': {
                        'last_name': '~{last_name}',
                        'first_name': '~{first_name}',
                        'middle_initial': '~{middle_initial}'
                    }
                }
            ],
            'organism': {
                'genus_species': '~{taxon}'
            }
        }
        with open('submol.yaml', 'w') as file:
            yaml.dump(submol, file, sort_keys=False, default_flow_style=False)

        

        controller = {
            'fasta': {
                'class': 'File',
                'location': 'assembly.fasta'
            },
            'submol': {
                'class': 'File',
                'location': 'submol.yaml'
            },
            'supplemental_data': {
                'class': 'Directory',
                'location': 'input-~{pgap_version}'
            },
            'report_usage': str_to_bool('~{report_usage}'),
            'ignore_all_errors': str_to_bool('~{ignore_errors}')
        }
        with open('controller.yaml', 'w') as file:
            yaml.dump(controller, file, sort_keys=False, default_flow_style=False)
        
        CODE

        cwltool --timestamps --debug --disable-color --preserve-entire-environment /pgap/pgap-2024-07-18.build7555/pgap.cwl ./controller.yaml
        bash -c 'rename_id="~{sample_id}_"; for file in "$@"; do [ -e "$file" ] && mv -- "$file" "${rename_id}${file}"; done; exit 0' -- ani-tax-report.txt ani-tax-report.xml annot.faa annot.fna annot.gbk annot.gff annot.sqn annot_cds_from_genomic.fna annot_translated_cds.faa annot_with_genomic_fasta.gff calls.tab checkm.txt cwltool.log
        
    >>>

    output {
        File? pgap_ani_tax_report = "~{sample_id}_ani-tax-report.txt"
        File? pgap_ani_tax_report_xml = "~{sample_id}_ani-tax-report.xml"
        File? pgap_annotated_proteins_fasta = "~{sample_id}_annot.faa"
        File? pgap_genomic_fasta = "~{sample_id}_annot.fna"
        File? pgap_gbk_annotation = "~{sample_id}_annot.gbk"
        File? pgap_gff3_annotation = "~{sample_id}_annot.gff"
        File? pgap_annotated_sequence = "~{sample_id}_annot.sqn"
        File? pgap_annotated_cds = "~{sample_id}_annot_cds_from_genomic.fna"
        File? pgap_translated_cds = "~{sample_id}_annot_translated_cds.faa"
        File? pgap_gff_with_sequence = "~{sample_id}_annot_with_genomic_fasta.gff"
        File? pgap_foreign_sequence = "~{sample_id}_calls.tab"
        File? pgap_completeness_contamination = "~{sample_id}_checkm.txt"
        File? pgap_log = "~{sample_id}_cwltool.log"
    } 

    runtime {
        docker: "danylmb/pgap:2024-07-18.build7555"
        memory:"~{memory} GB"
        cpu: "~{cpu}"
        disks: "local-disk ~{disk_size} SSD"
        preemptible:  0
        maxRetries: 0
  }
}
