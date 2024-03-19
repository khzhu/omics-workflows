version 1.0

## This WDL annotates a cohort VCF jointly called from human DNAs using VEP and
## subsets individual sample from the cohort VCF using GATK4
##
## Requirements/expectations :
## - Sample or cohort VCF produced by GATK4 and its index file
##
## Output :
## - VEP annotated VCF file and its index, containing variants jointly called from GATK4
## - Individual annotated VCF files and their index files
##
## Software version requirements (see recommended dockers in inputs JSON)
## - GATK 4 or later (see gatk docker)
## - Ensembl Variant Effect Predictor (VEP)

# WORKFLOW DEFINITION
workflow VEPAnnotationWF {
	input {
        File input_vcf
        File vcf_index
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File vep_cache
        File vep_tar_ball
        File sample_name_map
        String docker_image = "kzhupm/perl:5.36.1.3"
        String gatk_docker = "kzhupm/gatk:4.4.0.0"
	}
  Array[String] samples = read_lines(sample_name_map)
  call CallVEP {
      input:
        input_vcf = input_vcf,
        vcf_index = vcf_index,
        ref_fasta = ref_fasta,
        ref_fai = ref_fasta_index,
        docker_image = docker_image,
        vep_cache = vep_cache,
        vep_tar_ball = vep_tar_ball
  }

  scatter (sample in samples) {
      call SelectSampleVariants {
        input:
          input_vcf = CallVEP.output_vcf,
          input_vcf_idx = CallVEP.output_tbi,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          sample_name = sample,
          gatk_docker = gatk_docker
      }
  }

  output {
      File vep_vcf = CallVEP.output_vcf
      File vep_tbi = CallVEP.output_tbi
      Array[File] sample_vcf = SelectSampleVariants.output_vcf
      Array[File] sample_vcf_idx = SelectSampleVariants.output_vcf_idx
	}

  meta {
      author: [
          {
              name: "Kelsey Zhu",
          }
      ]
  }
}

# TASK DEFINITION
task CallVEP {
    input {
        File input_vcf
        File vcf_index
        File ref_fasta
        File ref_fai
        File vep_cache
        File vep_tar_ball
        String docker_image
        Int preemptible = 2
    }
    String callset_name = basename(input_vcf, ".vcf.gz")
    String vep_tar = basename(vep_tar_ball)
    String vep_cache_dir = basename(vep_cache)

    command {
        set -euox pipefail
        cp ~{vep_tar_ball} .
        tar -zxf ~{vep_tar}
        cp ~{vep_cache} .
        tar -zxf ~{vep_cache_dir}

        ./ensembl-vep-release-109.3/vep -i ~{input_vcf} --assembly GRCh38 --cache --offline --vcf --sift b --polyphen b \
        --ccds --uniprot --hgvs  --symbol --numbers --domains --regulatory --canonical --protein --biotype \
        --af --max_af --af_gnomadg --af_gnomade --mane \
        --variant_class --gene_phenotype --pubmed --shift_hgvs 0 --allele_number --total_length --no_escape --format vcf \
        --force --buffer_size 100000 --no_stats --dir ./ --fasta ~{ref_fasta} \
        --fork 12 -o ~{callset_name}.vep.vcf \
        --plugin pLI,./ensembl-vep-release-109.3/Plugins/pLI_values.txt \
        --dir_plugins ./ensembl-vep-release-109.3/Plugins/

        ./ensembl-vep-release-109.3/htslib/bgzip -c ~{callset_name}.vep.vcf > ~{callset_name}.vep.vcf.gz
        ./ensembl-vep-release-109.3/htslib/tabix -p vcf ~{callset_name}.vep.vcf.gz
    }

    runtime {
        docker: docker_image
        memory: "212 GiB"
        cpu: 12
        maxRetries: preemptible
    }

    output {
        File output_vcf = "~{callset_name}.vep.vcf.gz"
        File output_tbi = "~{callset_name}.vep.vcf.gz.tbi"
    }
}

task SelectSampleVariants {
  input {
    File input_vcf
    File input_vcf_idx
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    String sample_name
    String gatk_docker
  }

  parameter_meta {
    input_vcf: {
      localization_optional: true
    }
  }

  command <<<
    set -euo pipefail
    gatk --java-options "-Xms12000m -Xmx14000m" \
      SelectVariants \
      -V ~{input_vcf} \
      -R ~{ref_fasta} \
      -sn ~{sample_name} \
      --exclude-non-variants true \
      -O ~{sample_name}.vcf.gz
  >>>

  runtime {
    memory: "15000 MiB"
    cpu: 2
    maxRetries: 1
    docker: gatk_docker
  }

  output {
    File output_vcf = "~{sample_name}.vcf.gz"
    File output_vcf_idx = "~{sample_name}.vcf.gz.tbi"
  }
}
