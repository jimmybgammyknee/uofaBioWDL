
## Workflow

workflow standard_rnaseq_quant {

    # Tools
    File AdapterRemoval
    File STAR
    File featureCounts

    # Annotation and indexes
    File annotation_file
    String STARindexDir

    # Data
    # fastqs_R1: fastq.gz files for Read1 (only these if single-ended)
    Array[File] fastqs_R1
    # fastqs_R2: fastq.gz files for Read2 (omit if single-ended) in order
    # corresponding to fastqs_R1
    Array[File] fastqs_R2 = []

    String sample_name
    Array[Array[File]] fastqs_ = if length(fastqs_R2)>0 then transpose([fastqs_R1, fastqs_R2]) else transpose([fastqs_R1])

    scatter (i in range(length(fastqs_))) {
        call AdapterRemoval_Trim {
            input:
                AdapterRemoval = AdapterRemoval,
                sample_name = sample_name,
                fastqs = fastqs_[i]
            }

        call rnaseq_star_map {
            input:
                STAR = STAR,
                STARindexDir = STARindexDir,
                sample_name = sample_name,
                fastqs = AdapterRemoval_Trim.trimmed_merged_fastqs
            }

        call rnaseq_featureCounts_quant {
            input:
                featureCounts = featureCounts,
                annotation_file = annotation_file,
                sample_name = sample_name,
                in_bam = rnaseq_star_map.out_bam
            }
    }
}

## Tasks
# Trim reads using AdapterRemoval
task AdapterRemoval_Trim {
    File AdapterRemoval
    String sample_name
    Array[File] fastqs
    Int cpu=28

    command {
        ${AdapterRemoval} --file1 fastqs[0]
            --file2 fastqs[1]
            --output1 ${sample_name}_T1.fastq.gz \
            --output2 ${sample_name}_T2.fastq.gz \
            --trimns --trimqualities --gzip
        }

    output {
        # nice little hack from the ENCODE guys
        Array[File] trimmed_merged_fastqs = glob("merge_fastqs_R?_*.fastq.gz")
        }

    runtime {
                cpu: cpu
        }
}

# Map RNAseq reads:
task rnaseq_star_map {
    File STAR
    String STARindexDir
    String sample_name
    Array[File] fastqs
    #String suffix="Aligned.sortedByCoord.out"
    #This is STAR's suffix and cannot change:
    String suffix="Aligned.sortedByCoord.out"
    String genome="hg19"
    Int cpu=28

    command {
        ln -s ${STARindexDir} GenomeDir
        # Add alignment directory for outfiles?
        ${STAR} --readFilesIn fastqs[0] fastqs[1] \
                --readFilesCommand zcat \
                --outFilterType BySJout \
                --alignSJoverhangMin 8 \
                --alignSJDBoverhangMin 1 \
                --outFilterMismatchNmax 999 \
                --outSAMtype BAM SortedByCoordinate \
                --outFileNamePrefix ${sample_name} \
                --outSAMattrRGline ID:${sample_name} LB:library PL:illumina PU:machine SM:hg19 \
                --twopassMode Basic \
                --outSAMmapqUnique 60 \
                --runThreadN ${cpu}
        }

    output {
        File out_bam = glob("*.bam")[0]
    }
    runtime {
        cpu: cpu
  }
}

## featureCounts
task rnaseq_featureCounts_quant {
    File featureCounts
    File annotation_file
    String sample_name
    File in_bam
    Int cpu=28

    command {
        ${featureCounts} -a ${annotation_file} \
            -T ${cpu} \
            -s 1 \
            -o ${sample_name}.geneCounts.tsv \
            ${in_bam} > ${sample_name}.featureCounts.log
        }

    output {
        File out_counts = "${sample_name}.geneCounts.tsv"
        String out_log = "${sample_name}.featureCounts.log"
    }

    runtime {
        cpu: cpu
    }
}
