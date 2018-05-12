
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
    File rawdata_fastqR1
    File rawdata_fastqR2
    String sample_name

    call AdapterRemoval_Trim {
        input:
            AdapterRemoval = AdapterRemoval,
            sample_name = sample_name,
            in_fastqR1 = rawdata_fastqR1,
            in_fastqR2 = rawdata_fastqR2
    }

    call rnaseq_star_map {
        input:
            STAR = STAR,
            STARindexDir = STARindexDir,
            sample_name = sample_name,
            in_fastqR1 = in_fastqR1,
            in_fastqR2 = in_fastqR2
    }

    call rnaseq_featureCounts_quant {
        input:
            featureCounts = featureCounts,
            annotation_file = annotation_file,
            sample_name = sample_name,
            in_bam = in_bam
    }
}

## Tasks
# Trim reads using AdapterRemoval
task AdapterRemoval_Trim {
    File AdapterRemoval
    String sample_name
    File in_fastqR1
    File in_fastqR2
    Int cpu=28

    command {
        ${AdapterRemoval} --file1 ${in_fastqR1}
            --file2 ${in_fastqR2}
            --output1 ${sample_name}_T1.fastq.gz \
            --output2 ${sample_name}_T2.fastq.gz \
            --trimns --trimqualities --gzip
        }

    output {
        File out_R1 = "${sample_name}_T1.fastq.gz"
        File out_R2 = "${sample_name}_T2.fastq.gz"
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
    File in_fastqR1
    File in_fastqR2
    #String suffix="Aligned.sortedByCoord.out"
    #This is STAR's suffix and cannot change:
    String suffix="Aligned.sortedByCoord.out"
    String genome="hg19"
    Int cpu=28

    command {
        ln -s ${STARindexDir} GenomeDir
        # Add alignment directory for outfiles?
        ${STAR} --readFilesIn ${in_fastqR1} ${in_fastqR2} \
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
        File out_bam = "${sample_name}.${genome}.${suffix}.bam"
        String out_sample_name = "${sample_name}.${genome}.${suffix}"
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
