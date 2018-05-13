

# WDL Workflow for standard UofA RNAseq pipeline

# Jimmy Breen (jimmymbreen@gmail.com)

## Runs:
##  - Trimming using AdapterRemoval
##  - Mapping using STAR
##  - Pseudomapping using kallisto
##  - Quant using featureCounts

## Validated the workflow using wdltool
## - And made acyclic graph using:
# dot <(wdltool graph standard_rnaseq.wdl) -Tpdf > standard_rnaseq_wdl_validate.pdf

workflow standard_rnaseq_quant {

    # Tools
    File AdapterRemoval
    File STAR
    File kallisto
    File featureCounts

    # Annotation and indexes
    File annotation_file
    String STARindexDir
    File index

    # Data
    File sampleinfo
    String datadir
    Array[String] sample_name = read_lines(sampleinfo)
    File rawfastq_read1
    File rawfastq_read2
    File trim_fastqR1
    File trim_fastqR2

    scatter(idx in range(length(sample_name))) {
        call AdapterRemoval_Trim {
            input:
                AdapterRemoval = AdapterRemoval,
                sample_name = sample_name[idx],
                rawfastq_read1 = datadir+'/'+sample_name[idx]+"_1.fastq.gz",
                rawfastq_read2 = datadir+'/'+sample_name[idx]+"_2.fastq.gz",
            }
        call rnaseq_kallisto {
            input:
                kallisto = kallisto,
                index = index,
                sample_name = sample_name[idx],
                trim_fastqR1 = AdapterRemoval_Trim.out_trim_R1,
                trim_fastqR2 = AdapterRemoval_Trim.out_trim_R2
        }

        call rnaseq_star_map {
            input:
                STAR = STAR,
                STARindexDir = STARindexDir,
                sample_name = sample_name[idx],
                trim_fastqR1 = AdapterRemoval_Trim.out_trim_R1,
                trim_fastqR2 = AdapterRemoval_Trim.out_trim_R2
            }

        call rnaseq_featureCounts_quant {
            input:
                featureCounts = featureCounts,
                annotation_file = annotation_file,
                sample_name = sample_name[idx],
                in_bam = rnaseq_star_map.out_bam
            }
    }
}

## Tasks
# Trim reads using AdapterRemoval
task AdapterRemoval_Trim {
    File AdapterRemoval
    String sample_name
    File rawfastq_read1
    File rawfastq_read2
    Int cpu=1

    command {
        module load AdapterRemoval/2.2.0-foss-2016a
        ${AdapterRemoval} --file1 ${rawfastq_read1} \
            --file2 ${rawfastq_read2} \
            --output1 ${sample_name}_T1.fastq.gz \
            --output2 ${sample_name}_T2.fastq.gz \
            --trimns --trimqualities --gzip
        }

    output {
        # nice little hack from the ENCODE guys
        File out_trim_R1 = "${sample_name}_T1.fastq.gz"
        File out_trim_R2 = "${sample_name}_T2.fastq.gz"
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
    File trim_fastqR1
    File trim_fastqR1
    String suffix="Aligned.sortedByCoord.out"
    String genome="hg19"
    Int cpu=1

    command {
        ln -s ${STARindexDir} GenomeDir
        # Add alignment directory for outfiles?
        ${STAR} --readFilesIn ${trim_fastqR1} ${trim_fastqR2} \
                --readFilesCommand zcat \
                --outFilterType BySJout \
                --alignSJoverhangMin 8 \
                --alignSJDBoverhangMin 1 \
                --outFilterMismatchNmax 999 \
                --outSAMtype BAM SortedByCoordinate \
                --outFileNamePrefix ${sample_name} \
                --outSAMattrRGline ID:${sample_name} LB:library PL:illumina PU:machine SM:${genome} \
                --twopassMode Basic \
                --outSAMmapqUnique 60 \
                --runThreadN ${cpu}
        }

    output {
        File out_bam = "${sample_name}.${genome}.${suffix}.bam"
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
    Int cpu=1

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

task rnaseq_kallisto {
    File kallisto
    File trim_fastqR1
    File trim_fastqR2
    File index
    String sample_name
    Int cpu=4

    command {
    ${kallisto} quant \
      --index "${index}" \
      --output-dir . \
      --bootstrap-samples 100 \
      --threads ${cpu} \
      ${trim_fastqR1} ${trim_fastqR2}
    mv abundance.h5 "${sample_name}.abundance.h5"
    mv abundance.tsv "${sample_name}.abundance.tsv"
    mv run_info.json "${sample_name}.run_info.json"
    }

    runtime {
      cpu: cpu
    }

    output {
        File abundance_h5 = "${sample_name}.abundance.h5"
        File abundance_tsv = "${sample_name}.abundance.tsv"
        File log = "${sample_name}.run_info.json"
    }
}
