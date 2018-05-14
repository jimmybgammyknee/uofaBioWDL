workflow jointCallingGenotypes {
  File inputSamplesFile
  Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
  File refFasta
  File refIndex
  File refDict
  scatter (sample in inputSamples) {
    call HaplotypeCallerERC {
      input: bamFile=sample[1],
        bamIndex=sample[2],
        sampleName=sample[0],
        RefFasta=refFasta,
        RefIndex=refIndex,
        RefDict=refDict
    }
  }
  call GenotypeGVCFs {
    input: GVCFs=HaplotypeCallerERC.GVCF,
        sampleName="CEUtrio",
        RefFasta=refFasta,
        RefIndex=refIndex,
        RefDict=refDict
  }
}
task HaplotypeCallerERC {
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  File bamFile
  File bamIndex
  command {
    java -jar /usr/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -ERC GVCF \
        -R ${RefFasta} \
        -I ${bamFile} \
        -o ${sampleName}_rawLikelihoods.g.vcf
  }
  output {
    File GVCF = "${sampleName}_rawLikelihoods.g.vcf"
  }
  runtime {
    docker: "broadinstitute/gatk3:3.5-0"
    docker_name: "broadinstitute/gatk3:3.5-0"
  }
}
task GenotypeGVCFs {
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  Array[File] GVCFs
  command {
    java -jar /usr/GenomeAnalysisTK.jar \
        -T GenotypeGVCFs \
        -R ${RefFasta} \
        -V ${sep=" -V " GVCFs} \
        -o ${sampleName}_rawVariants.vcf
  }
  output {
    File rawVCF = "${sampleName}_rawVariants.vcf"
  }
  runtime {
    docker: "broadinstitute/gatk3:3.5-0"
    docker_name: "broadinstitute/gatk3:3.5-0"
  }
}
