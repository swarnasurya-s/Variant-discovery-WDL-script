task adapter_trimming{


File R1_input
File R2_input

File adapter_fasta

command {
java -jar /home/smc/Desktop/Swarna/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 \
${R1_input} ${R1_input} \
AT.R1.fastq /dev/null AT.R2.fastq /dev/null \
ILLUMINACLIP:${adapter_fasta}:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
}

output {
File R1_output = "AT.R1.fastq"
File R2_output = "AT.R2.fastq"
}
}

task count_lines{

File trimmed_R1
File trimmed_R2

command {
wc -l ${trimmed_R1} > R1_num_lines.txt ; wc -l ${trimmed_R1} > R2_num_lines.txt
}

}


task BWATest {

File input_AT_R1_fastq
File input_AT_R2_fastq
File reference_fasta
  File ref_index
  File ref_dict
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa

command {
/home/smc/Desktop/Swarna/bwa mem -aM -R "@RG\tID:Humanp\tSM:Human\tPL:ILLUMINA\tPI:330" ${reference_fasta} ${input_AT_R1_fastq} ${input_AT_R2_fastq} > output.sam
}              

output {
        File test_output = "output.sam"
}
}


task sambam {
File input_sam_file_name


command {
/home/smc/Desktop/Swarna/sambamba-0.7.1-linux-static view --format=bam -h -t 16 -S ${input_sam_file_name} -o output.bam
}

output {
File samtobam_output= "output.bam"
}
}

task sort_bam {
File bam_out


command {
/home/smc/Desktop/Swarna/sambamba-0.7.1-linux-static sort -t 16 -o sorted_bam.bam ${bam_out}
}
 
output {
File sorted_bam= "sorted_bam.bam"
File sorted_bam_index= "sorted_bam.bam.bai"
}
}



task haplotypeCaller {

  
  File RefFasta
  File RefIndex
  File RefDict
  File inputBAM
  File bamIndex
  command {
     java -jar /home/smc/Desktop/Swarna/gatk-4.1.6.0/gatk-package-4.1.6.0-local.jar  \
        HaplotypeCaller \
        -R ${RefFasta} \
        -I ${inputBAM} \
        -O human.raw.indels.snps.vcf 
  }
runtime {
  continueOnReturnCode: 2
}
  output {
    File rawVCF = "human.raw.indels.snps.vcf"
  }
}


task selectSNPs {

  File RefFasta
  File RefIndex
  File RefDict
  String type
  File rawVCF

  command {
    java -jar /home/smc/Desktop/Swarna/gatk-4.1.6.0/gatk-package-4.1.6.0-local.jar \
      SelectVariants \
      -R ${RefFasta} \
      -V ${rawVCF} \
      -select-type ${type} \
      -O human_raw.snps.vcf
  }
runtime {
  continueOnReturnCode: 0
}
  output {
    File rawSubset = "human_raw.snps.vcf"
  }
}



task selectIndels {


  File RefFasta
  File RefIndex
  File RefDict
  String type
  File rawVCF

  command {
    java -jar /home/smc/Desktop/Swarna/gatk-4.1.6.0/gatk-package-4.1.6.0-local.jar \
      SelectVariants \
      -R ${RefFasta} \
      -V ${rawVCF} \
      -select-type ${type} \
      -O human_raw.indels.vcf
  }
runtime {
  continueOnReturnCode: 0
}
  output {
    File rawSubset = "human_raw.indels.vcf"
  }
}



task hardFilterSNP {

 
  File RefFasta
  File RefIndex
  File RefDict
  File rawSNPs
  
  command {
    java -jar /home/smc/Desktop/Swarna/gatk-4.1.6.0/gatk-package-4.1.6.0-local.jar \
      VariantFiltration \
      -R ${RefFasta} \
      -V ${rawSNPs} \
      --filter-expression "FS > 60.0" \
      --filter-name "snp_filter" \
      -O human.filtered.snps.vcf
  }
runtime {
  continueOnReturnCode: 0
}
  output {
    File filteredSNPs = "human.filtered.snps.vcf"
  }
}

task hardFilterIndel {

 
  File RefFasta
  File RefIndex
  File RefDict
  File rawIndels

  command {
    java -jar /home/smc/Desktop/Swarna/gatk-4.1.6.0/gatk-package-4.1.6.0-local.jar \
       VariantFiltration \
      -R ${RefFasta} \
      -V ${rawIndels} \
      --filter-expression "FS > 200.0" \
      --filter-name "indel_filter" \
      -O human.filtered.indels.vcf
  }
runtime {
  continueOnReturnCode: 0
}
  output {
    File filteredIndels = "human.filtered.indels.vcf"
  }
}

task combine {

  File RefFasta
  File RefIndex
  File RefDict
  File filteredSNPs
  File filteredIndels

  command {
    /home/smc/Desktop/Swarna/jdk1.8.0_241/bin/java -jar /home/smc/Desktop/Swarna/GATK-3.5.0/GenomeAnalysisTK.jar \
      -T CombineVariants \
      -V ${filteredSNPs} \
      -R ${RefFasta} \
      -V ${filteredIndels} \
      --genotypemergeoption UNSORTED \
      -o human.filtered.snps.indels.vcf
  } 

  output {
    File filteredVCF = "human.filtered.snps.indels.vcf"
  }
}

task table{
File sorted_bam
File bamIndex
File RefFasta
File fastaindex
File refdict
File chr20_vcf
File human_vcf
command {
      /home/smc/Desktop/Swarna/jdk1.8.0_241/bin/java -jar /home/smc/Desktop/Swarna/GATK-3.5.0/GenomeAnalysisTK.jar  \
     -T BaseRecalibrator \
      -R ${RefFasta} \
      -I ${sorted_bam} \
      --knownSites ${chr20_vcf} \
      --knownSites ${human_vcf}  \
      -o recal.grp
}
runtime {
  continueOnReturnCode: 0
}
output {
 File table_output="recal.grp"
}
}



task bqsr {
File sorted_bam
File bamIndex
File RefFasta
File fastaindex
File refdict
File chr20_vcf
File human_vcf
File recal_data

command {
      /home/smc/Desktop/Swarna/jdk1.8.0_241/bin/java -jar /home/smc/Desktop/Swarna/GATK-3.5.0/GenomeAnalysisTK.jar  \
      -T PrintReads \
      -R ${RefFasta} \
      -I ${sorted_bam} \
      -filterNoBases \
     -BQSR ${recal_data} \
      -o output_bqsr.bam
}
runtime {
  continueOnReturnCode: 0
}
output {
File bqsr_output= "output_bqsr.bam"
}
}

workflow bqsr_recalibration{

  call adapter_trimming
  call count_lines {input: trimmed_R1=adapter_trimming.R1_output, trimmed_R2=adapter_trimming.R2_output }
  call BWATest {input: input_AT_R1_fastq=adapter_trimming.R1_output, input_AT_R2_fastq=adapter_trimming.R2_output } 
  call sambam {input: input_sam_file_name= BWATest.test_output }
  call sort_bam {input: bam_out= sambam.samtobam_output }

  call haplotypeCaller {input: bamIndex=sort_bam.sorted_bam, inputBAM=sambam.samtobam_output }
  call selectSNPs { input:  rawVCF=haplotypeCaller.rawVCF }
  call selectIndels { input:  rawVCF=haplotypeCaller.rawVCF }
  call hardFilterSNP { input:  rawSNPs=selectSNPs.rawSubset }
  call hardFilterIndel { input:  rawIndels=selectIndels.rawSubset }
  call combine { input: filteredSNPs=hardFilterSNP.filteredSNPs, filteredIndels=hardFilterIndel.filteredIndels }
  
  call table { input: sorted_bam=sort_bam.sorted_bam, bamIndex=sort_bam.sorted_bam_index, human_vcf=combine.filteredVCF }
  call bqsr { input: recal_data=table.table_output, sorted_bam=sort_bam.sorted_bam, bamIndex=sort_bam.sorted_bam_index, human_vcf=combine.filteredVCF }
}