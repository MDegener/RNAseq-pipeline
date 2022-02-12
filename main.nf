#!/usr/bin/env nextflow

// Set date as identifier
date = new java.util.Date().format("ddMMyyyyhhmmss")

// Show used parameters
log.info "RNAseq - N F  ~  version 0.1"
log.info "====================================="
log.info "inFiles                  : ${params.inFiles}"
log.info "outDir                 : ${params.outDir}"
log.info "species                : ${params.species}"
log.info "flow                   : ${params.flow}"
log.info "====================================="
log.info "\n"

/*
 * General inFiles parameters validation
 */
 
envDir                          = params.envDir
libDir                         	= params.libDir
outDir				            = params.outDir
multiqc_file                    = file(params.multiqc)
flow                            = params.flow

if(params.species == "human"){
    star_index                  = file(libDir+"/star_index_gencode_v26/")
    gtf                         = file(libDir+"/gencode.v26.annotation.gtf")
    miso_annotation             = file(libDir+"/miso_SE_indexed_hg38/")
    miso_settings               = file(baseDir+"/config/miso_settings.txt")
}
    
flow = flow.split(',').collect { it.trim() }

if(!("bam" in flow)){
    Channel
        .fromFilePairs(params.inFiles, size: -1)
        .ifEmpty {error"Cannot find any files matching: ${params.inFiles}"}
        .into { read_files_1; read_files_2 }
} else {
    Channel
        .fromFilePairs(params.inFiles, size: -1)
        .ifEmpty {error"Cannot find any files matching: ${params.inFiles}"}
        .into { bam_miso; bam_rnaseqc }
}

  if("fastqc" in flow){
      process fastqc {
          label "fastqc"
		  conda "$envDir"

          publishDir = [path: "$outDir/fastqc", mode: 'copy']
          tag "reads: $sample_id"
  
          input:
          set val(sample_id), file(reads) from read_files_1
  
          output:
          file("${sample_id}*fastqc.html") into fastqc_html
          file("${sample_id}*fastqc.zip") into fastqc_zip
  
          script:
          def single = reads instanceof Path
          
          if(single)
              """
              fastqc ${reads}
              """
          else
              """
              fastqc ${reads[0]} ${reads[1]}
              """
      }
  } else {
   fastqc_zip = ""
  }
  
  if("star" in flow){
      process star {
          label "mapping"
		  conda "$envDir"
          cpus 8
          publishDir = [path: "$outDir/star_output", mode:'copy']
          tag "reads: $sample_id"
  
          input:
          file(star_index) from star_index
          set val(sample_id), file(reads) from read_files_2
  
          output:
          set val(sample_id), file("${sample_id}.Aligned.out.bam") into bam_files
          file("${sample_id}*Log.final.out") into final_log
  
          script:
          def single = reads instanceof Path
          if(single)
              """
              STAR --runThreadN ${task.cpus} --genomeDir ${star_index} --readFilesIn ${reads} --readFilesCommand gunzip -c --outSAMtype BAM Unsorted --outFileNamePrefix "${sample_id}."
              """
          else
              """
              STAR --runThreadN ${task.cpus} --genomeDir ${star_index} --readFilesIn ${reads[0]} ${reads[1]} --readFilesCommand gunzip -c --outSAMtype BAM Unsorted --outFileNamePrefix "${sample_id}."
              """
      }
  } else {
   final_log = ""
  }
  
  if("sort" in flow){
      process sorting{
          label "counting"
		  conda "$envDir"
          cpus 8
          publishDir = [path: "$outDir/sort_index", mode:'copy']
          tag "bam: $sample_id"
          
          input:
            set val(sample_id), file(bam) from bam_files
           
          output:
            set val(sample_id), file("${sample_id}.Aligned.sorted.out.bam") into bam_sorted
            
          script:
              """
              samtools sort --threads ${task.cpus} ${bam} -o ${sample_id}.Aligned.sorted.out.bam
              """
      }
  }

  if("index" in flow){
      process index{
          label "counting"
		  conda "$envDir"
          cpus 8
          publishDir = [path: "$outDir/sort_index", mode:'copy']
          tag "bam: $sample_id"
          
          input:
            set val(sample_id), file(bam_s) from bam_sorted
            
          output:
            set val(sample_id), file("${sample_id}.Aligned.sorted.out.bam.bai") into bam_indexed
            
          script:
              """
              samtools index -@ ${task.cpus} ${bam_s}
              """
      }
  }

if("miso" in flow){

  process miso{
    label "miso"
	conda "$envDir"
    cpus 2
    publishDir = [path: "$outDir/miso/output", mode:'copy']
    tag "bam: $sample_id"
    
    input:
      set val(sample_id), file(bam) from bam_miso
 
    output:
      path "./${sample_id}" into miso_output
      
    script:
      """
      cp ${miso_settings} ./miso_settings.txt
      sed -i 's/num_processors = [0-9]*/num_processors = ${task.cpus}/g' ./miso_settings.txt
      miso --settings-filename="miso_settings.txt" \
      --run ${miso_annotation} ${sample_id}*.bam \
      --output-dir $sample_id \
      --read-len ${params.readLength} 
      """
  }

  process summarize_miso{
    label "summarize_miso"
	conda "$envDir"
    publishDir = [path: "$outDir/miso", mode:'copy']
    tag "bam: $sample_id"
    
    input:
      path sample_id from miso_output
      
    output:
      file "summary/${sample_id}.miso_summary"
      
    script:
      """
      summarize_miso --summarize-samples $sample_id ./
      """
  }
}  

if("rnaseqc" in flow){
    process rnaseqc{
	    label "rnaseqc"
		conda "$envDir"
        cpus 1
        tag "bam: $sample_id"
        publishDir = [path: "$outDir/rnaseqc_output", mode: 'copy']

        input:
        set val(sample_id), file(bam) from bam_rnaseqc
        file(gtf) from gtf

        output:
         path "./${sample_id}" into (rnaseqc_output)

        script:
          if(params.paired == "no")
            """
            rnaseqc --unpaired ${gtf} ${sample_id}*.bam ${sample_id}
            """
          else
            """
            rnaseqc ${gtf} ${sample_id}*.bam ${sample_id}
            """
    }
}

if("multiqc" in flow){
  process multiqc {
      label "report"       
	  conda "$envDir"
      publishDir = [path: "$outDir/multiqc_report", mode:'copy']
  
      input:
      file('*') from fastqc_zip.collect()
      file('*') from final_log.collect()
      file(config) from multiqc_file
  
      output:
      file('multiqc_report.html')
  
      script:
          """
          cp $config/* .
          multiqc .
          """
  }
}

workflow.onComplete {
    println ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}