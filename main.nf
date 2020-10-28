#!/usr/bin/env nextflow

// Set date as identifier
date = new java.util.Date().format("ddMMyyyyhhmmss")

// Show used parameters
log.info "RNAseq - N F  ~  version 0.1"
log.info "====================================="
log.info "inDir                  : ${params.inDir}"
log.info "outDir                 : ${params.outDir}"
log.info "species                : ${params.species}"
log.info "flow                   : ${params.flow}"
log.info "====================================="
log.info "\n"

/*
 * General inDir parameters validation
 */
 
libDir                         	= params.libDir
outDir				            = params.outDir
multiqc_file                    = file(params.multiqc)
flow                            = params.flow

if(params.species == "human"){
    star_index                  = file(libDir+"/star_index_gencode_v26/")
    gtf                         = file(libDir+"/gencode.v26.annotation.gtf")
    miso_annotation             = file(libDir+"/miso_SE_indexed_hg38/")
    miso_settings               = file(baseDir+"/miso_settings.txt")
}
    
flow = flow.split(',').collect { it.trim() }

if(!("bam" in flow)){
    Channel
        .fromFilePairs(params.inDir, size: -1)
        .ifEmpty {error"Cannot find any files matching: ${params.inDir}"}
        .into { read_files_1; read_files_2 }
} else {
    Channel
        .fromFilePairs(params.inDir, size: -1)
        .ifEmpty {error"Cannot find any files matching: ${params.inDir}"}
        .into { bam_miso; bam_rnaseqc }
}

  if("fastqc" in flow){
      process fastqc {
          label "fastqc"

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
              /opt/cmbi/bin/fastqc ${reads}
              """
          else
              """
              /opt/cmbi/bin/fastqc ${reads[0]} ${reads[1]}
              """
      }
  } else {
   fastqc_zip = ""
  }
  
  if("star" in flow){
      process star {
          label "mapping"
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
              /opt/cmbi/bin/STAR --runThreadN ${task.cpus} --genomeDir ${star_index} --readFilesIn ${reads} --readFilesCommand gunzip -c --outSAMtype BAM Unsorted --outFileNamePrefix "${sample_id}."
              """
          else
              """
              /opt/cmbi/bin/STAR --runThreadN ${task.cpus} --genomeDir ${star_index} --readFilesIn ${reads[0]} ${reads[1]} --readFilesCommand gunzip -c --outSAMtype BAM Unsorted --outFileNamePrefix "${sample_id}."
              """
      }
  } else {
   final_log = ""
  }
  
  if("sort" in flow){
      process sorting{
          label "counting"
          cpus 8
          publishDir = [path: "$outDir/sort_index", mode:'copy']
          tag "bam: $sample_id"
          
          input:
            set val(sample_id), file(bam) from bam_files
           
          output:
            set val(sample_id), file("${sample_id}.Aligned.sorted.out.bam") into bam_sorted
            
          script:
              """
              /opt/cmbi/bin/samtools sort --threads ${task.cpus} ${bam} -o ${sample_id}.Aligned.sorted.out.bam
              """
      }
  }

  if("index" in flow){
      process index{
          label "counting"
          cpus 8
          publishDir = [path: "$outDir/sort_index", mode:'copy']
          tag "bam: $sample_id"
          
          input:
            set val(sample_id), file(bam_s) from bam_sorted
            
          output:
            set val(sample_id), file("${sample_id}.Aligned.sorted.out.bam.bai") into bam_indexed
            
          script:
              """
              /opt/cmbi/bin/samtools index -@ ${task.cpus} ${bam_s}
              """
      }
  }

if("miso" in flow){

  process miso{
    label "miso"
    cpus 2
    publishDir = [path: "$outDir/miso/output", mode:'copy']
    tag "bam: $sample_id"
    
    input:
      set val(sample_id), file(bam) from bam_miso
 
    output:
      path "./${sample_id}" into miso_output
      
    script:
      """
      source /home/maxd/miniconda3/bin/activate /home/maxd/miso-env
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
    publishDir = [path: "$outDir/miso", mode:'copy']
    tag "bam: $sample_id"
    
    input:
      path sample_id from miso_output
      
    output:
      file "summary/${sample_id}.miso_summary"
      
    script:
      """
      source /home/maxd/miniconda3/bin/activate /home/maxd/miso-env
      summarize_miso --summarize-samples $sample_id ./
      """
  }
}  

if("rnaseqc" in flow){
    process rnaseqc{
        cpus 1
        label "rnaseqc"
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
            source /home/maxd/miniconda3/bin/activate /home/maxd/rnaseqc-env
            rnaseqc --unpaired ${gtf} ${sample_id}*.bam ${sample_id}
            """
          else
            """
            source /home/maxd/miniconda3/bin/activate /home/maxd/rnaseqc-env
            rnaseqc ${gtf} ${sample_id}*.bam ${sample_id}
            """
    }
}

if("multiqc" in flow){
  process multiqc {
      label "report"
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