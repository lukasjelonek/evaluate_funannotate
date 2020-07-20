/**
 * This workflow is based on the tutorial of funannotate:
 * https://funannotate.readthedocs.io/en/latest/tutorials.html
 *
 * In its default configuration it processes the data of Psylocybe cubensis
 * 
 * This workflow was created as an attempt to make this analysis reproducible. 
 * The plan is to extend it as more genomes are processes or other datatypes 
 * are incorporated.
 * 
 * AUTHOR: Lukas Jelonek (Lukas.Jelonek@computational.bio.uni-giessen.de)
 */

nextflow.preview.dsl = 2

params.contigs = "data/p.cubensis/assembly/Psicub1_1/Psicub1_1_AssemblyScaffolds.fasta.gz"
params.reads = "data/p.cubensis/rnaseq/**/*.fastq.gz"
params.species = "Psilocybe cubensis"
params.cpus = 10
params.contigs_minlen = 1000
params.resultdir = "results"

Channel.fromPath(params.contigs).set{ch_contigs}
Channel.fromPath(params.reads).set{ch_reads}

process gunzip_contigs { 

  input:
  path contigs

  output:
  path "${contigs.baseName}" 

  script:
  """
  zcat $contigs > ${contigs.baseName}
  """
}

process fix_ena_headers_for_trinity {
  input:
  path reads 

  output:
  path "${reads.simpleName}.corrected_headers.fastq.gz" 

  script:
  """
  zcat $reads | sed  -r 's#^(@.+)\\.[0-9]+ ([0-9]+)/([0-9]+)#\\1_\\2/\\3#' | pigz -c > ${reads.simpleName}.corrected_headers.fastq.gz
  """
}

process clean_contigs {

  input:
  path contigs 

  output:
  path "${contigs.simpleName}.cleaned.fa" 

  script:
  """
  funannotate clean -i $contigs --minlen ${params.contigs_minlen} -o ${contigs.simpleName}.cleaned.fa
  """
}

process sort_contigs {

  input:
  path contigs 

  output:
  path "${contigs.simpleName}.sorted.fa" 

  script:
  """
  funannotate sort -i $contigs -b scaffold -o ${contigs.simpleName}.sorted.fa
  """
}

process mask_contigs {

  input:
  path contigs 

  output:
  path "${contigs.simpleName}.masked.fa" 

  script:
  """
  funannotate mask -i $contigs --cpus ${params.cpus} -o ${contigs.simpleName}.masked.fa
  """
}

process train {

  input:
  path contigs 
  path reads 

  output:
  path "fun/training" 

  script:
  """
  funannotate train -i $contigs -o fun --single $reads --jaccard_clip --species  "${params.species}"  --cpus ${params.cpus} 
  """
}

process predict {

  input:
  path contigs 
  path "fun/training" 

  output:
  path "fun/{training,predict_misc,predict_results}" includeInputs true

  script:
  """
  funannotate predict -i $contigs -o fun --species  "${params.species}"  --cpus ${params.cpus} 
  """
}

process annotate_utrs {

  input:
  path "fun/" 

  output:
  path "fun/{training,predict_misc,predict_results,update_misc,update_results}" includeInputs true

  script:
  """
  funannotate update -i fun --cpus ${params.cpus} 
  """

}

process annotate_function {

  publishDir params.resultdir, mode: 'copy', saveAs: {x-> x.replaceAll("fun/annotate_results/","")}
  input:
  path "fun/" 

  output:
  path "fun/annotate_results/*"

  script:
  """
  funannotate annotate -i fun --cpus ${params.cpus} 
  """

}

workflow {

  main:
  gunzip_contigs(ch_contigs) | clean_contigs | sort_contigs | mask_contigs
  fix_ena_headers_for_trinity(ch_reads)
  train(mask_contigs.out, fix_ena_headers_for_trinity.out.toList())
  predict(mask_contigs.out, train.out) | annotate_utrs | annotate_function

}
