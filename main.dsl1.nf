/**
 * @deprecated Initial version of the workflow. Has been replaced with a workflow
 *             based on the nextflow DSL2. This workflow is still kept here for 
 *             reference and as the analysis is not managed in a VCS.
 *
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

params.contigs = "data/p.cubensis/assembly/Psicub1_1/Psicub1_1_AssemblyScaffolds.fasta.gz"
params.reads = "data/p.cubensis/rnaseq/**/*.fastq.gz"
params.species = "Psilocybe cubensis"
params.cpus = 10
params.contigs_minlen = 1000

Channel.fromPath(params.contigs).set{ch_contigs}
Channel.fromPath(params.reads).set{ch_reads}

process gunzip_contigs { 

  input:
  path contigs from ch_contigs

  output:
  path "${contigs.baseName}" into ch_unzipped_contigs

  script:
  """
  zcat $contigs > ${contigs.baseName}
  """
}

process fix_ena_headers_for_trinity {
  input:
  path reads from ch_reads

  output:
  path "${reads.simpleName}.corrected_headers.fastq.gz" into ch_corrected_reads

  script:
  """
  zcat $reads | sed  -r 's#^(@.+)\\.[0-9]+ ([0-9]+)/([0-9]+)#\\1_\\2/\\3#' | gzip -c > ${reads.simpleName}.corrected_headers.fastq.gz
  """
}

process clean_contigs {

  input:
  path contigs from ch_unzipped_contigs

  output:
  path "${contigs.simpleName}.cleaned.fa" into ch_cleaned_contigs

  script:
  """
  funannotate clean -i $contigs --minlen ${params.contigs_minlen} -o ${contigs.simpleName}.cleaned.fa
  """
}

process sort_contigs {

  input:
  path contigs from ch_cleaned_contigs

  output:
  path "${contigs.simpleName}.sorted.fa" into ch_sorted_contigs

  script:
  """
  funannotate sort -i $contigs -b scaffold -o ${contigs.simpleName}.sorted.fa
  """
}

process mask_contigs {

  input:
  path contigs from ch_sorted_contigs

  output:
  path "${contigs.simpleName}.masked.fa" into ch_masked_contigs_1, ch_masked_contigs_2

  script:
  """
  funannotate mask -i $contigs --cpus ${params.cpus} -o ${contigs.simpleName}.masked.fa
  """
}

process train {

  input:
  path contigs from ch_masked_contigs_1
  path reads from ch_corrected_reads.toList()

  output:
  path "fun/training" into ch_trained_data_1, ch_trained_data_2

  script:
  """
  funannotate train -i $contigs -o fun --single $reads --jaccard_clip --species  "${params.species}"  --cpus ${params.cpus} 
  """
}

process predict {

  input:
  path contigs from ch_masked_contigs_2
  path "fun/training" from ch_trained_data_1

  output:
  path "fun/*" into ch_prediction

  script:
  """
  funannotate predict -i $contigs -o fun --species  "${params.species}"  --cpus ${params.cpus} 
  """
}

process annotate_utrs {

  input:
  path "fun/" from ch_prediction
  path "fun/training" from ch_trained_data_2

  output:
  path "fun/*" into ch_annotated_utrs includeInputs true

  script:
  """
  funannotate update -i fun --cpus ${params.cpus} 
  """

}

process annotate_function {

  input:
  path "fun/" from ch_annotated_utrs

  output:

  script:
  """
  funannotate annotate -i fun --cpus ${params.cpus} 
  """


}
