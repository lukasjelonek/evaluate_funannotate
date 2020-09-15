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
params.proteins = "data/p.cubensis/p.cubensis_extrinsic/ref/*.faa"
params.species = "Psilocybe cubensis"
params.cpus = 10
params.contigs_minlen = 1000
params.resultdir = "results"
params.no_clean = false
params.no_repeat_masking = false
params.split_single_fastq = false

// remove stop codos from the fasta sequences
process cleanup_proteins {
  input:
  path seq

  output:
  path "${seq.simpleName}.c.faa"

  script:
  """
  sed 's/[.*]\$//' $seq > ${seq.simpleName}.c.faa
  """

}

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

process split_single_fastq {

  input:
  path reads

  output:
  path "${reads.simpleName}_*.fastq.gz"

  script:
  """
  split_mixed_fastq.groovy -f $reads -p ${reads.simpleName} -r
  """
}

process clean_contigs {
  
  //conda "${workflow.projectDir}/funannotate-env.yaml"

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
  //conda "${workflow.projectDir}/funannotate-env.yaml"

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
  //conda "${workflow.projectDir}/funannotate-env.yaml"

  input:
  path contigs 

  output:
  path "${contigs.simpleName}.masked.fa" 

  script:
  """
  funannotate mask -i $contigs --cpus ${params.cpus} -o ${contigs.simpleName}.masked.fa
  """
}

String getFilePrefix(def s) {
  return s.toString().replaceAll(/(_1|_2)?\.fastq(\.gz)?/, "")
}

// tests for getFilePrefix
assert "SRR7028479" == getFilePrefix("SRR7028479_1.fastq.gz")
assert "SRR7028479" == getFilePrefix("SRR7028479_2.fastq.gz")
assert "SRR7028479" == getFilePrefix("SRR7028479.fastq.gz")
assert "SRR7028479" == getFilePrefix("SRR7028479_1.fastq")
assert "SRR7028479" == getFilePrefix("SRR7028479_2.fastq")
assert "SRR7028479" == getFilePrefix("SRR7028479.fastq")

process train {
  //conda "${workflow.projectDir}/funannotate-env.yaml"

  input:
  // a fasta file with contigs/scaffolds/other nucleotide sequences
  path contigs 
  // A list of reads. They will be separated into single end files 
  // and paired-end file-pairs according to their prefix and the presence
  // of _1|_2 in their filename.
  path reads 

  output:
  path "fun/training" 

  when:
  reads.size() > 0

  script:
  // group reads on prefix
  all_reads = reads.groupBy{x -> getFilePrefix(x)}
  se_reads = all_reads.findAll{k,v -> v.size() == 1}.collect{k,v -> v}
  pe_reads = all_reads.findAll{k,v -> v.size() == 2}.collect{k,v -> v}

  pe1_reads = pe_reads.collect{x -> x[0]}
  pe2_reads = pe_reads.collect{x -> x[1]}

  read_params = '' 
  if (se_reads) read_params += "--single ${se_reads.join(' ')}"
  if (pe_reads) read_params += "--left ${pe1_reads.join(' ')} --right ${pe2_reads.join(' ')}"
  """
  funannotate train -i $contigs -o fun $read_params --jaccard_clip --species  '${params.species}'  --cpus ${params.cpus} 
  """
}

process predict {
  //conda "${workflow.projectDir}/funannotate-env.yaml"

  input:
  path contigs 
  path "fun/training" 
  path proteins

  output:
  path "fun/{training,predict_misc,predict_results}" includeInputs true

  script:
  options = "--protein_evidence "
  if (proteins) {
    options += proteins.join(" ")
  }
  options += ' $FUNANNOTATE_DB/uniprot_sprot.fasta'
  """
  echo "tes"
  funannotate predict -i $contigs -o fun --species "${params.species}" ${options} --cpus ${params.cpus} 
  """
}

process predict_protein_evidence_only {
  input:
  path contigs 
  path proteins

  output:
  path "fun/{training,predict_misc,predict_results}" includeInputs true

  script:
  options = "--protein_evidence "
  if (proteins) {
    options += proteins.join(" ")
  }
  options += ' $FUNANNOTATE_DB/uniprot_sprot.fasta'
  """
  funannotate predict -i $contigs -o fun --species "${params.species}" ${options} --cpus ${params.cpus} 
  """

}

process annotate_utrs {
  //conda "${workflow.projectDir}/funannotate-env.yaml"

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
  //conda "${workflow.projectDir}/funannotate-env.yaml"

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
  Channel.fromPath(params.contigs).set{ch_contigs}
  if (params.reads)
    Channel.fromPath(params.reads)
         .set{ch_reads}
  else
    Channel.empty()
         .set{ch_reads}
    
  Channel.fromPath(params.proteins)
         .set{ch_proteins}

  ch_contigs = ch_contigs.branch {
    gzipped: it.name.endsWith(".gz")
    uncompressed: !it.name.endsWith(".gz")
  }
  // prepare protein sequences
  ch_proteins = cleanup_proteins(ch_proteins)

  // prepare contigs
  gunzip_contigs(ch_contigs.gzipped) 
  ch_contigs = ch_contigs.uncompressed.concat(gunzip_contigs.out)
  if (! params.no_clean)
    ch_contigs = clean_contigs(ch_contigs)
  ch_contigs = sort_contigs(ch_contigs)
  if (! params.no_repeat_masking)
    ch_contigs = mask_contigs(ch_contigs)

  // prepare reads
  if (params.split_single_fastq)
    ch_reads = split_single_fastq(ch_reads)

  // train if rnaseq data is present
  train(ch_contigs, ch_reads.flatten().toList())

  // predict 
  predict(ch_contigs, 
          train.out.ifEmpty(file("training", type:"dir")), 
          ch_proteins.toList()) | annotate_utrs | annotate_function

}
