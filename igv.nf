process sam_to_bam_and_fasta_to_fai {
  cpus = '1'
  memory = '1G'
  queue = 'nextflow-ci'
   
  input:
  path input1

  output:
  path 'results/unsorted_bam.bam', emit: unsorted_bam
  path 'results/sorted_bam.bam', emit: sorted_bam
  path 'results/sorted_bam.bam.bai', emit: sorted_bam_bai

  script:
  """
  mkdir -p results
  samtools view -S -b -o results/unsorted_bam.bam $input1
  samtools sort -o results/sorted_bam.bam results/unsorted_bam.bam
  samtools index results/sorted_bam.bam
  mv results/sorted_bam.bam.bai results/sorted_bam_index.bai
  """
}

workflow {
  def sam = Channel.fromPath( '*.sam' )
  sam_to_bam(sam)
}