process test {
//   conda.enabled = true
//   conda = 'samtools'

  cpus = '1'
  memory = '1G'
  queue = 'nextflow-ci'
  
//   container "docker.img"
//   containerOptions ""
   
  input:
  path 'input.sam'
//   path ('input.sam', arity: '1')         // exactly one file is expected
//   path input2

  output:
//   path 'results/sorted_output.bam', emit: output1
//   path 'results/sorted_output.bam.bai', emit: output2

  script:
  """
  mkdir -p results
//   samtools view -S -b -o results/output.bam input.sam
//   samtools sort -o results/sorted_output.bam results/output.bam
//   samtools index results/sorted_output.bam
  """
}

workflow {
  def sam = Channel.fromPath( '*.sam' )
  thisWorkflow = test(sam)
//   thisWorkflow.view { it }
}