process {{id}} {
  cpus "${params.{{id}}Cpu}"
  memory "${params.{{id}}Memory}"
  queue "${params.{{id}}Queue}"

  errorStrategy { sleep(Math.pow(3, task.attempt) * 1000 as long); return 'retry' }
  maxRetries 3
  publishDir "${params.outputBucket}/${params.userid}/${params.jobid}/{{id}}"
  resourceLabels jobID: '{{jobid}}', userID: '{{userid}}', app: '{{apptag}}', service_tagged: 'batch', stage: "${params.stage}", processID: '{{id}}'
  container "public.ecr.aws/e9z9p5l1/biotailer-igv-amd64"
  containerOptions ""
   
  {{inputs}}

  output:
  path 'results/unsorted_bam.bam', emit: unsorted_bam
  path 'results/sorted_bam.bam', emit: sorted_bam
  path 'results/sorted_bam_bai.bai', emit: sorted_bam_bai
  path 'results/fasta.fa', emit: fasta
  path 'results/fasta_fai.fai', emit: fasta_fai

  script:
  tag_script = """
    # Obtain the IMDSv2 token
    TOKEN=\\$(wget --header="X-aws-ec2-metadata-token-ttl-seconds: 21600" --method=PUT -qO- "http://169.254.169.254/latest/api/token")

    # Use the token to fetch the instance ID
    INSTANCE_ID=\\$(wget --header="X-aws-ec2-metadata-token: \\$TOKEN" -qO- "http://169.254.169.254/latest/meta-data/instance-id")

    # Retrieve the EBS volume IDs
    VOLUME_IDS=\\$(aws ec2 describe-instances --instance-ids \\$INSTANCE_ID --query 'Reservations[].Instances[].BlockDeviceMappings[].Ebs.VolumeId' --region us-east-1 --output text)

    # Tag the EBS volumes
    if [ -n \\$VOLUME_IDS ]; then
        for VOLUME_ID in \\$VOLUME_IDS; do
            aws ec2 create-tags --resources \\$VOLUME_ID --tags Key=userID,Value=${params.userid} Key=jobID,Value=${params.jobid} Key=service_tagged,Value=ebs Key=stage,Value=${params.stage} --region us-east-1
        done
    else
        echo "Error: No EBS volume found for instance ID \\$INSTANCE_ID"
    fi
    """
  return ("mkdir -p results && " +
    tag_script +
  """
  samtools view -S -b -o results/unsorted_bam.bam ${{sam}} &&
  samtools sort -o results/sorted_bam.bam results/unsorted_bam.bam &&
  samtools index results/sorted_bam.bam &&
  mv results/sorted_bam.bam.bai results/sorted_bam.bai
  samtools faidx ${{ref}}
  """
  )
}