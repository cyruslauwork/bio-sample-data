export const configs = [
  {
    hasSEPE: true,
    outputs: [
      {
        fileName: "unsorted_bam.bam",
        id: "unsorted_bam",
        displayName: "Unsorted BAM",
      },
      {
        fileName: "sorted_bam.bam",
        id: "sorted_bam",
        displayName: "Sorted BAM",
      },
      {
        fileName: "sorted_bam.bai",
        id: "sorted_bam_bai",
        displayName: "Sorted BAI",
      },
      {
        fileName: "fasta.fa",
        id: "fasta",
        displayName: "FASTA",
      },
      {
        fileName: "fasta_fai.fai",
        id: "fasta_fai",
        displayName: "FASTA FAI",
      },
    ],
    defaultHardware: {
      cpu: 8,
      ramInGB: 32,
    },
    toolid: "igv",
    displayName: "IGV",
    category: "Visualisation",
    parameters: [
      {
        description: "A Sequence Alignment Map (SAM) file",
        id: "sam",
        type: "file",
        displayName: "SAM",
        validation: {
          required: true,
          acceptedInputTypes: ["sam"],
        },
      },
      {
        description:
          "A FASTA file use to prepare the index reference sequences",
        id: "ref",
        type: "file",
        displayName: "Reference Index",
        validation: {
          required: true,
          acceptedInputTypes: ["fasta"],
        },
      },
      {
        description:
          "Number of threads for parallel search. Threads will run on separate processors/cores and synchronize when parsing reads and outputting alignments. Searching for alignments is highly parallel, and speedup is close to linear. Increasing -p increases Bowtie 2's memory footprint. E.g. when aligning to a human genome index, increasing -p from 1 to 8 increases the memory footprint by a few hundred megabytes.",
        id: "--threads",
        type: "number",
        displayName: "Threads",
        defaultValue: 8,
      },
    ],
  },
];
