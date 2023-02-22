conditions = glob_wildcards("sequences/{condition}_R1_001.fastq.gz").condition

rule all:
    input:
        expand("trimmed_reads/{cond}_{read}_trimmed.fastq", cond = conditions, read = ["R1", "R2"]),
        expand("trimmed_reads/{cond}_{read}_trimmed_un.fastq", cond = conditions, read = ["R1", "R2"]),
        expand("assemblies/{cond}/{cond}.contigs.fa", cond = conditions)

rule trimmomatic:
    input:
        F = "sequences/{cond}_R1_001.fastq.gz",
        R = "sequences/{cond}_R2_001.fastq.gz",
        Adapter = "scripts/adapters.fasta"
    output:
        "trimmed_reads/{cond}_R1_trimmed.fastq",
        "trimmed_reads/{cond}_R1_trimmed_un.fastq",
        "trimmed_reads/{cond}_R2_trimmed.fastq",
        "trimmed_reads/{cond}_R2_trimmed_un.fastq"
    conda:
        "conda/trimmomatic.yaml"
    shell:
        "trimmomatic PE -threads 8  {input.F} {input.R} "
        "{output} ILLUMINACLIP:{input.Adapter}:2:30:10 MINLEN:40"

rule assembly:
    input:
        F = "trimmed_reads/{cond}_R1_trimmed.fastq",
        R = "trimmed_reads/{cond}_R2_trimmed.fastq"
    params:
        outdir = "assemblies/{cond}",
        prefix = "{cond}"
    output:
        "assemblies/{cond}/{cond}.contigs.fa"
    conda:
        "conda/megahit.yaml"
    shell:
        "megahit -1 {input.F} -2 {input.R} -f -o {params.outdir} --out-prefix " 
        "{params.prefix}"
        #"{params.prefix} --num-cpu-threads 6 -m 0.5 --mem-flag 1 --tmp-dir /tmp"