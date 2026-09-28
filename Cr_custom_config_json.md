## Custom config
```bash
process {

stageInMode = 'symlink'
clusterOptions = '--account=ga03714 --partition=genoa,milan'

    withName: 'BWAMEM2_MEM' {
        memory = 100.GB
   }
   withName: 'BWAMEM2_INDEX' {
        memory = 50.GB
   }
   withName: 'MANTA_GERMLINE' {
        time = 72.h
   }
   withName: 'GATK4_HAPLOTYPECALLER' {
        time = 72.h
   }
    withName: 'NFCORE_SAREK:SAREK:CRAM_QC_NO_MD:SAMTOOLS_STATS' {
        time = 48.h
   }
    withName: '.*:VCF_VARIANT_FILTERING_GATK:FILTERVARIANTTRANCHES' {
        time = 48.h
   }
    withName: 'NFCORE_SAREK:SAREK:FASTQC' {
        time = 48.h
    }

    withName: 'NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_GERMLINE_ALL:BAM_JOINT_CALLING_GERMLINE_GATK:GATK4_GENOTYPEGVCFS' {
        time = 180.h
        memory = 450.GB
        cpus = 8
        ext.args = {'--intervals /nesi/nobackup/ga03714/Cr_raw_lcWGS_data/Cr_all_sarek/Cr_all_sarek_out/reference/intervals/fai.bed'}
    }

    withName: 'NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_GERMLINE_ALL:BAM_JOINT_CALLING_GERMLINE_GATK:GATK4_GENOMICSDBIMPORT'{ 
        time = 72.h
        memory = 60.GB
    }

}

apptainer {
  pullTimeout = '4h'
}

errorStrategy = { task.exitStatus == 140 ? 'retry' : 'terminate' }
maxRetries = 5

executor {
    name = 'slurm'
    array = 100
}

trace {
    enabled = true 
}
cleanup = true
```

## Params JSON
```bash
input"/nesi/nobackup/ga03714/Cr_raw_lcWGS_data/Cr_all_sarek/Cr_all_samplesheet_sarek.csv"
outdir"/nesi/nobackup/ga03714/Cr_raw_lcWGS_data/Cr_all_sarek/Cr_all_sarek_out"
outdir_cache"/nesi/nobackup/ga03714/Cr_raw_lcWGS_data/Cr_all_sarek/Cr_all_sarek_out"
tools"haplotypecaller"
skip_tools"baserecalibrator"
trim_fastqtrue
save_mappedtrue
save_output_as_bamtrue
save_referencetrue
genomenull
igenomes_ignoretrue
fasta"/nesi/nobackup/ga03714/Cr_raw_lcWGS_data//Crod_v1.1_h1_masked.fa"
fasta_fai"/nesi/nobackup/ga03714/Cr_raw_lcWGS_data/Crod_v1.1_h1_masked.fa.fai"
email"mneh623@aucklanduni.ac.nz"
multiqc_title"Cr_all"
joint_germlinetrue
save_trimmedtrue
trim_nextseqtrue
detect_adapter_for_petrue
```
