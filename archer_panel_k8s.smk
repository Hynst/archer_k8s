### Archer panel pipeline for kubernetes ###
######
######
#

import re
import datetime
import boto3
from snakemake.utils import min_version
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

configfile: "./config/config_211204_A42.yaml"
#container: "cerit.io/snakemake:v6.10.0"

# check minimal version
min_version("5.18.0")

# S3 credentials
AWS_ID="bibs"
AWS_KEY="oIhPZYg3xHiE206N"
S3_BUCKET = "bibs"

# conect to remote S3
S3 = S3RemoteProvider(host="https://storage-elixir1.cerit-sc.cz", access_key_id=AWS_ID, secret_access_key=AWS_KEY)
#boto3.client('s3', aws_access_key_id=AWS_ID, aws_secret_access_key=AWS_KEY, region_name="US", endpoint_url="https://storage-elixir1.cerit-sc.cz")

#Global variables
REF_SEQ = S3.remote(S3_BUCKET + "/references/GRCh37-p13/seq/GRCh37-p13.fa")
REF_BWA = S3.remote(S3_BUCKET + "/references/GRCh37-p13/index/BWA/GRCh37-p13")

TARGET_BED = S3.remote(S3_BUCKET + "/src/myelo_proliferation/archer/beds/jana_archer_unique_plus2nt.bed")
TARGET_INTERVALS = S3.remote(S3_BUCKET + "/src/myelo_proliferation/archer/beds/jana_archer_unique_plus2nt.intervals")


CUTADAPT = "cutadapt"
BWA = "bwa"
SAMBLASTER= "samblaster"
SAMTOOLS= "samtools"
SAMBAMBA= "sambamba"
PICARD = S3.remote(S3_BUCKET + "/src/myelo_proliferation/archer/picard.jar")
UMI_TOOLS = "umi_tools"
MULTIQC = "multiqc"

GATK2= S3.remote(S3_BUCKET + "/src/myelo_proliferation/archer/GATK/gatk")

VEP_DATA = S3.remote(S3_BUCKET + "/references/GRCh37-p13/VEP/")
VEP= "vep"
BEDTOOLS= "bedtools"
VCF2TABLE= S3.remote(S3_BUCKET + "/src/myelo_proliferation/archer/reGenotype/vcf2table_MDS.py")
BCFTOOLS = "bcftools"
FLT3= S3.remote(S3_BUCKET + "/src/myelo_proliferation/archer/FLT3_ITD/FLT3_ITD_ext/FLT3_ITD_ext.pl")

# post processing
#Cover_Stat_noProc= "/mnt/nfs/shared/MedGen/MDS_panel_jana/src/coverage_stat_jana_noprocessed.R"
Cover_Stat = S3.remote(S3_BUCKET + "/src/myelo_proliferation/archer/scripts/coverage_stat_jana.R")

#COV_STAT = "/mnt/nfs/shared/999993-Bioda/projects/lymphopanel/src/coverage_lympho_stat.R"

def flagstatlist_second(wildcards):
    statlist = []
    run = config["Run"]
    samples = config["Samples"]
    for sample in samples:
            statlist = statlist + S3.remote([S3_BUCKET + "/data/myelo/" + run + "/samples/" + sample \
                       + "/stats/second_bam_qc/" + sample + ".flagstat"])
    return statlist

def covlist(wildcards):
    covlist = []
    run = config["Run"]
    samples = config["Samples"]
    for sample in samples:
            covlist = covlist + S3.remote([S3_BUCKET + "/data/myelo/" + run + "/samples/" + sample + "/coverage/" + sample \
            + ".perexon_stat.txt"])
    return covlist

def vcfvep_list(wildcards):
        vcfvep_list = []
        run = config["Run"]
        samples = config["Samples"]
        for sample in samples:
            vcfvep_list = vcfvep_list + S3.remote([S3_BUCKET + "/data/myelo/" + run + "/samples/" + sample + "/variants/mutect2/" + sample \
            + ".mutect2.filt.norm.vep.csv"])
        return vcfvep_list

def FLT3_list(wildcards):
        FLT3_list = []
        run = config["Run"]
        samples = config["Samples"]
        for sample in samples:
            FLT3_list = FLT3_list + S3.remote([S3_BUCKET + "/data/myelo/" + run + "/samples/" + sample + "/variants/flt3/" + sample \
            + ".second_FLT3_ITD_summary.txt"])
        return FLT3_list

def getFasta(wildcards):
    files = {'fwd': S3.remote(S3_BUCKET + config["Samples"][wildcards.sample].get('fwd'))}
    files['rev'] = S3.remote(S3_BUCKET + config["Samples"][wildcards.sample].get('rev'))
    return files

def multiqc(wildcards):
    run = config['Run']
    multiqc_f = S3.remote([S3_BUCKET + "/data/myelo/" + run + "/multiqc_reports/multiqc.second"])
    return multiqc_f

rule all:
    output:
        S3.remote(S3_BUCKET + "/src/myelo_proliferation/run_" + str(datetime.datetime.now()).replace(" ","-"))
    input:
        vcfvep_list,
        #multiqc,
        covlist,
        FLT3_list,
        S3.remote(S3_BUCKET + "/data/myelo/" + config["Run"] + "/" + config["Run"] + ".merged_variant.table_NEW.tsv")
    shell:
        """
        echo {input}
        touch {output}
        """

rule merge_tables:
     input:
         csv = S3.remote(expand("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.filt.norm.vep.csv", S3_BUCKET = "bibs", run = config["Run"], sample = config["Samples"])),
         mrg = S3.remote("{S3_BUCKET}/src/myelo_proliferation/archer/scripts/merge_tables.R")
     output:
         S3.remote("{S3_BUCKET}/data/myelo/{run}/{run}.merged_variant.table_NEW.tsv")
     conda: "./envs/erko.yml"
     shell:
         "Rscript --vanilla {input.mrg} {wildcards.run} {S3_BUCKET}/data/myelo/"
#         command = "Rscript --vanilla ./merge_tables.R " + config['Run']
#         shell(command)

rule coverage_cons_postprocess:
    input:
        txt = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/coverage/{sample}.PBcov.cons.txt"),
        rsc = S3.remote("{S3_BUCKET}/src/myelo_proliferation/archer/scripts/coverage_stat_jana.R")
    output:
        S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/coverage/{sample}.perexon_stat.txt")
    conda: "./envs/erko.yml"
    shell:
        "Rscript --vanilla {input.rsc} {input.txt}"
        # command = "/mnt/ssd/ssd_2/install/dir/anaconda/envs/myelopanel/bin/Rscript --vanilla " + Cover_Stat + " " + str(input)
        # print(command)
        # shell(command)
        # command2 = "cp " + str(output.cov1) + " " + str(output.cov2)
        # shell(command2)


rule coverage_cons:
    input:
        bam = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/bam/{sample}.second.bam"),
        bed = S3.remote("{S3_BUCKET}/src/myelo_proliferation/archer/beds/jana_archer_unique_plus2nt.bed")
    output:
        S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/coverage/{sample}.PBcov.cons.txt")
    conda: "./envs/bedtools.yml"
    shell:
        "{BEDTOOLS} coverage -abam {input.bed} -b {input.bam} -d > {output}"
        # command = "bedtools coverage -abam {TARGET_BED}" \
        #     + " -b {input} -d > {output}"
        # shell(command)

## mutation branch
rule FLT3_analysis:
    output:
        S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/variants/flt3/{sample}.second_FLT3_ITD_summary.txt")
    input:
        bam = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/bam/{sample}.second.bam"),
        bai = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/bam/{sample}.second.bam.bai"),
        flt = S3.remote("{S3_BUCKET}/src/myelo_proliferation/archer/FLT3_ITD/tmp/FLT3.tar.gz"),
        jar = S3.remote("{S3_BUCKET}/src/myelo_proliferation/archer/picard.jar"),
        fgb = S3.remote("{S3_BUCKET}/src/myelo_proliferation/archer/fgbio-0.8.1.jar")
    conda: "./envs/perl.yml"
    shell:
        """
        tar -C {S3_BUCKET}/src/myelo_proliferation/archer/FLT3_ITD/tmp/ -xf {input.flt}
        ls -al bibs/src/myelo_proliferation/archer/FLT3_ITD/tmp/FLT3_ITD_ext/
        perl {S3_BUCKET}/src/myelo_proliferation/archer/FLT3_ITD/tmp/FLT3_ITD_ext/FLT3_ITD_ext.pl --bam {input.bam} --output {S3_BUCKET}/data/myelo/{wildcards.run}/samples/{wildcards.sample}/variants/flt3/ --ngstype amplicon --genome hg19
        """
        # id = wildcards.sample
        # run = wildcards.run
        # out_fold = "../data/" + run + "/samples/" + id + "/variants/flt3/"
        # command = FLT3 + " --bam " + str(input) + " --output " \
        # + str(out_fold) + " --ngstype amplicon --genome hg19"
        # #command2 = "cp " + str(output) + " ../../Myelo/runs/{run}/variants/" + id + ".second_FLT3_ITD_summary.txt"
        # shell(command)
        # #shell(command2)

rule vcf2csv:
    output:
        csv = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.filt.norm.vep.csv")
    input:
        vcf = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.filt.norm.vep.vcf"),
        v2t = S3.remote("{S3_BUCKET}/src/myelo_proliferation/archer/reGenotype/vcf2table_MDS.py")
    conda: "./envs/vcf2csv.yml"
    shell:
        "{BCFTOOLS} view -f 'PASS,clustered_events,multiallelic' {input.vcf} | python {input.v2t} simple --build GRCh37 -i /dev/stdin -o {output.csv}"


rule annotate_mutect_cons:
    output:
        vcf = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.filt.norm.vep.vcf")
    input:
        vcf = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.filt.norm.vcf"),
        ref = S3.remote("{S3_BUCKET}/references/GRCh37-p13/seq/GRCh37-p13.fa"),
        vep_data = S3.remote("{S3_BUCKET}/references/GRCh37-p13/test/VEP_95.tar.gz"),
    conda: "./envs/vep.yml"
    shell:
        """
        tar -C {S3_BUCKET}/references/GRCh37-p13/test/ -xf {input.vep_data}
        {VEP} -i {input.vcf} --cache --cache_version 95 --dir_cache {S3_BUCKET}/references/GRCh37-p13/test/VEP_95 -fasta {input.ref} --merged --offline --vcf --everything -o {output.vcf}
        """
        # id = wildcards.sample
        # command1 = "{VEP} -i {input}" \
        #     + " --cache --cache_version 95 --dir_cache {VEP_DATA}" \
        #     + " --fasta {REF_SEQ}" \
        #     + " --merged --offline --vcf --everything" \
        #     + " -o {output.vcf}"
        # print(command1)
        # # clustered events zatim ponechano, protoze u Jirky mezi nimi byly true positive
        # # eventualne vypnout
        # command2 = "{BCFTOOLS} view -f 'PASS,clustered_events,multiallelic' {output.vcf} " \
        #      + "| python {VCF2TABLE} simple --build GRCh37 -i /dev/stdin -o {output.csv}"
        #
        # command3 = "cp " + str(output.csv) + " " + str(output.csv2)
        # shell(command1)
        # shell(command2)
        # shell(command3)

rule normalize_vcf_mutect_cons:
    output:
        S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.filt.norm.vcf")
    input:
        vcf = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.filt.vcf"),
        ref = S3.remote("{S3_BUCKET}/references/GRCh37-p13/seq/GRCh37-p13.fa")
    conda: "./envs/bcftools.yml"
    shell:
        "{BCFTOOLS} norm -f {input.ref} {input.vcf} -o {output}"
        # command = BCFTOOLS + " norm -f {REF_SEQ} {input} \
        #      -o {output}"
        # shell(command)

rule filter_mutect:
     output:
         S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.filt.vcf")
     input:
         vcf = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.cons.vcf"),
     conda: "./envs/gatk4.yml"
     shell:
         "gatk FilterMutectCalls -V {input.vcf} -O {output}"
        #  command = GATK2 + " FilterMutectCalls" \
        #      + " -V {input}" \
        #      + " -O {output}"
        #  shell(command)

rule mutect2_cons:
    output:
        S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.cons.vcf")
    input:
        bam = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/bam/{sample}.second.bam"),
        ref = S3.remote("{S3_BUCKET}/references/GRCh37-p13/seq/GRCh37-p13.fa"),
        fai = S3.remote("{S3_BUCKET}/references/GRCh37-p13/seq/GRCh37-p13.fa.fai"),
        dct = S3.remote("{S3_BUCKET}/references/GRCh37-p13/seq/GRCh37-p13.dict"),
        ivl = S3.remote("{S3_BUCKET}/src/myelo_proliferation/archer/beds/jana_archer_unique_plus2nt.intervals"),
        bai = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/bam/{sample}.second.bam.bai")
    conda: "./envs/gatk4.yml"
    shell:
        "gatk Mutect2 --reference {input.ref} --input {input.bam} --tumor-sample {wildcards.sample} --annotation StrandArtifact --min-base-quality-score 30 --intervals {input.ivl} --output {output}"
        # id = wildcards.sample
        # command = GATK2 + " Mutect2 --reference {REF_SEQ}" \
        #     + " --input {input}" \
        #     + " --tumor-sample " + id \
        #     + " --annotation StrandArtifact" \
        #     + " --min-base-quality-score 30" \
        #     + " --intervals {TARGET_BED} --output {output}"
        # shell(command)

## QC
# rule multiqc:
#     output:
#         report_second = S3.remote("{S3_BUCKET}/data/myelo/{run}/multiqc_reports/multiqc.second")
#     input:
#         flagstatlist_second
#     conda: "./envs/multiqc.yml"
#     shell:
#         "{MULTIQC} {S3_BUCKET}/data/myelo/{wildcards.run}/samples/*/stats/second_bam_qc -o {output.report_second}"
#         #command_second = "multiqc ../data/" + wildcards.run + "/samples/*/stats/second_bam_qc -o {output.report_second}"
#         #shell(command_second)

rule sedond_bam_QC:
    output:
        flagstat    = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/stats/second_bam_qc/{sample}.flagstat"),
        samstats    = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/stats/second_bam_qc/{sample}.samstats"),
        hs_metrics  = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/stats/second_bam_qc/{sample}.hs_metrics"),
        aln_metrics = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/stats/second_bam_qc/{sample}.aln_metrics")
    input:
        bam = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/bam/{sample}.second.bam"),
        jar = S3.remote("{S3_BUCKET}/src/myelo_proliferation/archer/picard.jar"),
        ivl = S3.remote("{S3_BUCKET}/src/myelo_proliferation/archer/beds/jana_archer_unique_plus2nt.intervallist"),
        ref = S3.remote("{S3_BUCKET}/references/GRCh37-p13/seq/GRCh37-p13.fa")
    conda: "./envs/samtools.yml"
    shell:
        """
        {SAMTOOLS} flagstat {input.bam} > {output.flagstat}
        {SAMTOOLS} stats {input.bam} > {output.samstats}
        java -jar {input.jar} CollectHsMetrics I={input.bam} BAIT_INTERVALS={input.ivl} TARGET_INTERVALS={input.ivl} R={input.ref} O={output.hs_metrics}"
        java -jar {input.jar} CollectAlignmentSummaryMetrics I={input.bam} R={input.ref} O={output.hs_metrics}
        """
        # command_flag = SAMTOOLS + " flagstat {input} > {output.flagstat}"
        # command_samtools = SAMTOOLS +  " stats {input.bam} > {output.samstats}"
        # command_hs_met   = "java -Djava.io.tmpdir={TEMP_DIR} -jar " + PICARD + \
        #     " CollectHsMetrics I={input.bam}"       + \
        #     " BAIT_INTERVALS="   + TARGET_INTERVALS + \
        #     " TARGET_INTERVALS=" + TARGET_INTERVALS + \
        #     " R=" + REF_SEQ + " O={output.hs_metrics}"
        # command_aln_met  = "java -Djava.io.tmpdir={TEMP_DIR} -jar {PICARD}" + \
        #     " CollectAlignmentSummaryMetrics  I={input.bam}" + \
        #     " R={REF_SEQ} O={output.aln_metrics}"
        # shell(command_flag)
        # shell(command_samtools)
        # shell(command_hs_met)
        # shell(command_aln_met)

rule dedup:
    output:
        bam = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/bam/{sample}.second.bam"),
        bai = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/bam/{sample}.second.bam.bai")
    input:
        bam = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/bam/{sample}.first.bam"),
        bai = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/bam/{sample}.first.bam.bai")
    conda: "./envs/umi_tools.yml"
    log:
        stats = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/bam/{sample}.dedup")
    shell:
           """
           {UMI_TOOLS} dedup -I {input} --paired --output-stats={log.stats} -S {output}
           {SAMTOOLS} index {output.bam}
           """
    # run:
    #     command = "umi_tools dedup -I " + str(input) + " --paired --output-stats=" + str(log.stats) \
    #         + " -S " + str(output)
    #     command2 = SAMTOOLS + " index {output}"
    #     shell(command)
    #     shell(command2)

rule first_align_bam:
    output:
        bam = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/bam/{sample}.first.bam"),
        bai = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/bam/{sample}.first.bam.bai")
    input:
        fwd = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/fastq_trimmed/{sample}.trimmed.UMI.R1.fastq"),
        rev = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/fastq_trimmed/{sample}.trimmed.UMI.R2.fastq"),
        ref_bwa = S3.remote(expand(S3_BUCKET + "/references/GRCh37-p13/index/BWA/GRCh37-p13{prip}", prip = ["", ".amb",".ann",".bwt",".pac",".sa"]))
    conda: "./envs/bwa.yml"
    log:
        run = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/bam/{sample}.log")
    shell:
           """
           {BWA} mem -t 20 -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:illumina' -v 1 {input.ref_bwa[0]} {input.fwd} {input.rev} 2>> {log.run} | {SAMTOOLS} view -Sb - | {SAMBAMBA} sort -o {output.bam} /dev/stdin 2>> {log.run}
           {SAMTOOLS} index {output.bam}
           """
#    run:
#        id = wildcards.sample
#        sm = wildcards.sample
#        command = BWA + " mem -t 20 -R '@RG\\tID:" + id + "\\tSM:" + sm + "\\tPL:illumina' \
#            -v 1 {REF_BWA} {input.fwd} {input.rev} 2>> {log.run}\
#            | {SAMTOOLS} view -Sb - | {SAMBAMBA} sort -o {output} /dev/stdin 2>> {log.run}"
#        command2 = SAMTOOLS + " index " + str(output)
#        shell(command)
#        shell(command2)

### UMI a read preprocessing ###
# adaptors to trim: AACCGCCAGGAGT
# UMI 1-8 bases in read position
rule trim_adaptors:
    output:
        fwd = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/fastq_trimmed/{sample}.trimmed.UMI.R1.fastq"),
        rev = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/fastq_trimmed/{sample}.trimmed.UMI.R2.fastq")
    input:
        fwd = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/fastq_trimmed/{sample}.trimmed.R1.fastq"),
        rev = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/fastq_trimmed/{sample}.trimmed.R2.fastq")
    conda: "./envs/cutadapt.yml"
    shell: "{CUTADAPT} -g AACCGCCAGGAGT -m 50 -o {output.fwd} -p {output.rev} {input.fwd} {input.rev} > trim.out"
#    run:
#        command = CUTADAPT + " -g AACCGCCAGGAGT -m 50 -o " + str(output.fwd) \
#            + " -p " + str(output.rev) + " " + str(input.fwd) + " " + str(input.rev) \
#            + " > trim.out"
#        shell(command)

rule UMI_extract:
    output:
        fwd = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/fastq_trimmed/{sample}.trimmed.R1.fastq"),
        rev = S3.remote("{S3_BUCKET}/data/myelo/{run}/samples/{sample}/fastq_trimmed/{sample}.trimmed.R2.fastq")
    input:
        unpack(getFasta)
    conda: "./envs/umi_tools.yml"
    shell: "{UMI_TOOLS} extract -I {input.fwd} --bc-pattern=NNNNNNNN --read2-in={input.rev} --stdout={output.fwd} --read2-out={output.rev}"
#    run:
#        command = "{UMI_TOOLS} extract -I " + str(input.fwd) \
#            + " --bc-pattern=NNNNNNNN --read2-in=" + str(input.rev) \
#            + " --stdout=" + str(output.fwd) + " --read2-out=" + str(output.rev)
#        shell(command)
