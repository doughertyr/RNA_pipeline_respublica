#!/miniconda3/bin/python3

import glob
import os


### PATHS ### 
respublica_path = "/mnt/isilon/cbmi/variome/Taylor/rnaedit/"
ref_genome_dir = "/mnt/isilon/cbmi/variome/rnaseq_workspace/refs/hg19/fasta/"
ref_genome = ref_genome_dir + "hg19.fa"
output_dir = respublica_path + "rnaedit_tmp/"
data_dir = respublica_path + "fastq/"
genomeDir = respublica_path + "genomeDir/"
STAR_output_dir = respublica_path + "STARsnakeOutput/"
GATK = "/mnt/isilon/cbmi/variome/bin/java -D java.io.tmpdir=/mnt/isilon/cbmi/variome/gatk_tmp -jar /mnt/isilon/cbmi/variome/bin/GATK/3.4-46/GenomeAnalysisTK.jar"
dbsnp = respublica_path + "vcf/00-All.vcf"
ref_gene = "ucsc_hg19.txt"
htslib_path = "/mnt/isilon/cbmi/variome/Taylor/htslib"

### TOOLS ###
STAR = "/mnt/isilon/cbmi/variome/rnaseq_workspace/tools/STAR-STAR_2.4.1c/bin/Linux_x86_64_static/STAR"
samtools = "/mnt/isilon/cbmi/variome/bin/samtools"
CreateSequenceDictionary = "/mnt/isilon/cbmi/variome/bin/java -jar /mnt/isilon/cbmi/variome/bin/CreateSequenceDictionary.jar"
AddOrReplaceReadGroups = "/mnt/isilon/cbmi/variome/bin/java -jar /mnt/isilon/cbmi/variome/bin/AddOrReplaceReadGroups.jar"
samtools = "/mnt/isilon/cbmi/variome/bin/samtools"
samtools_index = samtools + " index"
GATK = "/mnt/isilon/cbmi/variome/bin/java -jar /mnt/isilon/cbmi/variome/bin/GATK/3.4-46/GenomeAnalysisTK.jar"
SplitNCigarReads = GATK + " -T SplitNCigarReads"
ReorderSam = "/mnt/isilon/cbmi/variome/bin/java -jar /mnt/isilon/cbmi/variome/bin/ReorderSam.jar"
RealignerTargetCreator = GATK + " -T RealignerTargetCreator"
IndelRealigner = GATK + " -T IndelRealigner"
BaseRecalibrator = GATK + " -T BaseRecalibrator"
PrintReads_BQSR = GATK + " -T PrintReads"
UnifiedGenotyper = GATK + " -T UnifiedGenotyper"
HaplotypeCaller = GATK + " -T HaplotypeCaller"
giremi = "./giremi"

### CHECKPOINT FILES ###
chrName = genomeDir + "chrName.txt" #Created after STAR_index step

data_files = glob.glob("{}/*".format(data_dir))
samples = [os.path.basename(i).strip(".fastq") for i in data_files]

STAR_output = ["{}{}STAR_output.bam".format(STAR_output_dir,i) for i in samples]
AddOrReplaceReadGroups_output = ["{}{}AddOrReplaceReadGroups_output.bam".format(output_dir,i) for i in samples]
samtools_index_1_output = ["{}{}AddOrReplaceReadGroups_output.bam.bai".format(output_dir,i) for i in samples]
SplitNCigarReads_output = ["{}{}SplitNCigarReads_output.bam".format(output_dir,i) for i in samples]
ReorderSam_output = ["{}{}ReorderSam_output.bam".format(output_dir,i) for i in samples]
RealignerTargetCreator_output = ["{}{}RealignerTargetCreator_output.intervals".format(output_dir,i) for i in samples]
IndelRealigner_output = ["{}{}IndelRealigner_output.bam".format(output_dir,i) for i in samples]
BaseRecalibrator_output = ["{}{}BaseRecalibrator_output.table".format(output_dir,i) for i in samples]
PrintReads_BQSR_output = ["{}{}PrintReads_BQSR_output.bam".format(output_dir,i) for i in samples]
HaplotypeCaller_output = ["{}{}HaplotypeCaller_output.vcf".format(output_dir,i) for i in samples]
AnnotateSNVs_output = ["{}AnnotateSNVs_output.txt".format(i) for i in samples]
giremi_output = ["{}giremi_output".format(i) for i in samples]

ref_genome_fai = ref_genome + ".fai"
ref_genome_dict = ref_genome_dir + "hg19.dict"
ref_gene_uniq = "ucsc_hg19_uniq.txt"

### RULES ###


rule all:	
    input: genomeDir, chrName, STAR_output_dir, ref_genome_fai, ref_genome_dict, AddOrReplaceReadGroups_output, samtools_index_1_output, ReorderSam_output, SplitNCigarReads_output#,RealignerTargetCreator_output, IndelRealigner_output, BaseRecalibrator_output, PrintReads_BQSR_output,HaplotypeCaller_output, ref_gene_uniq, AnnotateSNVs_output, giremi_output

rule clean:
    shell: "rm -rf {STAR_output_dir}"

rule make_output_dir:
	shell: "mkdir -p {output_dir}"

rule make_genomeDir:
    output: genomeDir
    shell: "mkdir -p {output}"

rule STAR_index:
    output: chrName
    threads: 12 
    shell: "{STAR} --runMode genomeGenerate --runThreadN 12 --genomeDir {genomeDir} --genomeFastaFiles {ref_genome}" 

rule make_STAR_output_dir:
    output: STAR_output_dir        
    shell: "mkdir -p {output}"


rule STAR_mapping:
    input: "fastq/SRR{tag, \d+}.fastq"
    output: "STARsnakeOutput/SRR{tag, \d+}STAR_output.bam"
    threads: 12
    shell:
        "{STAR} --runThreadN 12 --genomeDir {genomeDir} --readFilesIn {input} --sjdbGTFfile /mnt/isilon/cbmi/variome/rnaseq_workspace/refs/hg19/ncbi.gff --twopassMode Basic --outFileNamePrefix STARsnakeOutput/{wildcards.tag} --outSAMtype BAM SortedByCoordinate"

rule make_ref_genome_fai: 
	output: ref_genome_fai
	shell: "{samtools} faidx {ref_genome}"	

rule make_ref_genome_dict:
	output: ref_genome_dict
	shell: "{CreateSequenceDictionary} R={ref_genome} O= {output}"


rule AddOrReplaceReadGroups:
	input: STAR_mapping_output = STAR_output_dir + "SRR{tag, \d+}STAR_output.bam"
	output: output_dir + "SRR{tag, \d+}AddOrReplaceReadGroups_output.bam"
	threads: 12
	shell: "{AddOrReplaceReadGroups} RGLB={data_dir}{wildcards.tag} RGPL=illumina RGPU=run RGSM={wildcards.tag} I={input.STAR_mapping_output} O={output}"

rule samtools_index_1:
        input: AddOrReplaceReadGroups_output = output_dir + "SRR{tag, \d+}AddOrReplaceReadGroups_output.bam"
        output: output_dir + "SRR{tag, \d+}AddOrReplaceReadGroups_output.bam.bai"
        shell: "{samtools_index} {input.AddOrReplaceReadGroups_output}"

rule ReorderSam:
        input: AddOrReplaceReadGroups_output = output_dir + "SRR{tag, \d+}AddOrReplaceReadGroups_output.bam", samtools_index_1_output = output_dir + "SRR{tag, \d+}AddOrReplaceReadGroups_output.bam.bai"
        output: output_dir + "SRR{tag, \d+}ReorderSam_output.bam"
        shell: "{ReorderSam} I={input.AddOrReplaceReadGroups_output} R={ref_genome} O={output}"

rule samtools_index_2:
        input: ReorderSam_output = output_dir + "SRR{tag, \d+}ReorderSam_output.bam"
        output: output_dir + "SRR{tag, \d+}ReorderSam_output.bam.bai"
        shell: "{samtools_index} {input.ReorderSam_output}"

rule SplitNCigarReads:
	input: ReorderSam_output = output_dir + "SRR{tag, \d+}ReorderSam_output.bam", samtools_index_2_output = output_dir + "SRR{tag, \d+}ReorderSam_output.bam.bai"
	output: output_dir + "SRR{tag, \d+}SplitNCigarReads_output.bam"
	shell: "{SplitNCigarReads} -R {ref_genome} -I {input.ReorderSam_output} --filter_reads_with_N_cigar -o {output}"

#rule IndelRealigner:
#	input: RealignerTargetCreator_output = output_dir + "SRR{tag, \d+}RealignerTargetCreator_output.intervals", ReorderSam_output = output_dir + "SRR{tag, \d+}ReorderSam_output.bam"
#	output: output_dir + "SRR{tag, \d+}IndelRealigner_output.bam"
#	threads: 12
#	shell: "{IndelRealigner} -R {ref_genome} --targetIntervals {input.RealignerTargetCreator_output} --filter_reads_with_N_cigar -I {input.ReorderSam_output} -o {output}"
#
#rule BaseRecalibrator:
#	input: IndelRealigner_output = output_dir + "SRR{tag, \d+}IndelRealigner_output.bam"	
#	output: output_dir + "SRR{tag, \d+}BaseRecalibrator_output.table"
#	threads: 12
#	shell: "{BaseRecalibrator} -R {ref_genome} -I {input.IndelRealigner_output} -knownSites {dbsnp} -o {output}"
#
#rule PrintReads_BQSR:
#	input: BaseRecalibrator_output = output_dir + "SRR{tag, \d+}BaseRecalibrator_output.table", IndelRealigner_output = output_dir + "SRR{tag, \d+}IndelRealigner_output.bam"
#	output: output_dir + "SRR{tag, \d+}PrintReads_BQSR_output.bam"
#	threads: 12
#	shell: "{PrintReads_BQSR} -R {ref_genome} -I {input.IndelRealigner_output} -BQSR {input.BaseRecalibrator_output} -o {output}"
#
#rule HaplotypeCaller:
#	input: PrintReads_BQSR_output = output_dir + "SRR{tag, \d+}PrintReads_BQSR_output.bam"
#	output: output_dir + "SRR{tag, \d+}HaplotypeCaller_output.vcf"
#	shell: "{HaplotypeCaller} -R {ref_genome} -I {input} --dbsnp {dbsnp} -stand_call_conf 20 -stand_emit_conf 20 -drf DuplicateRead -o {output} -mmq 20"
#
#rule UniqGenes:
#	input: {ref_gene}
#	output: {ref_gene_uniq}
#	script: "genename_uniquey.py"
#
#rule AnnotateSNVs:
#	input: HaplotypeCaller_output = output_dir + "SRR{tag, \d+}HaplotypeCaller_output.vcf", UniqGenes_output = {ref_gene_uniq}
#	output: "SRR{tag, \d+}AnnotateSNVs_output.txt"
#	shell: "python annotate_snvs.py {input.HaplotypeCaller_output}"
#
#rule giremi:
#	input: PrintReads_BQSR_output = output_dir + "SRR{tag, \+d}PrintReads_BQSR_output.bam", AnnotateSNVs_output = "SRR{tag, \d+}AnnotateSNVs_output.txt"
#	output: "SRR{tag, \d+}giremi_output"
#	shell: "LD_LIBRARY_PATH={htslib_path}; export LD_LIBRARY_PATH; {giremi} -f {ref_genome} -l {input.AnnotateSNVs_output} -o {output} {input.PrintReads_BQSR_output}"
