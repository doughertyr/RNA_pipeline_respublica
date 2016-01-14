#!/miniconda3/bin/python3

import glob
import os


### PATHS ### 
ref_genome_dir = "mnt/isilon/cbmi/variome/rnaseq_workspace/refs/hg19/fasta/"
ref_genome = ref_genome_dir + "hg19.fa"
output_dir = "rnaedit_tmp/"
data_dir = "fastq/"
genomeDir = "genomeDir/"
STAR_output_dir = "STARsnakeOutput/"
GATK = "java -jar mnt/isilon/cbmi/variome/bin/GATK/3.4-46/GenomeAnalysisTK.jar"

### TOOLS ###
STAR = "mnt/isilon/cbmi/variome/rnaseq_workspace/tools/STAR-STAR_2.4.1c/bin/Linux_x86_64_static/STAR"
samtools = "mnt/isilon/cbmi/variome/bin/samtools"
CreateSequenceDictionary = "java -jar mnt/isilon/cbmi/variome/bin/CreateSequenceDictionary.jar"
AddOrReplaceReadGroups = "java -jar mnt/isilon/cbmi/variome/bin/AddOrReplaceReadGroups.jar"
GATK = "java -jar mnt/isilon/cbmi/variome/bin/GATK/3.4-46/GenomeAnalysisTK.jar"
ReorderSam = "java -jar mnt/isilon/cbmi/variome/bin/ReorderSam.jar"
samtools = "mnt/isilon/cbmi/variome/bin/samtools"
samtools_index = samtools + " index"
RealignerTargetCreator = GATK + " -T RealignerTargetCreator"
IndelRealigner = GATK + " -T IndelRealigner"
BaseRecalibrator = GATK + " -T BaseRecalibrator"
PrintReads_BQSR = GATK + " -T PrintReads"
UnifiedGenotyper = GATK + " -T UnifiedGenotyper"
HaplotypeCaller = GATK + " -T HaplotypeCaller"
giremi = "mnt/isilon/cbmi/variome/Taylor/giremi/giremi"

### CHECKPOINT FILES ###
chrName = genomeDir + "chrName.txt" #Created after STAR_index step

data_files = glob.glob("{}/*".format(data_dir))
samples = [os.path.basename(i).strip(".fastq") for i in data_files]

STAR_output = ["{}{}STAR_output.bam".format(STAR_output_dir,i) for i in samples]
AddOrReplaceReadGroups_output = ["{}{}AddOrReplaceReadGroups_output.bam".format(output_dir,i) for i in samples]
samtools_index_1_output = ["{}{}AddOrReplaceReadGroups_output.bam.bai".format(output_dir,i) for i in samples]
ReorderSam_output = ["{}{}ReorderSam_output.bam".format(output_dir,i) for i in samples]
RealignerTargetCreator_output = ["{}{}RealignerTargetCreator_output.intervals".format(output_dir,i) for i in samples]
IndelRealigner_output = ["{}{}IndelRealigner_output.bam".format(output_dir,i) for i in samples]
BaseRecalibrator_output = ["{}{}BaseRecalibrator_output.table".format(output_dir,i) for i in samples]
PrintReads_BQSR_output = ["{}{}PrintReads_BQSR_output.bam".format(output_dir,i) for i in samples]
UnifiedGenotyper_output = ["{}{}UnifiedGenotyper_output.vcf".format(output_dir,i) for i in samples]
HaplotypeCaller_output = ["{}{}HaplotyperCaller_output.vcf".format(output_dir,i) for i in samples]
giremi_output = ["{}{}giremi_output.res".format(output_dir,i) for i in samples]

ref_genome_fai = ref_genome + ".fai"
ref_genome_dict = ref_genome_dir + "hg19.dict"


### RULES ###


rule all:	
    input: genomeDir, chrName, STAR_output_dir, ref_genome_fai, ref_genome_dict, AddOrReplaceReadGroups_output, samtools_index_1_output, ReorderSam_output, RealignerTargetCreator_output, IndelRealigner_output, BaseRecalibrator_output, PrintReads_BQSR_output, UnifiedGenotyper_output, HaplotypeCaller_output#, giremi_output

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
        "{STAR} --runThreadN 12 --genomeDir {genomeDir} --readFilesIn {input} --sjdbGTFfile mnt/isilon/cbmi/variome/rnaseq_workspace/refs/hg19/ncbi.gff --twopassMode Basic --outFileNamePrefix STARsnakeOutput/{wildcards.tag} --outSAMtype BAM SortedByCoordinate"

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
	threads: 12
        shell: "{samtools_index} {input.AddOrReplaceReadGroups_output}"

rule ReorderSam:
	input: AddOrReplaceReadGroups_output = output_dir + "SRR{tag, \d+}AddOrReplaceReadGroups_output.bam", samtools_index_1_output = output_dir + "SRR{tag, \d+}AddOrReplaceReadGroups_output.bam.bai"
	output: output_dir + "SRR{tag, \d+}ReorderSam_output.bam"
	threads: 12
	shell: "{ReorderSam} I={input.AddOrReplaceReadGroups_output} R={ref_genome} O={output}"

rule samtools_index_2:
	input: ReorderSam_output = output_dir + "SRR{tag, \d+}ReorderSam_output.bam"
	output: output_dir + "SRR{tag, \d+}ReorderSam_output.bam.bai"
	threads: 12
	shell: "{samtools_index} {input.ReorderSam_output}"

rule RealignerTargetCreator:
	input: ReorderSam_output = output_dir + "SRR{tag, \d+}ReorderSam_output.bam", samtools_index_2_output = output_dir + "SRR{tag, \d+}ReorderSam_output.bam.bai"
	output: output_dir + "SRR{tag, \d+}RealignerTargetCreator_output.intervals"	
	threads: 12
	shell: "{RealignerTargetCreator} -R {ref_genome} -I {input.ReorderSam_output} --filter_reads_with_N_cigar -o {output}"

rule IndelRealigner:
	input: RealignerTargetCreator_output = output_dir + "SRR{tag, \d+}RealignerTargetCreator_output.intervals", ReorderSam_output = output_dir + "SRR{tag, \d+}ReorderSam_output.bam"
	output: output_dir + "SRR{tag, \d+}IndelRealigner_output.bam"
	threads: 12
	shell: "{IndelRealigner} -R {ref_genome} --targetIntervals {input.RealignerTargetCreator_output} --filter_reads_with_N_cigar -I {input.ReorderSam_output} -o {output}"

rule BaseRecalibrator:
	input: IndelRealigner_output = output_dir + "SRR{tag, \d+}IndelRealigner_output.bam"	
	output: output_dir + "SRR{tag, \d+}BaseRecalibrator_output.table"
	threads: 12
	shell: "{BaseRecalibrator} -R {ref_genome} -I {input.IndelRealigner_output} -knownSites vcf/00-All.vcf -o {output}"

rule PrintReads_BQSR:
	input: BaseRecalibrator_output = output_dir + "SRR{tag, \d+}BaseRecalibrator_output.table", IndelRealigner_output = output_dir + "SRR{tag, \d+}IndelRealigner_output.bam"
	output: output_dir + "SRR{tag, \d+}PrintReads_BQSR_output.bam"
	threads: 12
	shell: "{PrintReads_BQSR} -R {ref_genome} -I {input.IndelRealigner_output} -BQSR {input.BaseRecalibrator_output} -o {output}"

rule HaplotypeCaller:
	input: PrintReads_BQSR_output = output_dir + "SRR{tag, +d}PrintReads_BQSR_output.bam"
	output: output_dir + "SRR{tag, +d}HaplotypeCaller_output.vcf"
	threads: 12
	shell: "{HaplotypeCaller} -R {ref_genome} -I {input} --dbsnp vcf/00-All.vcf -stand_call_conf 20 -stand_emit_conf 20 -drf DuplicateRead -rf ReassignOneMappingQuality --reassign_mapping_quality_from 255 --reassign_mapping_quality_to 60 -o {output}"
#rule giremi:
#	input: PrintReads_BQSR_output = output_dir + "SRR{tag, +d}PrintReads_BQSR_output.bam", UnifiedGenotyper_output = output_dir + "SRR{tag, +d}UnifiedGenotyper_output.vcf"
#	output: output_dir + "SRR{tag, +d}giremi_output.res"
#	threads: 12
#	shell: "LD_LIBRARY_PATH=mnt/isilon/cbmi/variome/Taylor/giremi/htslib/; export LD_LIBRARY_PATH; {giremi} -f {ref_genome} -l {input.UnifiedGenotyper_output} -o {output} {input.PrintReads_BQSR_output}"