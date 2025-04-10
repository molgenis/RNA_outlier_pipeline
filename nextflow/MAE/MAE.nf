#!/usr/bin/env nextflow
/**
Nextflow MAE workflow
Uses Gagneur-lab MAE derived implementation from the DROP pipeline.. 
Results will be stored in a resultsfile (*.tsv)
**/

nextflow.enable.dsl=2


process MAEreadCounting {
    time '12h'
    memory '8 GB'
    cpus 1

    publishDir "$params.output/MAE/ASEreadcounts", mode: 'copy'

    input:
        tuple val(sampleID), path(vcf), path(bamFile)
        tuple val(fasta), path(fastafolder)

    output:
        val "${sampleID}"
        path "${sampleID}_maecounts.tsv"

    script:
    """
    set -e
    set -u

    # From DROP:
    if samtools view -H "${bamFile}" | grep -q "@RG";then
    printf "BAM contains read groups (RG), continuing with ASEReadCounter...\n"
    else
    printf "%s\n" "" "ERROR: BAM file doesn't contain Read Group Tag (RG)" \
    " RG doesn't exist, it can be added using -" \
    "   gatk AddOrReplaceGroups -R /path/to/reference -I /your/input.bam -O /your/output.bam --QUIET true" \
    " https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-" \
    " Try rerunning this module using the BAM with RG tags"
    exit 1
    fi

    # remove INFO field, it sometimes contains whitespace or invalid line lengths. 
    bcftools annotate -x 'INFO' "$vcf" -o ${sampleID}_noinfo.vcf

    # create index
    gatk IndexFeatureFile -I "${sampleID}_noinfo.vcf"

    # Only select hetzyg variants and only SNP's
    gatk SelectVariants -V ${sampleID}_noinfo.vcf -O "${sampleID}_temp.vcf" -R "${fasta}" --restrict-alleles-to BIALLELIC -select 'vc.getHetCount()==1' --select-type-to-include SNP

    # remove duplicates
    bcftools norm --rm-dup all "${sampleID}_temp.vcf" > "${sampleID}_temp.vcf.tmp" && mv "${sampleID}_temp.vcf.tmp" "${sampleID}_temp.vcf"

    # create new index
    gatk IndexFeatureFile -I "${sampleID}_temp.vcf"

    gatk ASEReadCounter -I $bamFile -V "${sampleID}_temp.vcf" -O "${sampleID}_temp_counts.tsv" -R "${fasta}"

    # sorting the count file.
    (head -n 1 "${sampleID}_temp_counts.tsv" && tail -n +2 "${sampleID}_temp_counts.tsv" | sort -k1,1 -V -s ) > "${sampleID}_maecounts.tsv"
    
    """
}

process GetMAEresults {

    time '1h'
    memory '4 GB'
    cpus 1

    publishDir "$params.output/MAE/DeSEQresults", mode: 'copy'

    input:
        val sampleid
        path asecounts
        path resultsR

    output:
        
        path "${sampleid}_result_mae.tsv"

    script:
        """
        Rscript "${resultsR}" "${asecounts}" "${sampleid}_result_mae.tsv" "${sampleid}"
        """

}


