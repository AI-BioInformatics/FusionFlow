#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run federicacitarrella/pipelineGeneFusions --rnareads '*_{1,2}.fastq.gz' --dnareads_tumor '*_{3,4}.fastq.gz' --dnareads_normal '*_{5,6}.fastq.gz' -profile docker

    Principal arguments:
      --rnareads [file]               Path to input RNA data (must be surrounded with quotes)
                                      Format: fastq
      --dnareads_tumor [file]         Path to input tumor DNA data (must be surrounded with quotes)
                                      Format: fastq, bam
      --dnareads_normal [file]        Path to input normal DNA data (must be surrounded with quotes)
                                      Format: fastq, bam
      -profile [str]                  Configuration profile to use.
                                      Available: docker, local, test_docker, test_local

    Tool flags:
      --arriba [bool]                 Run Arriba (RNA)
      --ericscript [bool]             Run Ericscript (RNA)
      --fusioncatcher [bool]          Run FusionCatcher (RNA)
      --integrate [bool]              Run Integrate (RNA; RNA+DNA)
      --genefuse [bool]               Run GeneFuse (DNA)

    References:
      --ericscript_ref [file]         Path to EricScript reference
      --arriba_ref [dir]              Path to Arriba reference
      --fusioncatcher_ref [dir]       Path to FusionCatcher reference
      --integrate_ref [file]          Path to Integrate reference
      --genefuse_ref [file]           Path to GeneFuse reference

    Options:
      --dnabam [bool]                 Specifies that the wgs input has bam format
      --nthreads [int]                Specifies the number of threads [8]

    Other options:
      --outdir [dir]                  The output directory where the results will be saved

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/* DATA CHANNELS CREATION */

// rna
if (params.rnareads){
    Channel.fromFilePairs(params.rnareads, checkIfExists: true)
    .into{rna_reads_ericscript; rna_reads_arriba; rna_reads_fusioncatcher; rna_reads_integrate}
}

// dna-tumor
if (params.dnareads_tumor){
    Channel.fromFilePairs(params.dnareads_tumor, checkIfExists: true)
    .into{dna_reads_tumor_integrate; dna_reads_tumor_genefuse}
}

// dna-normal
if (params.dnareads_normal){
    Channel.fromFilePairs(params.dnareads_normal, checkIfExists: true)
    .into{dna_reads_normal_integrate}
}

// Variabili integrate che prima o poi capirò
if (params.dnareads_tumor) {
    integrateWGSt = true
    command1 = "dna.tumor.bam"
}
else {
    integrateWGSt = false
    command1 = ""
}

if (params.dnareads_normal) {
    integrateWGSn = true
    command2 = "dna.normal.bam"
}
else {
    integrateWGSn = false
    command2 = ""
}

/* REFERENCE GENOME */

// download
process referenceGenome_downloader{
    tag "Downloading"
    storeDir "${params.outdir}/reference_genome"
    label "low_ram"

    output:
    file "hg38.fa" into ref_gen_ch, ref_gen_integrate_ch, ref_gen_integrate_ch2, ref_gen_integrate_ch3, ref_gen_genefuse_ch

    when: params.integrate || params.genefuse || !(params.arriba || params.ericscript || params.fusioncatcher || params.genefuse || params.integrate)

    script:
    """
    #!/bin/bash

    curl -O -J -L https://osf.io/yevub/download
    """
}

// index
process referenceGenome_index{
    tag "Indexing"
    storeDir "${params.outdir}/reference_genome"
    label "low_ram"

    input:
    file 'hg38.fa' from ref_gen_ch

    output:
    file "hg38.fa.amb" into index_amb_ch
    file "hg38.fa.ann" into index_ann_ch
    file "hg38.fa.bwt" into index_bwt_ch
    file "hg38.fa.pac" into index_pac_ch
    file "hg38.fa.sa" into index_sa_ch

    when: params.integrate || !(params.arriba || params.ericscript || params.fusioncatcher || params.genefuse || params.integrate)

    shell:
    '''
    #!/bin/bash

    export PATH="!{params.envPath_integrate}:$PATH"
    bwa index -a bwtsw hg38.fa
    '''
}

/* PIPELINE */

/* ERICSCRIPT */

process ericscript_downloader{
    tag "Downloading"
    storeDir "${params.outdir}/ericscript/files"
    label "low_ram"

    output:
    file "ericscript_db_homosapiens_ensembl84" into ericscript_download

    when: params.ericscript || !(params.arriba || params.ericscript || params.fusioncatcher || params.genefuse || params.integrate)

    script:
    """
    #!/bin/bash

    curl -O -J -L https://osf.io/54s6h/download
    tar -xf ericscript_db_homosapiens_ensembl84.tar.bz2
    rm ericscript_db_homosapiens_ensembl84.tar.bz2
    """
}

process ericscript{
    tag "${pair_id}"
    publishDir "${params.outdir}/ericscript/output"
    label "low_ram"

    input:
    tuple pair_id, file(rna_reads) from rna_reads_ericscript
    file(ericscript_db) from ericscript_download

    output:
    file "${pair_id}"

    when: params.ericscript || !(params.arriba || params.ericscript || params.fusioncatcher || params.genefuse || params.integrate)

    script:
    """
    #!/bin/bash

    export PATH="${params.envPath_ericscript}:$PATH"
    ericscript.pl -o ${pair_id} -db ${params.ericscript_ref}/${ericscript_db} ${rna_reads}
    """
}


/* ARRIBA */

process arriba_downloader{
    tag "Downloading"
    publishDir "${params.outdir}/arriba/files"
    label "high_ram"

    output:
    file "ENSEMBL93.gtf" into ens_arriba
    file "GRCh38.fa" into gr38_arriba_ch
    file "STAR_index_GRCh38_ENSEMBL93" into star_index_ch

    when: params.arriba || !(params.arriba || params.ericscript || params.fusioncatcher || params.genefuse || params.integrate)

    script:
    """

    #!/bin/bash

    export PATH="${params.envPath_arriba}bin:$PATH"
    ${params.envPath_arriba}var/lib/arriba/download_references.sh GRCh38+ENSEMBL93

    """
}

process arriba{
    tag "${pair_id}"
    publishDir "${params.outdir}/arriba/output"
    label "low_ram"

    input:
    tuple pair_id, file(rna_reads) from rna_reads_arriba
    file(arriba_ref) from arriba_download

    output:
    file "${pair_id}" into arriba_fusions

    when: params.arriba || !(params.arriba || params.ericscript || params.fusioncatcher || params.genefuse || params.integrate)

    script:
    """
    #!/bin/bash

    export PATH="${params.envPath_arriba}/bin:$PATH"

    run_arriba.sh ${arriba_ref}/STAR_index_GRCh38_ENSEMBL93/ ${arriba_ref}/ENSEMBL93.gtf ${arriba_ref}/GRCh38.fa ${params.envPath_arriba}var/lib/arriba/blacklist_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz ${params.envPath_arriba}var/lib/arriba/known_fusions_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz ${params.envPath_arriba}var/lib/arriba/protein_domains_hg19_hs37d5_GRCh37_v2.1.0.gff3 ${params.nthreads} ${rna_reads}

    mkdir ${pair_id}
    mv *.out ${pair_id}
    mv *.tsv ${pair_id}
    mv *.out ${pair_id}
    mv *bam* ${pair_id}

    """
}

/* FUSIONCATCHER */

process fusioncatcher_downloader {
    tag "Downloading"
    storeDir "${params.outdir}/fusioncatcher/files"
    label "low_ram"

    output:
    file "human_v102.tar.gz.aa" into fusioncatcher_download_aa
    file "human_v102.tar.gz.ab" into fusioncatcher_download_ab
    file "human_v102.tar.gz.ac" into fusioncatcher_download_ac
    file "human_v102.tar.gz.ad" into fusioncatcher_download_ad

    when: params.fusioncatcher || !(params.arriba || params.ericscript || params.fusioncatcher || params.genefuse || params.integrate)

    shell:
    '''
    #!/bin/bash

    wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.aa
    wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ab
    wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ac
    wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ad

    '''
}

process fusioncatcher{
    tag "${pair_id}"
    publishDir "${params.outdir}/fusioncatcher/output"
    label "high_ram" // >24GB

    input:
    tuple pair_id, file(rna_reads)  from rna_reads_fusioncatcher
    file(fusioncatcher_db) from fusioncatcher_download_aa.mix(fusioncatcher_download_ab, fusioncatcher_download_ac, fusioncatcher_download_ad)

    output:
    file "${pair_id}" optional true into fusioncatcher_fusions

    when: params.fusioncatcher || !(params.arriba || params.ericscript || params.fusioncatcher || params.genefuse || params.integrate)

    script:

    """
    #!/bin/bash

    export PATH="${params.envPath_fusioncatcher}:$PATH"
    fusioncatcher -d ${params.fusioncatcher_ref}/${fusioncatcher_db}/human_v102 -i ${rna_reads} -o ${pair_id}
    """
}

/* INTEGRATE */

process integrate_downloader{
    tag "Downloading"
    storeDir "${params.outdir}/integrate/files"
    label "low_ram"

    output:
    file "GRCh38_noalt_as" into grch38_ch, grch38_ch2, grch38_ch3
    file "annot.refseq.txt" into annot_ch, annot_ch2, annot_ch3
    file "INTEGRATE_0_2_6" into int_ch, int_ch2, int_ch3

    when: params.integrate || !(params.arriba || params.ericscript || params.fusioncatcher || params.genefuse || params.integrate)

    shell:
    '''
    #!/bin/bash

    wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
    unzip GRCh38_noalt_as.zip
    rm GRCh38_noalt_as.zip

    curl -O -J -L https://osf.io/dgvcx/download

    curl -O -J -L https://osf.io/gv7sq/download
    tar -xvf INTEGRATE.0.2.6.tar.gz
    rm INTEGRATE.0.2.6.tar.gz
    '''
}

// process integrate_converter{
//     tag "${pair_id}"
//     publishDir "${params.outdir}/integrate"
//
//     input:
//     tuple pair_id, file(rna_reads) from rna_reads_integrate
//     file(integrate_db) from grch38_ch2.mix(annot_ch2, int_ch2)
//     file(refgen) from ref_gen_integrate_ch
//     file(index) from index_amb_ch.mix(index_ann_ch, index_bwt_ch, index_pac_ch, index_sa_ch)
//     file(wgstinput) from dna_reads_tumor_integrate
//     file(wgsninput) from dna_reads_normal_integrate
//
//     output:
//     tuple pair_id, file("input/${pair_id}") into integrate_input
//
//     when: params.integrate || !(params.arriba || params.ericscript || params.fusioncatcher || params.genefuse || params.integrate)
//
//     script:
//     """
//     #!/bin/bash
//
//     export PATH="${params.envPath_integrate}:$PATH"
//
//     tophat --no-coverage-search ${integrate_db}/GRCh38_noalt_as/GRCh38_noalt_as ${rna_reads}
//
//     mkdir input && mkdir input/${pair_id}
//
//     cp tophat_out/accepted_hits.bam input/${pair_id}
//     cp tophat_out/unmapped.bam input/${pair_id}
//
//     if ${params.dnabam}; then
//       if ${integrateWGSt}; then
//         cp ${wgstinput} input/${pair_id}/dna.tumor.bam
//       fi
//       if ${integrateWGSn}; then
//         cp ${wgsninput} input/${pair_id}/dna.normal.bam
//       fi
//     elif ${integrateWGSt} || ${integrateWGSn}; then
//       mkdir index_dir
//       cp ${index}/* index_dir
//       if ${integrateWGSt}; then
//         bwa mem index_dir/hg38.fa ${wgstinput} | samtools sort -o input/${pair_id}/dna.tumor.bam
//       fi
//       if ${integrateWGSn}; then
//         bwa mem index_dir/hg38.fa ${wgsninput} | samtools sort -o input/${pair_id}/dna.normal.bam
//       fi
//     fi
//
//     """
// }
//
// process integrate{
//     tag "${pair_id}"
//     publishDir "${params.outdir}/integrate", mode: 'copy'
//
//     input:
//     tuple pair_id, file(input) from integrate_input
//     file(integrate_db) grch38_ch3.mix(annot_ch3, int_ch3)
//     file(refgen) from ref_gen_integrate_ch3
//     file(bwts) from bwts_integrate
//
//     output:
//     file "output/${pair_id}" optional true into integrate_fusions
//
//     when: params.integrate || !(params.arriba || params.ericscript || params.fusioncatcher || params.genefuse || params.integrate)
//
//     shell:
//     '''
//     #!/bin/bash
//
//     export PATH="!{params.envPath_integrate}:$PATH"
//
//     cp !{input}/* .
//
//     parallel samtools index ::: *.bam
//
//     LD_LIBRARY_PATH=/usr/local/lib
//     LD_LIBRARY_PATH=$LD_LIBRARY_PATH:!{integrate_db}/INTEGRATE_0_2_6/INTEGRATE-build/vendor/src/libdivsufsort-2.0.1-build/lib/
//     export LD_LIBRARY_PATH
//
//     !{integrate_db}/INTEGRATE_0_2_6/INTEGRATE-build/bin/Integrate fusion !{refgen} !{integrate_db}/annot.refseq.txt !{bwts} accepted_hits.bam unmapped.bam !{command1} !{command2}
//
//     mkdir output && mkdir output/!{pair_id}
//     cp *.tsv output/!{pair_id}
//     cp *.txt output/!{pair_id}
//
//     '''
// }

/* GENEFUSE */

process genefuse_downloader{
    tag "Downloading"
    storeDir "${params.outdir}/genefuse/files"
    label "low_ram"

    output:
    file "druggable.hg38.csv" into druggable_ch
    file "genefuse" into genefuse_download

    when: params.genefuse || !(params.arriba || params.ericscript || params.fusioncatcher || params.genefuse || params.integrate)

    shell:
    '''
    #!/bin/bash

    curl -O -J -L https://osf.io/8r9fh/download
    chmod a+x ./genefuse

    curl -O -J -L https://osf.io/jqywz/download
    '''
}


// check if low ram when converter is on
process genefuse_converter{
    tag "${pair_id}"
    publishDir "${params.outdir}/genefuse/input"
    label "low_ram"

    input:
    tuple pair_id, file(wgstinput) from dna_reads_tumor_genefuse

    output:
    tuple pair_id, file("${pair_id}") into genefuse_input

    when: params.genefuse || !(params.arriba || params.ericscript || params.fusioncatcher || params.genefuse || params.integrate)

    script:
    """
    #!/bin/bash

    export PATH="${params.envPath_genefuse}:$PATH"

    mkdir ${pair_id}

    if ${params.dnabam} && ${integrateWGSt}; then
        samtools sort -n  -o ${wgstinput}
        samtools fastq -@ ${params.nthreads} ${wgstinput} -1 ${pair_id}_3.fq.gz -2 ${pair_id}_4.fq.gz -0 /dev/null -s /dev/null -n
        cp *.fq.gz ${pair_id}
    elif ${integrateWGSt}; then
        cp ${wgstinput} ${pair_id}
    fi
    """

}

process genefuse{
    tag "${pair_id}"
    publishDir "${params.outdir}/genefuse/output"
    label "high_ram"

    input:
    tuple pair_id, file(input) from genefuse_input
    file(refgen) from ref_gen_genefuse_ch
    file(genefuse) from genefuse_download
    file(druggable) from druggable_ch


    when: params.genefuse || !(params.arriba || params.ericscript || params.fusioncatcher || params.genefuse || params.integrate)

    script:
    """
    #!/bin/bash

    export PATH="${params.envPath_genefuse}:$PATH"
    ${params.genefuse_ref}/${genefuse} -r ${params.referenceGenome}/${refgen} -f ${params.genefuse_ref}/${druggable} -1 ${pair_id}_3.fq.gz -2 ${pair_id}_4.fq.gz -h report.html > result
    """
}