#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Kinnex preprocessing processes. meta per raw pool BAM:
 *   pool, movie, readset, bcM, inputBase, id (= <movie>.<bcM>), runfolder
 */


// skera: MAS array reads -> segmented reads (S-reads). One BAM in, one out.
process skera_split {
    tag   "${meta.id}"
    label 'high'
    conda "${params.condaEnvs.skera}"

    publishDir "${params.demuxDir}/QC/skera", mode: 'copy', pattern: '*.summary.csv'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.inputBase}.skera.bam"),         emit: bam
    tuple val(meta), path("${meta.inputBase}.skera.summary.csv"), emit: summary

    script:
    """
    skera split \\
        --num-threads ${task.cpus} \\
        ${bam} \\
        ${params.kinnexMas8Adapters} \\
        ${meta.inputBase}.skera.bam
    """

    stub:
    """
    touch ${meta.inputBase}.skera.bam ${meta.inputBase}.skera.summary.csv
    """
}


// lima: S-reads -> one BAM per IsoSeqX barcode pair found in this pool.
// Prefix carries movie AND bcM, so two pools on one SMRT cell never collide.
process lima {
    tag   "${meta.id}"
    label 'high'
    conda "${params.condaEnvs.lima}"

    publishDir "${params.demuxDir}/QC/lima", mode: 'copy', pattern: '*.lima.{summary,counts}'

    input:
    tuple val(meta), path(sreads)

    output:
    tuple val(meta), path("${meta.movie}.${meta.readset}.${meta.bcM}.lima.*--*.bam"),          emit: bams
    tuple val(meta), path("${meta.movie}.${meta.readset}.${meta.bcM}.lima.lima.summary"),
                     path("${meta.movie}.${meta.readset}.${meta.bcM}.lima.lima.counts"),       emit: summary

    script:
    def prefix = "${meta.movie}.${meta.readset}.${meta.bcM}.lima"
    """
    lima \\
        --isoseq \\
        --peek-guess \\
        --overwrite-biosample-names \\
        --num-threads ${task.cpus} \\
        ${sreads} \\
        ${params.kinnexBarcodes12plex} \\
        ${prefix}.bam
    """

    stub:
    // params.stubLimaBarcodes (test only, -params-file): bcM -> barcodes lima "found"
    def prefix = "${meta.movie}.${meta.readset}.${meta.bcM}.lima"
    def bcs    = ((params.stubLimaBarcodes ?: [:])[meta.bcM] ?: ['IsoSeqX_bc01_5p']).collect { "${prefix}.${it}--IsoSeqX_3p.bam" }.join(' ')
    """
    touch ${bcs} ${prefix}.lima.summary ${prefix}.lima.counts
    """
}


// Per-sample uBAM: copy of the lima BAM with every @RG SM set to the NPN.
// samtools is system-wide (no conda env).
process set_sample_name {
    tag   "${meta.npn}.${meta.bcM}.${meta.barcode5}"
    label 'low'

    publishDir "${params.demuxDir}", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.outName}"), emit: bam

    script:
    """
    samtools view -H ${bam} \\
        | awk -v sm='${meta.npn}' 'BEGIN{FS=OFS="\\t"} /^@RG/ {for(i=1;i<=NF;i++) if(\$i ~ /^SM:/) \$i="SM:" sm} {print}' \\
        | samtools reheader - ${bam} > ${meta.outName}
    samtools quickcheck ${meta.outName}
    n=\$(samtools view -H ${meta.outName} \\
        | awk -F'\\t' -v sm='SM:${meta.npn}' '/^@RG/ { ok = 0; for (i = 1; i <= NF; i++) if (\$i == sm) ok = 1; if (!ok) bad++ } END { print bad + 0 }')
    if [ "\$n" -ne 0 ]; then
        echo "ERROR: \$n @RG line(s) in ${meta.outName} without SM:${meta.npn}" >&2
        exit 1
    fi
    """

    stub:
    """
    touch ${meta.outName}
    """
}


// Run-level Kinnex QC report (skera + lima + barcode balance), one per raw pool BAM.
// barcodeSheet: "<npn>\t<IsoSeqX_bcNN_5p>" lines for this pool, built from the
// LabWare sheet in main.nf (kinnex_run_summary.py's samplesheet format).
process kinnex_run_summary {
    tag   "${meta.id}"
    label 'low'
    conda "${params.condaEnvs.somaticSummaryEnv}"

    publishDir "${params.demuxDir}/QC", mode: 'copy'
    // Shared QC report folder (lrsShared kinnexQcReportDir), layout as the old script:
    // <dir>/<runfolder name>/*.html and <dir>/<runfolder name>/yaml_json/*.{yaml,json}.
    // Never for --outdir (test) runs.
    publishDir { "${params.kinnexQcReportDir ?: params.demuxDir}/${meta.runfolder}" },
               mode: 'copy', overwrite: true, pattern: '*.kinnex_run_summary.html',
               enabled: (params.kinnexQcReportDir && !params.outdir) as boolean
    publishDir { "${params.kinnexQcReportDir ?: params.demuxDir}/${meta.runfolder}/yaml_json" },
               mode: 'copy', overwrite: true, pattern: '*.kinnex_run_summary.{yaml,json}',
               enabled: (params.kinnexQcReportDir && !params.outdir) as boolean

    input:
    tuple val(meta), path(skeraSummary), path(limaSummary), path(limaCounts), val(barcodeSheet)

    output:
    tuple val(meta), path("${meta.id}.barcodes.tsv"),               emit: sheet
    tuple val(meta), path("${meta.id}.kinnex_run_summary.yaml"),    emit: yaml
    tuple val(meta), path("${meta.id}.kinnex_run_summary.json"),    emit: json
    tuple val(meta), path("${meta.id}.kinnex_run_summary.html"),    emit: html

    script:
    """
    cat > ${meta.id}.barcodes.tsv <<'EOF'
${barcodeSheet}
EOF

    python3 ${params.kinnexSummaryPy} \\
        --skera-summary ${skeraSummary} \\
        --lima-summary  ${limaSummary} \\
        --lima-counts   ${limaCounts} \\
        --samplesheet   ${meta.id}.barcodes.tsv \\
        --movie         ${meta.movie} \\
        --array         ${meta.pool} \\
        --mas-barcode   ${meta.bcM} \\
        --runfolder     ${meta.runfolder} \\
        --html-template ${params.kinnexSummaryHtml} \\
        --output        ${meta.id}.kinnex_run_summary
    """

    stub:
    """
    cat > ${meta.id}.barcodes.tsv <<'EOF'
${barcodeSheet}
EOF
    touch ${meta.id}.kinnex_run_summary.yaml ${meta.id}.kinnex_run_summary.json ${meta.id}.kinnex_run_summary.html
    """
}


// Fails the run when a sheet sample got no demultiplexed BAM (unless
// --allowMissingSamples). The manifest itself is written before this runs,
// so it is on disk whatever the outcome.
process demux_check {
    tag   "${params.ssBase}"
    label 'low'
    executor 'local'

    input:
    path(manifest)

    output:
    path("demux_check.ok")

    script:
    def allow = params.allowMissingSamples ? 1 : 0
    """
    missing=\$(awk -F'\\t' 'NR > 1 && \$5 == "NO_READS"' ${manifest} | wc -l)
    unassigned=\$(awk -F'\\t' 'NR > 1 && \$5 == "UNASSIGNED"' ${manifest} | wc -l)
    echo "missing=\$missing unassigned=\$unassigned" > demux_check.ok
    if [ "\$unassigned" -gt 0 ]; then
        echo "WARNING: \$unassigned lima barcode(s) have no sheet sample (see ${manifest})" >&2
    fi
    if [ "\$missing" -gt 0 ] && [ ${allow} -eq 0 ]; then
        echo "ERROR: \$missing sheet sample(s) have no demultiplexed BAM (NO_READS in ${manifest})." >&2
        echo "Rerun with --allowMissingSamples to accept this." >&2
        exit 1
    fi
    """
}
