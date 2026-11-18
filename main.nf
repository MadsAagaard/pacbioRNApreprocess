#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * =============================================================================
 * pacbioRNApreprocess — KG Vejle PacBio Kinnex RNA preprocessing
 * =============================================================================
 * One Revio run in, demultiplexed per-sample RNA uBAMs out:
 *
 *   <runfolder>/hifi_reads/<pool>.<movie>.hifi_reads.bcM000N.bam      (raw MAS arrays)
 *     skera split  -> S-reads
 *     lima         -> one BAM per IsoSeqX barcode
 *     SM = NPN     -> <demuxDir>/<npn>.<movie>.hifi_reads.<bcM>.lima.<IsoSeqX_bcNN_5p>.SMfixed.bam
 *
 * demuxDir defaults to <runfolder>/rna_demultiplex, where pacbioSomatic (RNA arm)
 * and pacbioGermline rna.nf find the uBAMs (glob **\/*.lima.*.bam, NPN = first
 * '.' token).
 *
 * Sample sheet: the run's LabWare metadata sheet (same file as germline DNA).
 * See SAMPLESHEET below for which rows are RNA and how they map to barcodes.
 *
 * Run:
 *   nextflow run MadsAagaard/pacbioRNApreprocess -r main -profile slurm \
 *       --samplesheet <yymmdd>_<runName>_<NN>_metadata.txt \
 *       --runfolder   <raw data archive>/<year>/<runLabel>_<RUN_NAME>
 *
 * Every server path (references, conda envs, report folders) comes from the
 * local lrsShared config (nextflow.config). This repository holds none.
 * =============================================================================
 */

include { skera_split; lima; set_sample_name; kinnex_run_summary; demux_check } from './modules/preprocessModules.nf'


// =============================================================================
// SHARED CONFIG KEYS  (server paths live in lrsShared only; see nextflow.config)
// =============================================================================
// kinnexQcReportDir may be set to false/null in lrsShared to switch the copy off,
// but the key must exist, so a shared version without it is caught here.
def requiredShared = ['kinnexMas8Adapters', 'kinnexBarcodes12plex', 'kinnexQcReportDir']
def missingShared  = requiredShared.findAll { !params.containsKey(it) } +
                     ['skera', 'lima'].findAll { !(params.condaEnvs instanceof Map && params.condaEnvs[it]) }
                                      .collect { "condaEnvs.${it}" }
if (missingShared) {
    error("pacbioRNApreprocess: ${params.lrsSharedDir}/lrsShared.config lacks: ${missingShared.join(', ')}. " +
          "Add the Kinnex block to that lrsShared version (see README, Config).")
}
if (!workflow.stubRun) {
    file(params.kinnexMas8Adapters,   checkIfExists: true)
    file(params.kinnexBarcodes12plex, checkIfExists: true)
}

if (params.help) {
    log.info """\
    pacbioRNApreprocess ${workflow.manifest.version}
      --samplesheet <file>      LabWare metadata sheet of the run (required)
      --runfolder <dir>         archived run folder with hifi_reads/*.hifi_reads.bcM*.bam (required)
      --outdir <dir>            output folder (default <runfolder>/rna_demultiplex; no completion flag in the run folder)
      --allowMissingSamples     do not fail when a sheet sample gets no demultiplexed BAM
      -profile slurm
    """.stripIndent()
    exit 0
}


// =============================================================================
// SAMPLESHEET  (LabWare metadata, parsed at launch)
// =============================================================================
// LabWare writes an RNA (Kinnex) sample with its IsoSeqX barcode in the proband
// slot, e.g. (seen 2026-11-17):
//   DNA  0000000001_100000000001_78_SL-TESTLIST_M_F_100000000091_3
//   RNA  0000000002_100000000002_71_SL-TESTLIST-RNA_M_IsoSeqX-bc03-5p_100000000092
// A row is RNA when the line contains an IsoSeqX barcode token (KINNEX_RE:
// case-insensitive, '-', '_' or ' ' separators, matched anywhere, so a later
// move to its own column still works). Every other row (DNA) is ignored.
//
// NPN = field 2 of the LabWare line (rekv_npn_material_...), underscore-joined in
// column 1 or tab-exploded over the columns - the same rule as pacbioGermline.
// intRef = the field after the IsoSeqX token, kept in the manifest only.
// It is not the person's DNA NPN (LabWare never writes that on an RNA row).
//
// Pool: a bcM0001..bcM9999 token on the row ties the sample to that raw BAM.
// LabWare does not write one today, so the run must then have exactly one bcM
// pool (one Kinnex array per run so far). Mixed (some rows with, some
// without) is refused.

def KINNEX_RE = ~/(?i)IsoSeqX[_ -]?bc0*(\d{1,2})(?![0-9])/
def MAS_RE    = ~/(?<![A-Za-z0-9])(bcM\d{4})(?![0-9])/
def MOVIE_RE  = ~/^m\d+[a-z]?_\d{6}_\d{6}(_s\d+)?$/
def MAX_BC    = 12                        // 12-plex IsoSeqX kit

def problems = []
def warnings = []

def sheetFile = file(params.samplesheet, checkIfExists: true)
def rnaRows   = []

sheetFile.readLines().eachWithIndex { raw, i ->
    def line = raw.replaceAll(/\r$/, '')
    if (!line.trim() || line.startsWith('#')) return

    def bcNums = (line =~ KINNEX_RE).collect { it[1] as int }.unique()
    if (!bcNums) return                                   // DNA row
    def where  = "sheet line ${i + 1}"

    if (bcNums.size() > 1) { problems << "${where}: more than one IsoSeqX barcode (${bcNums}): ${line}"; return }
    if (bcNums[0] < 1 || bcNums[0] > MAX_BC) { problems << "${where}: IsoSeqX barcode ${bcNums[0]} outside 1-${MAX_BC}: ${line}"; return }

    def masList = (line =~ MAS_RE).collect { it[1] }.unique()
    if (masList.size() > 1) { problems << "${where}: more than one bcM pool adapter (${masList}): ${line}"; return }

    def cols   = line.split('\t') as List
    def fields = cols[0].contains('_') ? cols[0].tokenize('_') : cols
    def npn    = fields.size() > 1 ? fields[1].toString().trim() : ''
    if (!(npn ==~ /[A-Za-z0-9-]+/)) { problems << "${where}: cannot read an NPN from field 2 ('${npn}'): ${line}"; return }

    def iBc    = fields.findIndexOf { it.toString() =~ /(?i)^\s*IsoSeqX/ }
    def intRef = (iBc >= 0 && iBc + 1 < fields.size()) ? fields[iBc + 1].toString().trim() : 'NA'

    rnaRows << [
        npn       : npn,
        rekv      : fields[0].toString().trim(),
        intRef    : intRef ?: 'NA',
        barcode5  : String.format('IsoSeqX_bc%02d_5p', bcNums[0]),
        bcM       : masList ? masList[0] : null,
        sheetLine : i + 1,
    ]
}

if (!rnaRows) {
    problems << "no RNA rows in ${sheetFile.name}: no line contains an IsoSeqX barcode"
}


// =============================================================================
// RAW KINNEX BAMs  (one listing at launch)
// =============================================================================
// Non-recursive on purpose: <runfolder>/rna_demultiplex/ also holds files with
// '.hifi_reads.bcM' in their names (this pipeline's output).

def runDir    = file(params.runfolder, checkIfExists: true)
def rawFiles  = files("${runDir}/hifi_reads/*.hifi_reads.bcM*.bam")
def poolMetas = []

rawFiles.sort { it.name }.each { bam ->
    def tokens = bam.name.replaceFirst(/\.bam$/, '').tokenize('.')
    def iMovie = tokens.findIndexOf { it ==~ MOVIE_RE }
    def bcM    = tokens.find { it ==~ /bcM\d{4}/ }
    if (iMovie < 1 || !bcM) { problems << "raw BAM name not <pool>.<movie>.hifi_reads.bcM000N.bam: ${bam}"; return }
    poolMetas << [
        pool      : tokens[0..<iMovie].join('.'),
        movie     : tokens[iMovie],
        readset   : 'hifi_reads',
        bcM       : bcM,
        inputBase : bam.name.replaceFirst(/\.bam$/, ''),
        id        : "${tokens[iMovie]}.${bcM}",
        runfolder : runDir.name,
        path      : bam.toString(),
    ]
}

if (!rawFiles) {
    problems << "no ${runDir}/hifi_reads/*.hifi_reads.bcM*.bam (no Kinnex pool in this run)"
}

def runBcMs = poolMetas*.bcM.unique().sort()


// =============================================================================
// SHEET <-> RUN CONSISTENCY  (everything here stops the run before any task)
// =============================================================================
def withBcM = rnaRows.count { it.bcM }

if (rnaRows && withBcM == 0) {
    if (runBcMs.size() == 1) {
        rnaRows.each { it.bcM = runBcMs[0] }
    } else if (runBcMs.size() > 1) {
        problems << "sheet rows name no bcM pool, but the run has ${runBcMs.size()} pools (${runBcMs.join(', ')}): " +
                    "barcodes cannot be assigned. Add the bcM adapter to the RNA rows."
    }
} else if (rnaRows && withBcM < rnaRows.size()) {
    problems << "${rnaRows.size() - withBcM} of ${rnaRows.size()} RNA rows name no bcM pool; either all or none must"
}

rnaRows.groupBy { it.npn }.findAll { k, v -> v.size() > 1 }.each { npn, v ->
    problems << "NPN ${npn} on more than one RNA row (sheet lines ${v*.sheetLine.join(', ')})"
}
rnaRows.findAll { it.bcM }.groupBy { "${it.bcM}|${it.barcode5}" }.findAll { k, v -> v.size() > 1 }.each { key, v ->
    problems << "pool/barcode ${key.replace('|', ' ')} assigned to more than one sample (${v*.npn.join(', ')})"
}
rnaRows.findAll { it.bcM && !(it.bcM in runBcMs) }.each { r ->
    problems << "NPN ${r.npn} (sheet line ${r.sheetLine}): pool ${r.bcM} has no raw BAM in ${runDir}/hifi_reads"
}

def sheetBcMs = rnaRows*.bcM.findAll().unique()
(runBcMs - sheetBcMs).each { b ->
    warnings << "pool ${b} has raw BAM(s) but no RNA row in the sheet: not processed"
}

if (problems) {
    error("pacbioRNApreprocess: cannot start\n  - " + problems.join("\n  - "))
}
warnings.each { log.warn it }

def usedPools = poolMetas.findAll { it.bcM in sheetBcMs }


log.info """\
======================================================
pacbioRNApreprocess ${workflow.manifest.version}  (Kinnex skera + lima)
======================================================
Samplesheet    : ${params.samplesheet}  (${rnaRows.size()} RNA rows)
Runfolder      : ${runDir}
Pools          : ${usedPools.collect { "${it.bcM} (${it.movie})" }.join(', ')}
Output         : ${params.demuxDir}
QC report copy : ${params.outdir ? 'none (--outdir)' : (params.kinnexQcReportDir ? "${params.kinnexQcReportDir}/${runDir.name}" : 'none')}
workDir        : ${workflow.workDir}
lrsShared      : ${params.lrsSharedDir}
lrsCustomTools : ${params.lrsCustomToolsBase}
error mode     : ${params.errorMode}
"""


// =============================================================================
// WORKFLOW
// =============================================================================
workflow {

    Channel.fromList(usedPools)
        | map { m -> tuple(m.findAll { k, v -> k != 'path' }, file(m.path)) }
        | skera_split

    lima(skera_split.out.bam)

    // one item per lima BAM, keyed by pool + 5' barcode
    lima.out.bams
        | flatMap { meta, bams ->
            (bams instanceof List ? bams : [bams]).collect { bam ->
                def m = bam.name =~ /\.lima\.([^.]+?)--([^.]+)\.bam$/
                if (!m.find()) error("unexpected lima output name: ${bam.name}")
                def bc5 = m.group(1)
                tuple("${meta.bcM}|${bc5}", meta + [barcode5: bc5, barcode3: m.group(2)], bam)
            }
        }
        | set { limaSplit }

    // combine (not join): a pool sequenced on two SMRT cells gives one lima BAM
    // per movie for the same key, and each must reach the sample.
    Channel.fromList(rnaRows)
        | map { r -> tuple("${r.bcM}|${r.barcode5}", r) }
        | combine(limaSplit, by: 0)
        | map { key, r, meta, bam ->
            def m = meta + [npn: r.npn, rekv: r.rekv, sheetLine: r.sheetLine]
            tuple(m + [outName: "${m.npn}.${m.movie}.${m.readset}.${m.bcM}.lima.${m.barcode5}.SMfixed.bam"], bam)
        }
        | set_sample_name

    // ---- run-level QC: one report per raw pool BAM ---------------------------
    def sheetByBcM = rnaRows.groupBy { it.bcM }.collectEntries { b, rs ->
        [b, rs.sort { it.barcode5 }.collect { "${it.npn}\t${it.barcode5}" }.join('\n')]
    }

    skera_split.out.summary
        | map { meta, s -> tuple(meta.id, meta, s) }
        | join(lima.out.summary.map { meta, s, c -> tuple(meta.id, s, c) }, failOnMismatch: true)
        | map { id, meta, s, ls, lc -> tuple(meta, s, ls, lc, sheetByBcM[meta.bcM]) }
        | kinnex_run_summary

    // ---- manifest: every sheet sample, every lima barcode ---------------------
    // status OK | NO_READS (sheet sample, no lima BAM) | UNASSIGNED (lima BAM, no sheet sample)
    def demuxDir = params.demuxDir

    set_sample_name.out.bam
        | map { m, bam -> [key: "${m.bcM}|${m.barcode5}", npn: m.npn, movie: m.movie, file: "${demuxDir}/${bam.name}"] }
        | collect | ifEmpty([])
        | map { l -> tuple('all', l) }
        | set { namedList }

    limaSplit
        | map { key, meta, bam -> [key: key, movie: meta.movie] }
        | collect | ifEmpty([])
        | map { l -> tuple('all', l) }
        | set { limaList }

    namedList
        | join(limaList)
        | map { k, named, limas ->
            def rows   = ["npn\tbcM\tbarcode5\tmovie\tstatus\tbam\tintRef"]
            def sheetKeys = rnaRows.collect { "${it.bcM}|${it.barcode5}" } as Set
            rnaRows.sort { [it.bcM, it.barcode5] }.each { r ->
                def key  = "${r.bcM}|${r.barcode5}"
                def hits = named.findAll { it.key == key }.sort { it.movie }
                if (hits) hits.each { h -> rows << "${r.npn}\t${r.bcM}\t${r.barcode5}\t${h.movie}\tOK\t${h.file}\t${r.intRef}" }
                else      rows << "${r.npn}\t${r.bcM}\t${r.barcode5}\tNA\tNO_READS\tNA\t${r.intRef}"
            }
            limas.findAll { !(it.key in sheetKeys) }.sort { [it.key, it.movie] }.each { l ->
                def (b, bc) = l.key.tokenize('|')
                rows << "NA\t${b}\t${bc}\t${l.movie}\tUNASSIGNED\tNA\tNA"
            }
            rows.join('\n') + '\n'
        }
        | collectFile(name: "${params.ssBase}.demux_manifest.tsv", storeDir: "${params.demuxDir}/QC")
        | demux_check
}


// =============================================================================
// COMPLETION
// =============================================================================
// The flag tells downstream automation that <runfolder>/rna_demultiplex is final.
// Only for a successful run writing into the run folder itself (not --outdir).
workflow.onComplete {

    if (!workflow.success) {
        log.error "pacbioRNApreprocess failed - no completion flag written."
        return
    }
    if (params.outdir) {
        log.info "--outdir given: no completion flag written into the run folder."
        return
    }

    def flag = file("${params.runfolder}/kgvejle_rna_demultiplex_complete.txt")
    flag.text = """\
        pipeline    : ${workflow.manifest.name} ${workflow.manifest.version} (${workflow.revision ?: 'local'} ${workflow.commitId ?: ''})
        runfolder   : ${params.runfolder}
        samplesheet : ${params.samplesheet}
        output      : ${params.demuxDir}
        manifest    : ${params.demuxDir}/QC/${params.ssBase}.demux_manifest.tsv
        completed   : ${workflow.complete}
        duration    : ${workflow.duration}
        """.stripIndent()
    log.info "Completion flag written: ${flag}"
}
