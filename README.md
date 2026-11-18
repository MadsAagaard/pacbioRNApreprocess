# pacbioRNApreprocess

KG Vejle, Klinisk Genetik. Preprocessing of PacBio Kinnex (MAS-seq) RNA from a Revio run into
per-sample uBAMs for the RNA arms of pacbioSomatic and pacbioGermline (`rna.nf`).

```
<runfolder>/hifi_reads/<array>.<movie>.hifi_reads.bcM000N.bam      raw MAS arrays, one per Kinnex pool
  skera split   (mas8 primers)          -> S-reads
  lima --isoseq (IsoSeqX 12-plex)       -> one BAM per barcode
  @RG SM = NPN                          -> <runfolder>/rna_demultiplex/<npn>.<movie>.hifi_reads.<bcM>.lima.<IsoSeqX_bcNN_5p>.SMfixed.bam
  kinnex_run_summary.py (lrsCustomTools)-> rna_demultiplex/QC/<movie>.<bcM>.kinnex_run_summary.{html,json,yaml}
```

Downstream pipelines find the uBAMs with `**/*.lima.*.bam` and take the NPN from the first `.` token.

## Run

```
nextflow run MadsAagaard/pacbioRNApreprocess -r main -profile slurm \
    --samplesheet <yymmdd>_<runName>_<NN>_metadata.txt \
    --runfolder   <raw data archive>/<year>/<runLabel>_<RUN_NAME>
```

Options: `--outdir <dir>` (testing: output elsewhere; no completion flag in the run folder and no
copy to the shared QC report folder), `--allowMissingSamples` (do not fail when a sheet sample gets
no BAM).

## Sample sheet

The run's LabWare metadata sheet, the same file the germline DNA pipeline reads. Only RNA rows are
used: a row is RNA when it contains an IsoSeqX barcode. LabWare writes it in the proband slot:

```
rekv_npn_material_testlist_gender_IsoSeqX-bcNN-5p_intRef
```

- NPN = field 2.
- Barcode: `IsoSeqX[-_ ]bcNN` anywhere on the line, case-insensitive; must be 1-12.
- Pool: a `bcM000N` token on the row assigns the sample to that pool. LabWare does not write one
  today, so a run must then contain exactly one bcM pool.
- intRef (field after the barcode) is carried into the manifest only.

## Checks

At launch, before any task: at least one RNA row; at least one `hifi_reads/*.hifi_reads.bcM*.bam`;
no NPN twice; no pool/barcode twice; every named pool present; pools without RNA rows are reported
and skipped.

After lima: `rna_demultiplex/QC/<sheet>.demux_manifest.tsv` lists every sheet sample (OK or
NO_READS) and every lima barcode without a sheet sample (UNASSIGNED). NO_READS fails the run
unless `--allowMissingSamples`.

On success, `<runfolder>/kgvejle_rna_demultiplex_complete.txt` is written last. Automation
downstream should key on this file.

## Config

`lrsShared.config` -> `conf/rnaPreprocess.config` -> `finalize.config`, as in the other LRS
pipelines. The lrsShared folder is local to the server and is the only place server paths live;
this repository holds none. Keys this pipeline needs from lrsShared (v2 with the Kinnex block added), checked at launch:

- `kinnexMas8Adapters`, `kinnexBarcodes12plex`: Kinnex kit FASTAs
- `condaEnvs.skera`, `condaEnvs.lima`: conda envs (also `condaEnvs.somaticSummaryEnv`)
- `kinnexQcReportDir`: shared folder for the run QC report, copied to
  `<kinnexQcReportDir>/<runfolder name>/` (html) and `.../yaml_json/` (yaml, json); `false` switches the copy off

Nextflow reports: `rna_demultiplex/runInfo/<yyMMdd>_<sheet>/NextflowReports/`.

## Testing

`-stub` runs without skera/lima/samtools. `params.stubLimaBarcodes` (via `-params-file`) sets which
barcodes the lima stub "finds" per pool, e.g.

```yaml
stubLimaBarcodes:
  bcM0004: [IsoSeqX_bc01_5p, IsoSeqX_bc03_5p, IsoSeqX_bc07_5p]
```

No patient data, NPNs or real outputs in this repository.
