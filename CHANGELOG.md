# Changelog

## 2.0.2 (2026-11-18)
Deployed: not deployed.

- condaEnvs.skera/lima back in conf/rnaPreprocess.config (built from the shared root
  params.condaLocal): pipeline-only envs do not belong in lrsShared. lrsShared needs only the
  Kinnex keys kinnexMas8Adapters, kinnexBarcodes12plex, kinnexQcReportDir.

## 2.0.1 (2026-11-18)
Deployed: not deployed.

- lrsShared pin v3 -> v2: the Kinnex block is added to v2 (additive, no existing v2 key changed)
  instead of a new version. Launch-check message updated.

## 2.0.0 (2026-11-17)
First version in git. Deployed: not deployed.

- Sample sheet: the run's LabWare metadata sheet replaces the two-column NPN/barcode sheet.
  RNA rows are those with an IsoSeqX barcode (LabWare puts it in the proband slot).
- Config: lrsShared -> conf/rnaPreprocess.config -> finalize.config. Replaces the old
  single shared config. All server paths (Kinnex FASTAs, skera/lima envs, QC report folder) come
  from lrsShared; none in this repository. Missing keys stop the run at launch.
- Pool-safe: samples are matched on pool (bcM) + 5' barcode, not barcode alone. Output names carry
  the bcM token: <npn>.<movie>.hifi_reads.<bcM>.lima.<bc5>.SMfixed.bam (still matches the
  downstream glob **/*.lima.*.bam). lima output prefix carries the bcM too.
- A pool sequenced on several SMRT cells gives one uBAM per movie per sample (combine, not join).
- Launch checks (see README) stop the run before any task.
- demux manifest (OK / NO_READS / UNASSIGNED); NO_READS fails the run unless --allowMissingSamples.
- Completion flag only for a successful run writing into the run folder (not with --outdir).
- QC report copied to <kinnexQcReportDir>/<runfolder name>/ (html) and yaml_json/ (yaml, json), as
  the old script did with summaryBase; not with --outdir.
- kinnex_run_summary: --mas-barcode now gets the pool barcode (was meta.barcode5, which was unset
  for pool BAMs); the per-pool NPN/barcode sheet is built from the LabWare sheet.
- skera split gets --num-threads (speed only).
- SM rewrite: header split on tabs; every @RG SM checked after reheader.
- Removed: params.outdir2 (RNASEQstaging), unused.
