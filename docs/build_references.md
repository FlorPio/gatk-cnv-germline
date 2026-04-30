# Building the references Docker image

`build_references_image.sh` packages the **reference genome**, the **capture BED**, and a previously generated **Panel of Normals (PoN)** into a single Docker image (`germlinecnv-references:<tag>`). Once built, the image can be mounted in any environment and the pipeline run against the bundled paths — no need to re-share large files between machines or collaborators.

> The image only contains reference data. The pipeline itself is still launched with `nextflow run main.nf ...`.

## When to use it

- You finished generating a PoN (mode `pon`) and want to ship it together with the genome and BED so other people can run case mode without copying a dozen files.
- You want a reproducible, versioned snapshot of the references used by a given clinical batch (`hg38-v2.0`, `hg38-trusight-v1`, etc.).
- You want to avoid path drift between workstations / WSL / HPC.

If you only need the pipeline source code, you do **not** need this image — clone the repo and run Nextflow directly.

## Prerequisites

- A Linux or WSL shell (the script uses bash + `cp -r` + `mktemp`).
- Docker daemon running and your user with permission to build images.
- Local copies of:
  - `hg38.fa`, `hg38.fa.fai`, `hg38.dict`
  - the capture BED (e.g. `Regions.bed`)
  - a PoN directory containing the subdirectories `cnv_calls/`, `ploidy/`, `intervals/` (i.e. the `outdir` produced by a `--mode pon` run)

## What the script does

1. Creates a temporary staging directory (`mktemp -d`).
2. Copies into it:
   - `genome/` — FASTA + `.fai` + `.dict`
   - `bed/` — capture BED
   - `pon/cnv_calls/`, `pon/ploidy/`, `pon/intervals/` — PoN model + intervals
   - `assets/contig_ploidy_priors.tsv`, `assets/samplesheet.csv`
   - `docs/pon_documentation.md` (optional, if it exists)
   - The `Dockerfile.references` itself (renamed to `Dockerfile`)
3. Runs `docker build -t germlinecnv-references:<tag> <staging>`.
4. Cleans up the staging directory on exit.

## Configuring it for your setup

Open `build_references_image.sh` and adjust the variables at the top:

```bash
IMAGE_NAME="germlinecnv-references"
IMAGE_TAG="hg38-v2.0"             # bump on every new PoN / panel

GENOME_DIR="/home/hg38"           # contains hg38.fa, hg38.fa.fai, hg38.dict
BED_FILE="/home/Regions.bed"      # your capture BED
PON_DIR="/home/pon"               # outdir of a previous --mode pon run
ASSETS_DIR="assets"               # repo-relative
DOCS_DIR="docs"                   # repo-relative
DOCKERFILE="Dockerfile.references"
```

`PON_DIR` must contain at least:

```
${PON_DIR}/
├── cnv_calls/
│   └── cohort-cnv-model/cohort-model/
├── ploidy/
│   └── cohort-model/
└── intervals/
    └── intervals.interval_list
```

These are exactly the directories produced by `--mode pon` under `--outdir`.

## Building the image

From the repo root:

```bash
bash build_references_image.sh
```

Expected output (last lines):

```
SUCCESS: Image built as germlinecnv-references:hg38-v2.0

Image size:
  germlinecnv-references:hg38-v2.0  ~3.5GB
```

Verify the contents:

```bash
docker run --rm germlinecnv-references:hg38-v2.0
docker run --rm germlinecnv-references:hg38-v2.0 ls -la /references/genome/
docker run --rm germlinecnv-references:hg38-v2.0 ls /references/pon/cnv_calls/cohort-cnv-model/
```

## Using it in the pipeline

Once the image is built, you can run the pipeline pointing the parameters at the in-image paths. The simplest pattern is to launch a container that mounts your input/output dirs and runs Nextflow inside:

```bash
docker run --rm -it \
    -v $PWD:/work \
    -v /path/to/bams:/data/bams \
    -w /work \
    germlinecnv-references:hg38-v2.0 \
    nextflow run main.nf \
        --input /data/bams/samplesheet.csv \
        --mode case \
        --fasta       /references/genome/hg38.fa \
        --bed         /references/bed/Regions.bed \
        --pon_model   /references/pon/cnv_calls/cohort-cnv-model/cohort-model \
        --ploidy_model /references/pon/ploidy/cohort-model \
        --intervals   /references/pon/intervals/intervals.interval_list \
        --outdir results_case \
        -profile docker
```

Or, if you launch Nextflow on the host but want the references mounted inside processes, add a volume mount via `docker.runOptions` or use the image as a base for your processes (advanced).

### Pushing / sharing

```bash
docker tag  germlinecnv-references:hg38-v2.0  myregistry.example.org/germlinecnv-references:hg38-v2.0
docker push myregistry.example.org/germlinecnv-references:hg38-v2.0
```

Collaborators just `docker pull` the tag and run with the same `--fasta /references/...` paths.

## Versioning

Every time you regenerate the PoN or change the BED, bump `IMAGE_TAG` (e.g. `hg38-v2.1`, `hg38-trusight-2026-04`). Combined with the new `pon_manifest.json` (which contains the md5 of the FASTA, BED and intervals), this gives you a fully traceable reference snapshot — the manifest tells you *what* is inside, the tag tells you *which release* you are using.

## Troubleshooting

- **`docker build` fails copying `pon/...`** — `PON_DIR` does not have the expected `cnv_calls/`, `ploidy/`, `intervals/` subdirectories. Re-run a `--mode pon` job and point `PON_DIR` to its `outdir`.
- **Image is huge (>10 GB)** — you are bundling extra files. Check that `PON_DIR` is the PoN outdir and not the whole results directory (which also has counts and case calls).
- **Path mismatch when running** — paths inside the image always live under `/references/...`. The original host paths are irrelevant once the image is built.
- **Permissions on copied files** — the script uses `cp` (not `cp -p`); if you need to preserve ownership use `cp -p` and rebuild.
