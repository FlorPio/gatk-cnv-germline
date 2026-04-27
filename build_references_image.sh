#!/bin/bash
# =============================================================================
# Build script for germlinecnv-references Docker image
#
# This script stages all reference files into a temporary build context and
# builds the Docker image. Must be run from a Linux/WSL environment where
# the reference file paths are accessible.
#
# Usage: bash build_references_image.sh
# =============================================================================

set -euo pipefail

IMAGE_NAME="germlinecnv-references"
IMAGE_TAG="hg38-v2.0"
FULL_TAG="${IMAGE_NAME}:${IMAGE_TAG}"

# Source paths (adjust if your paths differ)
GENOME_DIR="/home/hg38"
BED_FILE="/home/Regions.bed"
PON_DIR="/home/pon"
ASSETS_DIR="assets"
DOCS_DIR="docs"
DOCKERFILE="Dockerfile.references"

# Temporary staging directory
STAGING_DIR=$(mktemp -d)
echo "============================================"
echo "Building Docker image: ${FULL_TAG}"
echo "Staging directory: ${STAGING_DIR}"
echo "============================================"

cleanup() {
    echo ""
    echo "Cleaning up staging directory..."
    rm -rf "${STAGING_DIR}"
}
trap cleanup EXIT

# --- Stage genome reference ---
echo "[1/6] Copying genome reference (hg38)..."
mkdir -p "${STAGING_DIR}/genome"
cp "${GENOME_DIR}/hg38.fa"     "${STAGING_DIR}/genome/"
cp "${GENOME_DIR}/hg38.fa.fai" "${STAGING_DIR}/genome/"
cp "${GENOME_DIR}/hg38.dict"   "${STAGING_DIR}/genome/"

# --- Stage BED file ---
echo "[2/6] Copying BED file..."
mkdir -p "${STAGING_DIR}/bed"
cp "${BED_FILE}" "${STAGING_DIR}/bed/"

# --- Stage PoN  ---
echo "[3/6] Copying PoN model..."
mkdir -p "${STAGING_DIR}/pon"
cp -r "${PON_DIR}/cnv_calls"  "${STAGING_DIR}/pon/"
cp -r "${PON_DIR}/ploidy"     "${STAGING_DIR}/pon/"
cp -r "${PON_DIR}/intervals"  "${STAGING_DIR}/pon/"

# --- Stage assets ---
echo "[4/6] Copying pipeline assets..."
mkdir -p "${STAGING_DIR}/assets"
cp "${ASSETS_DIR}/contig_ploidy_priors.tsv"         "${STAGING_DIR}/assets/"
cp "${ASSETS_DIR}/samplesheet.csv"  "${STAGING_DIR}/assets/"

# --- Stage documentation ---
echo "[5/6] Copying documentation..."
mkdir -p "${STAGING_DIR}/docs"
cp "${DOCS_DIR}/pon_documentation.md" "${STAGING_DIR}/docs/"

# --- Copy Dockerfile ---
cp "${DOCKERFILE}" "${STAGING_DIR}/Dockerfile"

# --- Build image ---
echo "[6/6] Building Docker image..."
echo ""
docker build -t "${FULL_TAG}" "${STAGING_DIR}"

echo ""
echo "============================================"
echo "SUCCESS: Image built as ${FULL_TAG}"
echo "============================================"
echo ""
echo "Image size:"
docker image ls "${IMAGE_NAME}" --format "  {{.Repository}}:{{.Tag}}  {{.Size}}"
echo ""
echo "Verify with:"
echo "  docker run --rm ${FULL_TAG}"
echo "  docker run --rm ${FULL_TAG} ls -la /references/genome/"
echo "  docker run --rm ${FULL_TAG} ls -la /references/pon/cnv_calls/cohort-cnv-model/ | head"
echo ""
echo "Use in pipeline:"
echo "  --fasta /references/genome/hg38.fa"
echo "  --bed /references/bed/Regions.bed"
echo "  --pon_model /references/pon/cnv_calls/cohort-cnv-model"
echo "  --ploidy_model /references/pon/ploidy/cohort-model"
