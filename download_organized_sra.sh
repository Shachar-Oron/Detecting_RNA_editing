#!/usr/bin/env bash
set -euo pipefail

########################################
# CONFIG
########################################
THREADS="${THREADS:-16}"

# Base output for organized downloads
OUTPUT_BASE="${OUTPUT_BASE:-/private11/oronsha/Raw_Data_Benchmarking_RNAe_Tools}"

# Inputs (update if you keep them elsewhere)
SRR_LIST="${SRR_LIST:-/private11/oronsha/human_sra/SRR_Acc_List.txt}"
RUN_TABLE="${RUN_TABLE:-/private11/oronsha/human_sra/SraRunTable.csv}"

# Containers
declare -A CONTAINERS
CONTAINERS["sra-tools"]="${SRA_TOOLS_IMAGE:-ncbi/sra-tools:latest}"

########################################
# SANITY
########################################
die(){ echo "ERROR: $*" >&2; exit 1; }

command -v docker >/dev/null || die "Docker not found in PATH"
docker info >/dev/null 2>&1 || die "Docker daemon not accessible"

[ -f "$SRR_LIST" ] || die "SRR list not found: $SRR_LIST"
[ -f "$RUN_TABLE" ] || die "Run table not found: $RUN_TABLE"

mkdir -p "$OUTPUT_BASE"

########################################
# UTILS
########################################
docker_run() {
  local tool="$1"; shift
  local image="${CONTAINERS[$tool]:-}"
  [ -n "$image" ] || die "Unknown tool: $tool"
  docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v /private11/oronsha:/private11/oronsha \
  -w /private11/oronsha \
  "$image" "$@"
	
  #docker run --rm \
    #-u "$(id -u):$(id -g)" \
    #-v "$OUTPUT_BASE:$OUTPUT_BASE" \
    #-v "$(dirname "$SRR_LIST"):$(dirname "$SRR_LIST")" \
    #-v "$(dirname "$RUN_TABLE"):$(dirname "$RUN_TABLE")" \
    #"$image" "$@"
}

clean_dirname(){ echo "$1" | sed 's/[^A-Za-z0-9._-]/_/g; s/__*/_/g; s/^_//; s/_$//' ; }

# CSV line parser (handles simple quoted fields)
parse_csv_line() {
  local line="$1"; local -a out=(); local f=""; local q=0 c
  while IFS= read -r -n1 c; do
    if [[ "$c" == '"' ]]; then q=$((1-q))
    elif [[ "$c" == ',' && $q -eq 0 ]]; then out+=("$f"); f=""
    else f+="$c"
    fi
  done <<< "$line"
  out+=("$f")
  printf "%s\n" "${out[@]}"
}

########################################
# PARSE RUN TABLE → maps: RUN -> Organism, Assay
########################################
echo "[*] Parsing run table header…"
mapfile -t HDR < <(parse_csv_line "$(head -n1 "$RUN_TABLE")")
RUN_COL=-1; ORG_COL=-1; ASSAY_COL=-1
for i in "${!HDR[@]}"; do
  col="$(echo "${HDR[$i]}" | tr -d '"' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
  case "$col" in
    RUN|Run|run) RUN_COL=$i ;;
    Organism|organism|ORGANISM|ScientificName|scientific_name) ORG_COL=$i ;;
    "Assay Type"|Assay_Type|assay_type|ASSAY_TYPE|AssayType) ASSAY_COL=$i ;;
  esac
done
(( RUN_COL>=0 && ORG_COL>=0 && ASSAY_COL>=0 )) || die "Could not find required columns (Run/Organism/Assay) in $RUN_TABLE"

declare -A RUN2ORG RUN2ASSAY
lnum=0
while IFS= read -r line; do
  lnum=$((lnum+1)); [ $lnum -eq 1 ] && continue
  mapfile -t row < <(parse_csv_line "$line")
  [ ${#row[@]} -gt "$ASSAY_COL" ] || continue
  run="$(echo "${row[$RUN_COL]}" | tr -d '"' | xargs)"
  org="$(echo "${row[$ORG_COL]}" | tr -d '"' | xargs)"
  asy="$(echo "${row[$ASSAY_COL]}" | tr -d '"' | xargs)"
  [ -n "$run" ] && [ -n "$org" ] && [ -n "$asy" ] || continue
  RUN2ORG["$run"]="$org"
  RUN2ASSAY["$run"]="$asy"
done < "$RUN_TABLE"

echo "[*] Loaded ${#RUN2ORG[@]} metadata entries"

# Resolve the exact SRA path for a given SRR in a target dir (flat or nested)
resolve_sra_file() {
  local base="$1" acc="$2"
  local f1="$base/$acc.sra"
  local f2="$base/$acc/$acc.sra"
  if   [ -s "$f1" ]; then echo "$f1"
  elif [ -s "$f2" ]; then echo "$f2"
  else return 1
  fi
}

# Compress with pigz if available in host (faster), else gzip via container/host
compress_fastq_gz() {
  local fq
  if command -v pigz >/dev/null; then
    pigz -f -- "$@"
  else
    for fq in "$@"; do gzip -f -- "$fq"; done
  fi
}

########################################
# CORE: download + validate + convert
########################################
download_sra() {
  local SRR="$1" ORG="$2" ASSAY="$3"
  local ORG_CLEAN="$(clean_dirname "$ORG")"
  local ASSAY_CLEAN="$(clean_dirname "$ASSAY")"
  local TDIR="$OUTPUT_BASE/$ORG_CLEAN/$ASSAY_CLEAN"
  local FQDIR="$TDIR/fastq"
  mkdir -p "$FQDIR/tmp"

  local FQ1="$FQDIR/${SRR}_1.fastq.gz"
  local FQ2="$FQDIR/${SRR}_2.fastq.gz"

  echo ">>> $SRR"
  echo "    Organism: $ORG  →  $ORG_CLEAN"
  echo "    Assay   : $ASSAY →  $ASSAY_CLEAN"
  echo "    Target  : $TDIR"

  # If both mates already present and non-empty, skip
  if [ -s "$FQ1" ] && [ -s "$FQ2" ]; then
    echo "    ✓ FASTQs already present; skipping."
    return 0
  fi

  # Prefetch (resume) if the SRR.sra not present
  if ! resolve_sra_file "$TDIR" "$SRR" >/dev/null 2>&1; then
    echo "    … prefetching SRA (resume enabled)"
    docker_run "sra-tools" prefetch "$SRR" \
      --output-directory "$TDIR" --max-size 500G --progress --verify yes --resume yes
  fi

  # Re-resolve after prefetch
  local SRA
  if ! SRA="$(resolve_sra_file "$TDIR" "$SRR")"; then
    echo "    ✗ SRA not found for $SRR after prefetch"; return 1
  fi
  echo "    SRA: $SRA"

  # Validate SRA (remove partials if broken, then re-download once)
  if ! docker_run "sra-tools" vdb-validate "$SRA" >/dev/null 2>&1; then
    echo "    ✗ vdb-validate failed → removing partials and retrying prefetch"
    rm -f "$SRA" "${SRA}.prf" "${SRA}.tmp" 2>/dev/null || true
    docker_run "sra-tools" prefetch "$SRR" \
      --output-directory "$TDIR" --max-size 500G --progress --verify yes
    SRA="$(resolve_sra_file "$TDIR" "$SRR")" || { echo "    ✗ still missing $SRR.sra"; return 1; }
    docker_run "sra-tools" vdb-validate "$SRA" >/dev/null 2>&1 || { echo "    ✗ vdb-validate failed again"; return 1; }
  fi

  # Convert SRA → FASTQ (write temp within mounted disk to avoid container /tmp)
  echo "    … fasterq-dump ($THREADS threads)"
  docker_run "sra-tools" fasterq-dump "$SRA" \
    --outdir "$FQDIR" --split-files --threads "$THREADS" --temp "$FQDIR/tmp"

  # Compress any produced .fastq (handle PE or SE)
  shopt -s nullglob
  mapfile -t NEWFQ < <(find "$FQDIR" -maxdepth 1 -type f -name "${SRR}_*.fastq" -size +0c | sort)
  shopt -u nullglob
  if [ ${#NEWFQ[@]} -eq 0 ]; then
    echo "    ✗ fasterq-dump produced no FASTQ"; return 1
  fi

  echo "    … compressing FASTQ (${#NEWFQ[@]} files)"
  compress_fastq_gz "${NEWFQ[@]}"

  # Clean uncompressed remnants (if any)
  find "$FQDIR" -maxdepth 1 -type f -name "${SRR}_*.fastq" -delete 2>/dev/null || true
  rm -rf "$FQDIR/tmp" 2>/dev/null || true

  # Report PE/SE status
  if [ -s "$FQ1" ] && [ -s "$FQ2" ]; then
    echo "    ✓ Completed (paired-end)"
  elif [ -s "$FQ1" ] && [ ! -s "$FQ2" ]; then
    echo "    ✓ Completed (single-end)"
  else
    echo "    ⚠ Converted but mates not standard; check: $(ls -1 "$FQDIR"/${SRR}_*.fastq.gz || true)"
  fi

  echo ""
}

########################################
# MAIN LOOP
########################################
processed=0 failed=0 skipped=0
echo "[*] Starting SRR loop…"

while IFS= read -r SRR; do
  SRR="$(echo "$SRR" | tr -d '\r' | xargs)"
  [ -n "$SRR" ] || continue

  ORG="${RUN2ORG[$SRR]:-}"
  ASY="${RUN2ASSAY[$SRR]:-}"
  if [ -z "$ORG" ] || [ -z "$ASY" ]; then
    echo ">>> $SRR : ⚠ missing Organism/Assay in run table → skipping"
    skipped=$((skipped+1))
    continue
  fi

  if download_sra "$SRR" "$ORG" "$ASY"; then
    processed=$((processed+1))
  else
    failed=$((failed+1))
  fi
done < "$SRR_LIST"

echo "======================================"
echo "Processed : $processed"
echo "Failed    : $failed"
echo "Skipped   : $skipped"
echo "Output dir: $OUTPUT_BASE"
echo "======================================"




# #!/bin/bash

# # Thread count parameter (tune as needed)
# THREADS=16

# PIPELINE_HOME="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# OUTPUT_BASE="/private8/Projects/Ahiad/Raw_Data_Benchmarking_RNAe_Tools/"
# SRR_LIST="/private10/Projects/Ahiadchenzi/Big_Project_RNApipe_Diatom/rna_editing_full_pipe/Raw_Data_Benchmarking/sra_downloads/SRR_Acc_List whole_benchmark.txt"
# RUN_TABLE="/private10/Projects/Ahiadchenzi/Big_Project_RNApipe_Diatom/rna_editing_full_pipe/Raw_Data_Benchmarking/sra_downloads/SraRunTable whole_benchmark.csv"

# declare -A CONTAINERS
# CONTAINERS["sra-tools"]="ncbi/sra-tools:latest"

# die() { echo "ERROR: $*" >&2; exit 1; }

# docker_run() {
    # local tool="$1"
    # shift
    # local container="${CONTAINERS[$tool]}"
    # [[ -z "$container" ]] && die "Unknown tool for docker_run: $tool"
    # docker run --rm \
        # -u "$(id -u):$(id -g)" \
        # -v "${PIPELINE_HOME}:${PIPELINE_HOME}" \
        # -v "${OUTPUT_BASE}:${OUTPUT_BASE}" \
        # -v "$(dirname "$SRR_LIST"):$(dirname "$SRR_LIST")" \
        # -v "$(dirname "$RUN_TABLE"):$(dirname "$RUN_TABLE")" \
        # -w "${PIPELINE_HOME}" \
        # "$container" "$@"
# }

# [ ! -f "$SRR_LIST" ] || [ ! -f "$RUN_TABLE" ] && die "Required files not found!\n  SRR List: $SRR_LIST\n  Run Table: $RUN_TABLE"
# ! command -v docker &> /dev/null && die "Docker is not installed or not available in PATH"
# ! docker info &> /dev/null && die "Cannot access Docker. Please check Docker is running and you have permissions"

# mkdir -p "$OUTPUT_BASE"

# clean_dirname() {
    # echo "$1" | sed 's/[^A-Za-z0-9._-]/_/g' | sed 's/__*/_/g' | sed 's/^_//g' | sed 's/_$//g'
# }

# # Function to parse CSV line properly (handles quoted fields with commas)
# parse_csv_line() {
    # local line="$1"
    # local -a result=()
    # local field=""
    # local in_quotes=0
    # local char
    
    # while IFS= read -r -n1 char; do
        # if [[ "$char" == '"' ]]; then
            # if [[ $in_quotes -eq 0 ]]; then
                # in_quotes=1
            # else
                # in_quotes=0
            # fi
        # elif [[ "$char" == ',' && $in_quotes -eq 0 ]]; then
            # result+=("$field")
            # field=""
        # else
            # field+="$char"
        # fi
    # done <<< "$line"
    # result+=("$field")
    
    # # Print array elements
    # for element in "${result[@]}"; do
        # echo "$element"
    # done
# }

# # Parse header to find column indices
# echo "Parsing CSV header..."
# HEADER=$(head -n 1 "$RUN_TABLE")
# mapfile -t HEADER_ARRAY < <(parse_csv_line "$HEADER")

# RUN_COL=-1; ORGANISM_COL=-1; ASSAY_COL=-1

# echo "Available columns in CSV:"
# for i in "${!HEADER_ARRAY[@]}"; do
    # col_name=$(echo "${HEADER_ARRAY[$i]}" | tr -d '"' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    # echo "  Column $((i+1)): '$col_name'"
    
    # case "$col_name" in
        # "Run"|"run"|"RUN") 
            # RUN_COL=$((i+1))
            # echo "    -> Found RUN column"
            # ;;
        # "Organism"|"organism"|"ORGANISM"|"ScientificName"|"scientific_name") 
            # ORGANISM_COL=$((i+1))
            # echo "    -> Found ORGANISM column"
            # ;;
        # "Assay Type"|"Assay_Type"|"assay_type"|"ASSAY_TYPE"|"AssayType") 
            # ASSAY_COL=$((i+1))
            # echo "    -> Found ASSAY TYPE column"
            # ;;
    # esac
# done

# [ $RUN_COL -eq -1 ] && die "Could not find Run column"
# [ $ORGANISM_COL -eq -1 ] && die "Could not find Organism column"
# [ $ASSAY_COL -eq -1 ] && die "Could not find Assay Type column"

# echo ""
# echo "Column indices found:"
# echo "  Run: Column $RUN_COL"
# echo "  Organism: Column $ORGANISM_COL"
# echo "  Assay Type: Column $ASSAY_COL"
# echo ""

# declare -A ORGANISM_MAP
# declare -A ASSAY_MAP

# echo "Parsing CSV metadata..."
# line_count=0
# while IFS= read -r line; do
    # line_count=$((line_count + 1))
    # [ $line_count -eq 1 ] && continue  # Skip header
    
    # mapfile -t LINE < <(parse_csv_line "$line")
    
    # if [ ${#LINE[@]} -ge $ORGANISM_COL ]; then
        # RUN_ID=$(echo "${LINE[$((RUN_COL-1))]}" | tr -d '"' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        # ORGANISM=$(echo "${LINE[$((ORGANISM_COL-1))]}" | tr -d '"' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        # ASSAY=$(echo "${LINE[$((ASSAY_COL-1))]}" | tr -d '"' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        
        # if [ -n "$RUN_ID" ] && [ -n "$ORGANISM" ] && [ -n "$ASSAY" ]; then
            # ORGANISM_MAP["$RUN_ID"]="$ORGANISM"
            # ASSAY_MAP["$RUN_ID"]="$ASSAY"
            
            # # Debug: Print first few entries
            # if [ ${#ORGANISM_MAP[@]} -le 5 ]; then
                # echo "  Entry $RUN_ID -> Organism: '$ORGANISM', Assay: '$ASSAY'"
            # fi
        # fi
    # fi
# done < "$RUN_TABLE"

# echo ""
# echo "Loaded ${#ORGANISM_MAP[@]} metadata entries from CSV"
# [ ${#ORGANISM_MAP[@]} -eq 0 ] && die "No metadata entries loaded!"

# docker_run "sra-tools" prefetch --help &> /dev/null || die "Cannot run SRA tools in Docker."

# download_sra() {
    # local SRR_ID="$1"
    # local ORGANISM="$2"
    # local ASSAY="$3"
    # local CLEAN_ORGANISM=$(clean_dirname "$ORGANISM")
    # local CLEAN_ASSAY=$(clean_dirname "$ASSAY")
    # local TARGET_DIR="$OUTPUT_BASE/$CLEAN_ORGANISM/$CLEAN_ASSAY"
    # local FASTQ_DIR="$TARGET_DIR/fastq"
    # mkdir -p "$TARGET_DIR" "$FASTQ_DIR"

    # local SRA_FILE_PATH="$TARGET_DIR/$SRR_ID.sra"
    # local FASTQ_GZ_1="$FASTQ_DIR/${SRR_ID}_1.fastq.gz"

    # # Skip completed runs (SRA and FASTQ R1 must exist)
    # # if [ -f "$SRA_FILE_PATH" ] && [ -f "$FASTQ_GZ_1" ]; then
        # # echo "✓ $SRR_ID: SRA & FASTQ already present at $TARGET_DIR, skipping."
        # # return 0
    # # fi
	# # At the beginning of download_sra function, after defining paths:
	# local FASTQ_GZ_1="$FASTQ_DIR/${SRR_ID}_1.fastq.gz"
	# local FASTQ_GZ_2="$FASTQ_DIR/${SRR_ID}_2.fastq.gz"
	# local SRA_FILE_ACTUAL="$TARGET_DIR/$SRR_ID/$SRR_ID.sra"

	# # Update skip condition:
	# if [ -f "$SRA_FILE_ACTUAL" ] && [ -f "$FASTQ_GZ_1" ]; then
		# echo "✓ $SRR_ID: SRA & FASTQ already present at $TARGET_DIR, skipping."
		# return 0
	# fi

    # echo "Downloading $SRR_ID..."
    # echo "  Organism: $ORGANISM -> $CLEAN_ORGANISM"
    # echo "  Assay Type: $ASSAY -> $CLEAN_ASSAY"
    # echo "  Target: $TARGET_DIR"
    
    # if [ ! -f "$SRA_FILE_PATH" ]; then
        # if ! docker_run "sra-tools" prefetch "$SRR_ID" --output-directory "$TARGET_DIR" --max-size 500G --progress --verify yes; then
            # echo "✗ SRA download failed for $SRR_ID"
            # return 1
        # fi
    # else
        # echo "  SRA file already exists: $SRA_FILE_PATH"
    # fi

    # local SRA_FILE="$SRA_FILE_PATH"
    # [ -f "$SRA_FILE" ] || SRA_FILE=$(find "$TARGET_DIR" -name "*.sra" | head -1)
    # [ -f "$SRA_FILE" ] || { echo "✗ SRA file not found for $SRR_ID"; return 1; }

    # if [ ! -f "$FASTQ_GZ_1" ]; then
        # echo "  Converting SRA to FASTQ ($THREADS threads)..."
        # if docker_run "sra-tools" fasterq-dump "$SRA_FILE" --outdir "$FASTQ_DIR" --split-files --threads "$THREADS"; then
            # echo "  FASTQ conversion successful"
            # for fq in "$FASTQ_DIR"/*.fastq; do
                # [ -f "$fq" ] || continue
                # [ -f "${fq}.gz" ] && continue
                # gzip -f "$fq"
            # done
            # echo "  Compression complete"
        # else
            # echo "✗ FASTQ conversion failed for $SRR_ID"
            # return 1
        # fi
    # else
        # echo "  FASTQ files already exist, skipping conversion and compression"
    # fi

    # echo "✓ $SRR_ID: Completed"
    # echo ""
    # return 0
# }

# processed_count=0
# skipped_count=0
# failed_count=0

# echo "=== Starting Download Loop ==="
# while read -r SRR_ID; do
    # [ -z "$SRR_ID" ] && continue
    # SRR_ID=$(echo "$SRR_ID" | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    # ORGANISM="${ORGANISM_MAP[$SRR_ID]}"
    # ASSAY="${ASSAY_MAP[$SRR_ID]}"
    # if [ -z "$ORGANISM" ] || [ -z "$ASSAY" ]; then
        # echo "⚠ Missing metadata for $SRR_ID, skipping."
        # skipped_count=$((skipped_count+1))
        # continue
    # fi
    # if download_sra "$SRR_ID" "$ORGANISM" "$ASSAY"; then
        # processed_count=$((processed_count+1))
    # else
        # failed_count=$((failed_count+1))
    # fi
# done < "$SRR_LIST"

# echo "=== Processing Summary ==="
# echo "Processed: $processed_count"
# echo "Failed: $failed_count"
# echo "Skipped: $skipped_count"
# echo "Done!"