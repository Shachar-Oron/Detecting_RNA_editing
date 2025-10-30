#!/usr/bin/env bash
###############################################################################
# rnae: RNA Editing Pipeline (Docker-enabled, tool-specific containers)
# Version: 3.5.0-docker
# Author: KranzlerLab
#
# Docker containers for each tool - safer and more reproducible
# Usage: ./rnae_docker.sh --pipeline-home /path/to/project --config config.env --steps redi,jacusa2 --threads 24
###############################################################################

set -Eeuo pipefail
export LC_ALL=C
IFS=$'\n\t'

# ----------------------------- Docker Setup ----------------------------------
# Container images (latest stable versions from BioContainers)
declare -A CONTAINERS=(
	["fastp"]="quay.io/biocontainers/fastp:0.18.0--hd28b015_0"
	["star"]="quay.io/biocontainers/star:2.7.11b--h5ca1c30_7"
	["samtools"]="quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"
	["bedtools"]="quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0"
	["bwa"]="quay.io/biocontainers/bwa:0.7.17--ha92aebf_3"
	["python"]="quay.io/biocontainers/python:3.9--h4d93c1b_0_cpython" #python:3.9.23-slim-trixie
	["reditools"]="quay.io/biocontainers/reditools3:3.4--pyhdfd78af_0"
	#["reditools"]="chiaenu/reditools3:3.4"
	["jacusa"]="quay.io/biocontainers/jacusa2:2.0.4--hdfd78af_0"
	["sprint"]="sprint:py2-0.1.8"
	["rseqc"]="quay.io/biocontainers/rseqc:5.0.3--py39hf95cd2a_0"
)


# Pipeline home directory (contains all data, references, outputs)
PIPELINE_HOME="${PIPELINE_HOME:-$PWD}"

# GTF Check Function
has_gtf(){
	[[ -n "${GENOME_GTF:-}" && -s "${GENOME_GTF}" ]]
}

# Docker run command builder with better error handling
# docker_run() {
	# local tool="$1"
	# shift
	# local container="${CONTAINERS[$tool]}"
	
	# # Check if container exists in our mapping
	# if [[ -z "$container" ]]; then
		# die "Unknown tool for docker_run: $tool"
	# fi
	
	# docker run --rm \
		# -u "$(id -u):$(id -g)" \
		# -v "${PIPELINE_HOME}:${PIPELINE_HOME}" \
		# -w "${PIPELINE_HOME}" \
		# "$container" "$@"
# }
docker_run() {
	local tool="$1"
	shift
	local container="${CONTAINERS[$tool]}"
	
	if [[ -z "$container" ]]; then
		die "Unknown tool for docker_run: $tool"
	fi
	
	# Build volume mount flags as an array
	local clean_pipeline_home="${PIPELINE_HOME%/}"
	local volume_flags=(-v "${clean_pipeline_home}:${clean_pipeline_home}")
	
	# Add EXTRA_MOUNTS if not empty (colon-separated paths)
	if [[ -n "${EXTRA_MOUNTS:-}" ]]; then
		IFS=':' read -ra MOUNTS <<< "$EXTRA_MOUNTS"
		for mount_path in "${MOUNTS[@]}"; do
			if [[ -n "$mount_path" ]]; then
				local clean_mount="${mount_path%/}"
				volume_flags+=(-v "${clean_mount}:${clean_mount}")
			fi
		done
	fi
	
	docker run --rm \
		-u "$(id -u):$(id -g)" \
		"${volume_flags[@]}" \
		-w "$PIPELINE_HOME" \
		"$container" \
		"$@"
}




# ----------------------------- Logging & Traps -------------------------------
log(){ echo "[$(date +'%F %T')] $*"; }
logv(){ if [[ "${VERBOSE:-0}" == "1" ]]; then log "[VERBOSE] $*"; fi; }
warn(){ log "[WARN] $*"; }
die(){ log "[ERROR] $*"; exit 1; }
error(){ log "[ERROR] $*" >&2; }

SCRIPT_START_TIME=$(date +%s)
cleanup(){
	local rc=$?
	local runtime=$(( $(date +%s) - SCRIPT_START_TIME ))
	log "Runtime: $((runtime/3600))h $(((runtime%3600)/60))m $((runtime%60))s"
	exit $rc
}
trap cleanup EXIT
trap 'rc=$?; [[ $rc -ne 0 ]] && log "[ERROR] line ${LINENO}: ${BASH_COMMAND} (exit $rc)"' ERR

run_cmd(){ if [[ "${DRY_RUN:-0}" -eq 1 ]]; then log "[DRY-RUN] $*"; else logv "Exec: $*"; "$@"; fi; }
run_sh(){ if [[ "${DRY_RUN:-0}" -eq 1 ]]; then log "[DRY-RUN] $*"; else logv "Exec(sh): $*"; bash -o pipefail -lc "$*"; fi; }
try_cmd(){ set +e; "$@"; local rc=$?; set -e; ((rc==0)) || warn "Command failed (continuing): $* [exit $rc]"; }

# ----------------------------- Built-in Defaults -----------------------------
THREADS="${THREADS:-16}"
RESUME="${RESUME:-1}"
FORCE="${FORCE:-0}"
DRY_RUN="${DRY_RUN:-0}"
VERBOSE="${VERBOSE:-0}"
OUT_DIR="${OUT_DIR:-$PIPELINE_HOME/rnae_out}"
RAW_DIR="${RAW_DIR:-$OUT_DIR/raw/RNA}"
RAW_PREFIX="${RAW_PREFIX:-}"
RAW_SUFFIX="${RAW_SUFFIX:-.fastq}"
RAW_DNA_DIR="${RAW_DNA_DIR:-}"
RAW_DNA_PATTERN="${RAW_DNA_PATTERN:-*.fastq*}"
LAYOUT="${LAYOUT:-auto}"
FILTERED_PREFIX="${FILTERED_PREFIX:-filtered_}"
BAM_SUFFIX_DEFAULT="${BAM_SUFFIX_DEFAULT:-_Aligned.sortedByCoord.out.bam}"
BAM_SUFFIX="${BAM_SUFFIX:-_Aligned.sortedByCoord.out.bam}"
UNMAPPED_STUB="${UNMAPPED_STUB:-Unmapped.out.mate}"
READ_LENGTH="${READ_LENGTH:-75}"
READ_LENGTH_AUTO_DETECT="${READ_LENGTH_AUTO_DETECT:-1}"

# Reference files
STAR_INDEX_DIR="${STAR_INDEX_DIR:-$OUT_DIR/indexes/STAR}"
BWA_INDEX_DIR="${BWA_INDEX_DIR:-$OUT_DIR/indexes/BWA}"
FASTA_INDEX_DIR="${FASTA_INDEX_DIR:-$OUT_DIR/indexes/FASTA}"
SPRINT_INDEX_DIR="${SPRINT_INDEX_DIR:-$OUT_DIR/indexes/SPRINT}"
GENOME_FA="${GENOME_FA:-${FASTA_INDEX_DIR}/genomic.fa}"
GENOME_GTF="${GENOME_GTF:-${FASTA_INDEX_DIR}/genomic.gtf}"
GENOME_GFF="${GENOME_GFF:-${FASTA_INDEX_DIR}/genomic.gff}"

# Annotation settings
ANNOTATE_REDI="${ANNOTATE_REDI:-1}"
ANNOTATOR_BIN="${ANNOTATOR_BIN:-enhanced_gene_annotator.py}"
ANNOTATOR_THREADS="${ANNOTATOR_THREADS:-16}"
ANNOTATOR_BATCH_SIZE="${ANNOTATOR_BATCH_SIZE:-10000}"
ANNOTATOR_CALC_REI="${ANNOTATOR_CALC_REI:-0}"

# Optional tools
REI_BIN="${REI_BIN:-RNAEditingIndex}"
JACUSA_JAR="${JACUSA_JAR:-}"
HE_RUNNER="${HE_RUNNER:-}"
HE_ANALYZE_ENABLE="${HE_ANALYZE_ENABLE:-1}"
HE_UNCOMPRESS="${HE_UNCOMPRESS:-gz}"
HE_FILES_SUFFIX="${HE_FILES_SUFFIX:-gz}"
# Strand handling (0=unstranded, 1=second-strand, 2=first-strand)
STRAND_LIB="${STRAND_LIB:-auto}"			  # used by REDItools analyze and JACUSA2
REDI_INDEX_STRAND="${REDI_INDEX_STRAND:-*}"	 # '*'=no strand filter in reditools index
REDI_WINDOW="${REDI_WINDOW:-100000}"	  # sensible default window for REDItools3
REDI_MIN_FREQ="${REDI_MIN_FREQ:-0.01}"	# default minimum editing frequency (10%)

# Metadata
GROUPS_CSV="${GROUPS_CSV:-}"
VIRUS_PATTERN="${VIRUS_PATTERN:-}"

# ----------------------------- SPRINT Settings --------------------------------
# Docker image tag (can be overridden with --sprint-image or env var)
SPRINT_IMAGE="${SPRINT_IMAGE:-sprint:py2-0.1.8}"
SPRINT_RP="${SPRINT_RP:-${FASTA_INDEX_DIR}/repeats.min5.bed}"		  # Repeat annotation file
SPRINT_CD="${SPRINT_CD:-}"	   # Cluster distance cutoff (-cd)
SPRINT_ARGS="${SPRINT_ARGS:-}"	  # Free-form extra args, e.g. "-csrg 5"
# SPRINT executable paths inside the container
SPRINT_BWA_PATH="/usr/local/bin/bwa"
SPRINT_SAMTOOLS_PATH="/usr/local/bin/samtools"

# ----------------------------- Utility Functions -----------------------------
samples_dir(){ echo "${OUT_DIR}/samples"; }
sample_dir(){ echo "$(samples_dir)/$1"; }
fastp_dir(){ echo "$(sample_dir "$1")/fastp"; }
star_dir(){ echo "$(sample_dir "$1")/star"; }
ensure_dir(){ mkdir -p "$1"; }

mate_of(){
	local r1="$1" cand
	local pairs=("_R1|_R2" ".R1|.R2" "_1|_2" ".1|.2")
	for p in "${pairs[@]}"; do
		local a="${p%%|*}" b="${p##*|}"
		cand="${r1//$a/$b}"
		if [[ "$cand" != "$r1" && -f "$cand" ]]; then
			echo "$cand"
			return
		fi
	done
	echo ""
}

decide_layout_for_file(){
	case "${LAYOUT^^}" in
		SE) echo SE;;
		PE) echo PE;;
		*) [[ -n "$(mate_of "$1")" ]] && echo PE || echo SE;;
	esac
}

sample_from_prefix_suffix(){
	local p="$1" b; b="$(basename "$p")"
	
	# Remove file extensions first
	b="${b%.gz}"
	b="${b%.fastq}"
	b="${b%.fq}"
	
	# Remove prefix and suffix if specified in config
	local ok_start=1 ok_end=1
	[[ -n "$RAW_PREFIX" && "${b}" != ${RAW_PREFIX}* ]] && ok_start=0
	[[ -n "$RAW_SUFFIX" && "${b}" != *${RAW_SUFFIX} ]] && ok_end=0
	if (( ok_start && ok_end )); then
		[[ -n "$RAW_PREFIX" ]] && b="${b#${RAW_PREFIX}}"
		[[ -n "$RAW_SUFFIX" ]] && b="${b%${RAW_SUFFIX}}"
	fi
	
	# Remove filtered prefix if present
	[[ -n "$FILTERED_PREFIX" && "$b" == ${FILTERED_PREFIX}* ]] && b="${b#${FILTERED_PREFIX}}"
	
	# CRITICAL FIX: Keep everything UP TO _R1/_R2, remove _R1/_R2 and everything after
	# Example: 09042025_1_S1_L001_R1_001 -> 09042025_1_S1_L001
	if [[ "$b" =~ ^(.+)_R[12](_.*)?$ ]]; then
		b="${BASH_REMATCH[1]}"
	elif [[ "$b" =~ ^(.+)\.R[12](\..*)?$ ]]; then
		b="${BASH_REMATCH[1]}"
	elif [[ "$b" =~ ^(.+)_[12](_.*)?$ ]]; then
		b="${BASH_REMATCH[1]}"
	elif [[ "$b" =~ ^(.+)\.[12](\..*)?$ ]]; then
		b="${BASH_REMATCH[1]}"
	fi
	
	echo "$b"
}

# Convert GTF to BED format for RSeQC (simple awk method)
ensure_genome_bed() {
	local gtf_file="${GENOME_GTF_LINK:-$GENOME_GTF}"
	# Create BED file in FASTA_INDEX_DIR (writable), not next to the GTF
	local gtf_basename="$(basename "${gtf_file%.gtf}")"
	local bed_file="${FASTA_INDEX_DIR}/${gtf_basename}.bed"
	
	# DEBUG LOGGING
	log "DEBUG: ensure_genome_bed called"
	log "DEBUG: GENOME_GTF_LINK=${GENOME_GTF_LINK:-NOT_SET}"
	log "DEBUG: GENOME_GTF=${GENOME_GTF}"
	log "DEBUG: gtf_file=$gtf_file"
	log "DEBUG: bed_file=$bed_file"
	log "DEBUG: FASTA_INDEX_DIR=$FASTA_INDEX_DIR"
	
	if [[ -s "$bed_file" && "$FORCE" -eq 0 ]]; then
		log "DEBUG: BED file already exists: $bed_file"
		echo "$bed_file"
		return 0
	fi
	
	if [[ ! -s "$gtf_file" ]]; then
		warn "GTF file not found or empty: $gtf_file"
		log "DEBUG: GTF check failed - file does not exist or is empty"
		return 1
	fi
	
	log "Converting GTF to BED format for strand detection..."
	log "DEBUG: GTF file exists and is readable"
	ensure_dir "$FASTA_INDEX_DIR"
	log "DEBUG: FASTA_INDEX_DIR created/verified"
	
	# Run awk directly on host (no Docker needed for this simple operation)
	awk -F'\t' '$3 == "exon" {print $1"\t"$4-1"\t"$5"\t.\t0\t"$7}' "$gtf_file" > "$bed_file" || {
		warn "Failed to convert GTF to BED"
		log "DEBUG: awk command failed"
		return 1
	}
	
	log "DEBUG: awk command completed"
	
	if [[ -s "$bed_file" ]]; then
		log "BED file created: $bed_file"
		log "DEBUG: BED file size: $(stat -c%s "$bed_file" 2>/dev/null || stat -f%z "$bed_file" 2>/dev/null) bytes"
		echo "$bed_file"
	else
		log "DEBUG: BED file is empty or does not exist after awk"
		return 1
	fi
}

# Infer strand orientation using RSeQC
infer_strandedness() {
	local bam="$1"
	local bed="$2"
	
	[[ ! -s "$bam" || ! -s "$bed" ]] && { echo ""; return 0; }
	
	# Write logs to stderr (>&2), NOT stdout
	log "Inferring strand orientation from: $(basename "$bam")" >&2
	
	local result
	result=$(docker_run rseqc infer_experiment.py -i "$bam" -r "$bed" 2>&1) || { echo ""; return 0; }
	
	local firststrand secondstrand
	firststrand=$(echo "$result" | grep "1+-,1-+,2++,2--" | grep -oP '[0-9.]+$')
	secondstrand=$(echo "$result" | grep "1++,1--,2+-,2-+" | grep -oP '[0-9.]+$')
	
	firststrand=${firststrand:-0}
	secondstrand=${secondstrand:-0}
	
	log "RSeQC results: fr-firststrand=${firststrand}, fr-secondstrand=${secondstrand}" >&2
	
	local first_pct second_pct
	first_pct=$(awk "BEGIN {printf \"%.0f\", $firststrand * 100}")
	second_pct=$(awk "BEGIN {printf \"%.0f\", $secondstrand * 100}")
	
	log "RSeQC percentages: fr-firststrand=${first_pct}%, fr-secondstrand=${second_pct}%" >&2
	
	if (( first_pct > 80 )); then
		log "Detected: fr-firststrand (>80%) → STRAND_LIB=2" >&2
		echo "2"
	elif (( second_pct > 80 )); then
		log "Detected: fr-secondstrand (>80%) → STRAND_LIB=1" >&2
		echo "1"
	else
		log "Detected: unstranded (~50/50) → STRAND_LIB=0" >&2
		echo "0"
	fi
}


bam_for_sample(){
	local s="$1"; local bam="$(star_dir "$s")/${s}${BAM_SUFFIX_DEFAULT}"
	[[ -f "$bam" ]] && echo "$bam" || echo ""
}

has_dna_inputs(){
	[[ -n "${RAW_DNA_DIR:-}" && -d "$RAW_DNA_DIR" ]] || return 1
	shopt -s nullglob
	local arr=( "$RAW_DNA_DIR"/$RAW_DNA_PATTERN )
	shopt -u nullglob
	(( ${#arr[@]} > 0 ))
}

get_read_length(){
	local sample_file="$1"
	local length_raw
	local length

	# Check file exists and is readable	 
	if [[ ! -f "$sample_file" ]]; then
		die "File not found: $sample_file"
	fi

	if [[ ! -r "$sample_file" ]]; then
		die "File not readable: $sample_file"
	fi

	# Get sequence line length
	if [[ "$sample_file" == *.gz ]]; then
		length_raw=$(zcat "$sample_file" 2>/dev/null | sed -n '2p' | wc -c)
	else
		length_raw=$(sed -n '2p' "$sample_file" | wc -c)
	fi

	# Clean up whitespace
	length_raw=$(echo "$length_raw" | tr -d ' \t\n\r')

	# Validate
	if ! [[ "$length_raw" =~ ^[0-9]+$ ]]; then
		die "Failed to get numeric read length from $sample_file. Got: '$length_raw'"
	fi

	# Calculate final length
	if [[ $length_raw -gt 1 ]]; then
		length=$((length_raw - 1))
	else
		die "Invalid read length: $length_raw (too short)"
	fi

	# CRITICAL: Only echo the result, no log statements in this function
	echo "$length"
}


auto_detect_read_length(){
	if [[ "${READ_LENGTH_AUTO_DETECT:-0}" -eq 1 ]]; then
		log "Attempting to auto-detect read length..."
		
		shopt -s nullglob
		local first_file=("${RAW_DIR}/${RAW_PREFIX}"*"${RAW_SUFFIX}")
		shopt -u nullglob
		
		if [[ ${#first_file[@]} -gt 0 ]]; then
			log "Found ${#first_file[@]} input files, using: $(basename "${first_file[0]}")"
			
			local detected_length
			if detected_length=$(get_read_length "${first_file[0]}"); then
				# Log AFTER we have the clean result
				log "Successfully auto-detected read length: $detected_length bases"
				READ_LENGTH="$detected_length"
			else
				warn "Auto-detection failed, using configured read length: $READ_LENGTH"
			fi
		else
			warn "No input files found for auto-detection"
			warn "Using configured read length: $READ_LENGTH"
		fi
	else
		log "Using configured read length: $READ_LENGTH"
	fi
}

# ----------------------------- DNA Merge Helpers -----------------------------
dna_is_r1(){
	local f="$1"
	[[ "$f" =~ (_R1[^/]*\.f(ast)?q(\.gz)?$|\.R1[^/]*\.f(ast)?q(\.gz)?$|_1[^/]*\.f(ast)?q(\.gz)?$|\.1[^/]*\.f(ast)?q(\.gz)?$) ]]
}

dna_is_r2(){
	local f="$1"
	[[ "$f" =~ (_R2[^/]*\.f(ast)?q(\.gz)?$|\.R2[^/]*\.f(ast)?q(\.gz)?$|_2[^/]*\.f(ast)?q(\.gz)?$|\.2[^/]*\.f(ast)?q(\.gz)?$) ]]
}

dna_natural_sort(){
	# Print arguments one per line and natural-sort them
	printf '%s\n' "$@" | sort -V
}

dna_find_fastqs(){
	# Finds DNA FASTQs under RAW_DNA_DIR matching RAW_DNA_PATTERN.
	# Emits absolute paths, one per line.
	shopt -s nullglob
	local arr=( "$RAW_DNA_DIR"/$RAW_DNA_PATTERN )
	shopt -u nullglob
	if (( ${#arr[@]} == 0 )); then
		warn "DNA merge: no files in RAW_DNA_DIR='${RAW_DNA_DIR}' with pattern '${RAW_DNA_PATTERN}'"
		return 1
	fi
	# Normalize to absolute paths
	local p
	for p in "${arr[@]}"; do
		[[ -f "$p" ]] && readlink -f "$p"
	done
}

dna_merge_stream(){
	# $1 = "R1" or "R2"
	# $2.. = list of files (gz or plain) to merge into stdout (as plain FASTQ)
	local mate="$1"; shift
	local f
	for f in "$@"; do
		if [[ "$f" == *.gz ]]; then
			zcat -- "$f"
		else
			cat -- "$f"
		fi
	done
}

ensure_dna_merged(){
    [[ -n "${RAW_DNA_DIR:-}" && -d "$RAW_DNA_DIR" ]] || { warn "DNA merge: RAW_DNA_DIR not set or missing"; return 1; }

    local outd="${OUT_DIR}/DNA_seq_bwa_aln"; ensure_dir "$outd"
    local out_r1="${outd}/DNA_R1.merged.fastq"
    local out_r2="${outd}/DNA_R2.merged.fastq"
    local qc="${outd}/merge.qc.txt"

    # Resume guard
    if [[ "${RESUME:-1}" -eq 1 && -s "$out_r1" && ( -s "$out_r2" || ! -e "$out_r2" ) && "${FORCE:-0}" -eq 0 ]]; then
        log "DNA merge: resume – found merged outputs"
        return 0
    fi

    # If forcing, clean previous outputs
    if [[ "${FORCE:-0}" -eq 1 ]]; then
        rm -f "$out_r1" "$out_r2" "$qc"
    fi

    # Discover files
    shopt -s nullglob
    local all_files=( "$RAW_DNA_DIR"/$RAW_DNA_PATTERN )
    shopt -u nullglob
    
    if (( ${#all_files[@]} == 0 )); then
        warn "DNA merge: no FASTQs discovered"
        return 1
    fi

    # Split into R1/R2 buckets WITH VALIDATION
    local -A r1_map r2_map
    local base_name
    
    for f in "${all_files[@]}"; do
        if dna_is_r1 "$f"; then
            # Extract base name (everything before _R1/_1)
            base_name=$(basename "$f" | sed -E 's/[._](R)?1[._].*//')
            r1_map["$base_name"]="$f"
        elif dna_is_r2 "$f"; then
            base_name=$(basename "$f" | sed -E 's/[._](R)?2[._].*//')
            r2_map["$base_name"]="$f"
        fi
    done

    # Verify pairing
    local unpaired=0
    for base in "${!r1_map[@]}"; do
        if [[ -z "${r2_map[$base]:-}" ]]; then
            error "DNA merge: R1 file has no R2 mate: ${r1_map[$base]}"
            unpaired=1
        fi
    done
    
    for base in "${!r2_map[@]}"; do
        if [[ -z "${r1_map[$base]:-}" ]]; then
            error "DNA merge: R2 file has no R1 mate: ${r2_map[$base]}"
            unpaired=1
        fi
    done
    
    if [[ $unpaired -eq 1 ]]; then
        error "DNA merge: Unpaired files detected - cannot proceed"
        return 1
    fi

    # Get sorted list of base names
    local sorted_bases=()
    while IFS= read -r base; do
        sorted_bases+=("$base")
    done < <(printf '%s\n' "${!r1_map[@]}" | sort -V)

    log "DNA merge: Found ${#sorted_bases[@]} paired file sets"
    
    # Merge in order
    for base in "${sorted_bases[@]}"; do
        local r1="${r1_map[$base]}"
        local r2="${r2_map[$base]}"
        
        log "DNA merge: Adding pair $(basename "$r1") + $(basename "$r2")"
        
        # Append R1
        if [[ "$r1" == *.gz ]]; then
            zcat "$r1" >> "$out_r1"
        else
            cat "$r1" >> "$out_r1"
        fi
        
        # Append R2
        if [[ "$r2" == *.gz ]]; then
            zcat "$r2" >> "$out_r2"
        else
            cat "$r2" >> "$out_r2"
        fi
    done

    # Validate merged files
    local r1_reads r2_reads
    r1_reads=$(($(wc -l < "$out_r1") / 4))
    r2_reads=$(($(wc -l < "$out_r2") / 4))
    
    if [[ $r1_reads -ne $r2_reads ]]; then
        error "DNA merge: Read count mismatch after merge!"
        error "  R1: $r1_reads reads"
        error "  R2: $r2_reads reads"
        return 1
    fi

    # QC report
    {
        echo "# DNA merge QC ($(date +'%F %T'))"
        echo "Merged_pairs: ${#sorted_bases[@]}"
        echo "Total_reads: $r1_reads"
        echo "Files_merged:"
        for base in "${sorted_bases[@]}"; do
            echo "  - $(basename "${r1_map[$base]}") + $(basename "${r2_map[$base]}")"
        done
    } > "$qc"

    log "DNA merge: Complete – $r1_reads paired reads"
}

ensure_dna_merged_OLD(){
	# Creates ${OUT_DIR}/DNA_seq_bwa_aln/DNA_R1.merged.fastq and (optional) R2
	# from files in RAW_DNA_DIR matching RAW_DNA_PATTERN.
	[[ -n "${RAW_DNA_DIR:-}" && -d "$RAW_DNA_DIR" ]] || { warn "DNA merge: RAW_DNA_DIR not set or missing"; return 1; }

	local outd="${OUT_DIR}/DNA_seq_bwa_aln"; ensure_dir "$outd"
	local out_r1="${outd}/DNA_R1.merged.fastq"
	local out_r2="${outd}/DNA_R2.merged.fastq"
	local qc="${outd}/merge.qc.txt"

	# Resume guard
	if [[ "${RESUME:-1}" -eq 1 && -s "$out_r1" && ( -s "$out_r2" || ! -e "$out_r2" ) && "${FORCE:-0}" -eq 0 ]]; then
		log "DNA merge: resume — found merged outputs (${out_r1} $( [[ -s "$out_r2" ]] && echo 'and R2' ))"
		return 0
	fi

	# If forcing, clean previous outputs
	if [[ "${FORCE:-0}" -eq 1 ]]; then
		rm -f "$out_r1" "$out_r2" "$qc"
	fi

	# Discover files
	mapfile -t all < <( dna_find_fastqs ) || return 1
	if (( ${#all[@]} == 0 )); then
		warn "DNA merge: no FASTQs discovered"
		return 1
	fi

	# Split into R1/R2 buckets
	local r1_list=() r2_list=() other_list=()
	local f
	for f in "${all[@]}"; do
		if dna_is_r1 "$f"; then
			r1_list+=("$f")
		elif dna_is_r2 "$f"; then
			r2_list+=("$f")
		else
			other_list+=("$f")
		fi
	done

	# If neither R1 nor R2 matched, but we have "other", treat as SE
	if (( ${#r1_list[@]}==0 && ${#r2_list[@]}==0 && ${#other_list[@]} > 0 )); then
		warn "DNA merge: could not classify mates by name; assuming SINGLE-END from ${#other_list[@]} files"
		r1_list=( "${other_list[@]}" )
		other_list=()
	fi

	if (( ${#r1_list[@]} == 0 )); then
		warn "DNA merge: no R1 files; skipping DNA merge"
		return 1
	fi

	# Natural sort within each mate
	if (( ${#r1_list[@]} > 1 )); then
		mapfile -t r1_list < <( dna_natural_sort "${r1_list[@]}" )
	fi
	if (( ${#r2_list[@]} > 1 )); then
		mapfile -t r2_list < <( dna_natural_sort "${r2_list[@]}" )
	fi

	log "DNA merge: R1 files = ${#r1_list[@]}$( (( ${#r2_list[@]} > 0 )) && echo "; R2 files = ${#r2_list[@]}")"
	if (( ${#r2_list[@]} > 0 )); then
		log "DNA merge: detected PAIRED-END DNA"
	else
		log "DNA merge: detected SINGLE-END DNA"
	fi

	# Merge to uncompressed .fastq (bwa mem reads plain fastq fast)
	log "DNA merge: writing ${out_r1}"
	dna_merge_stream "R1" "${r1_list[@]}" > "${out_r1}"

	if (( ${#r2_list[@]} > 0 )); then
		log "DNA merge: writing ${out_r2}"
		dna_merge_stream "R2" "${r2_list[@]}" > "${out_r2}"
	else
		rm -f "${out_r2}" 2>/dev/null || true
	fi

	# Quick QC (line counts; FASTQ lines should be multiple of 4)
	{
		echo "# DNA merge QC ($(date +'%F %T'))"
		echo "R1_input_files: ${#r1_list[@]}"
		printf 'R1_lines\t'; wc -l < "${out_r1}" || true
		if [[ -s "${out_r2}" ]]; then
			echo "R2_input_files: ${#r2_list[@]}"
			printf 'R2_lines\t'; wc -l < "${out_r2}" || true
		fi
	} > "${qc}"

	log "DNA merge: done → ${out_r1} $( [[ -s "${out_r2}" ]] && echo "and ${out_r2}" )"
}


# ----------------------------- Index Functions (GTF-Optional) ---------------
# ensure_fasta_index(){
	# ensure_dir "$FASTA_INDEX_DIR"
  
  # # Create fixed-name symlinks
	# ln -sf "$GENOME_FA" "$FASTA_INDEX_DIR/genomic.fa"
	# if [ -n "$GENOME_GTF" ] && [ -f "$GENOME_GTF" ]; then
	  # ln -sf "$GENOME_GTF" "$FASTA_INDEX_DIR/genomic.gtf"
	# fi
  
  # # Create fasta index if not present
	# if [ ! -f "$FASTA_INDEX_DIR/genomic.fa.fai" ]; then
	  # run_cmd docker_run samtools samtools faidx "$FASTA_INDEX_DIR/genomic.fa"
	# fi
# }
ensure_fasta_index(){
	ensure_dir "$FASTA_INDEX_DIR"
	# Create fixed-name symlinks (force overwrite)
	ln -sf "$GENOME_FA" "$FASTA_INDEX_DIR/genomic.fa"
	if [ -n "$GENOME_GTF" ] && [ -f "$GENOME_GTF" ]; then
		ln -sf "$GENOME_GTF" "$FASTA_INDEX_DIR/genomic.gtf"
	fi
	# Export symlink paths for use throughout pipeline
	export GENOME_FA_LINK="$FASTA_INDEX_DIR/genomic.fa"
	export GENOME_GTF_LINK="$FASTA_INDEX_DIR/genomic.gtf"
	# Create fasta index if not present
	if [ ! -f "$GENOME_FA_LINK.fai" ]; then
		run_cmd docker_run samtools samtools faidx "$GENOME_FA_LINK"
	fi
}


ensure_star_index(){
	ensure_dir "$STAR_INDEX_DIR"
	
	if [[ -s "${STAR_INDEX_DIR}/SAindex" || -s "${STAR_INDEX_DIR}/Genome" || -s "${STAR_INDEX_DIR}/SA" ]]; then
		log "STAR index present: $STAR_INDEX_DIR (skipping build)"
	else
		log "STAR index missing → building at: $STAR_INDEX_DIR"
		
		# Auto-detect read length if needed
		auto_detect_read_length
		log "Using read length: $READ_LENGTH (sjdbOverhang: $((READ_LENGTH - 1)))"
		
		if has_gtf; then
			log "Building STAR index WITH GTF annotation (memory optimized)"
			run_cmd docker_run star STAR --runMode genomeGenerate --genomeDir "$STAR_INDEX_DIR" \
				--genomeFastaFiles "${GENOME_FA_LINK:-$GENOME_FA}" --sjdbGTFfile "${GENOME_GTF_LINK:-$GENOME_GTF}" --sjdbOverhang $((READ_LENGTH - 1)) \
				--outFileNamePrefix "$STAR_INDEX_DIR/" \
				--runThreadN "${THREADS}" \
				--limitGenomeGenerateRAM 31000000000 \
				--genomeSAsparseD 3 \
				--genomeSAindexNbases 10 \
				--genomeChrBinNbits 12
		else
			log "Building STAR index WITHOUT GTF annotation (memory optimized)"
			run_cmd docker_run star STAR --runMode genomeGenerate --genomeDir "$STAR_INDEX_DIR" \
				--genomeFastaFiles "${GENOME_FA_LINK:-$GENOME_FA}" --sjdbOverhang $((READ_LENGTH - 1)) \
				--outFileNamePrefix "$STAR_INDEX_DIR/" \
				--runThreadN "${THREADS}" \
				--limitGenomeGenerateRAM 31000000000 \
				--genomeSAsparseD 3 \
				--genomeSAindexNbases 10 \
				--genomeChrBinNbits 12
		fi
	fi
}


ensure_bwa_index() {
  # Honor explicit prefix if valid
  if [[ -n "${BWA_INDEX_PREFIX-}" && -f "${BWA_INDEX_PREFIX}.bwt" ]]; then
	export BWA_INDEX_PREFIX
	log "BWA index: using preset BWA_INDEX_PREFIX='${BWA_INDEX_PREFIX}'"
	return 0
  fi

  # Prefer ${BWA_INDEX_DIR} if provided, else default to .../indexes/BWA
  local bwa_dir="${BWA_INDEX_DIR:-$(dirname "${FASTA_INDEX_DIR}")/BWA}"
  if [[ ! -d "$bwa_dir" ]]; then
	error "BWA index: expected directory not found: '$bwa_dir'"
	error "Create it or pass --bwa-index-dir / --bwa-index-prefix."
	return 1
  fi

  local first_bwt
  first_bwt="$(ls -1 "${bwa_dir}"/*.bwt 2>/dev/null | head -n1 || true)"
  if [[ -z "$first_bwt" ]]; then
	# Auto-build if the reference is available
	local ref="${GENOME_FA_LINK:-${FASTA_INDEX_DIR}/genomic.fa}"
	if [[ -s "$ref" ]]; then
	  log "BWA index: not found in '$bwa_dir' → building from ${ref}"
	  ensure_dir "$bwa_dir"
	  run_cmd docker_run bwa bwa index -p "${bwa_dir}/genomic.fa" "$ref"
	  first_bwt="${bwa_dir}/genomic.fa.bwt"
	else
	  error "BWA index: no .bwt and missing reference '$ref'"
	  return 1
	fi
  fi

  BWA_INDEX_PREFIX="${first_bwt%.bwt}"
  export BWA_INDEX_PREFIX
  log "BWA index: using BWA_INDEX_PREFIX='${BWA_INDEX_PREFIX}'"
}


JACUSA_REF=""
ensure_jacusa_ref(){
	ensure_dir "$FASTA_INDEX_DIR"
	# Use symlink if available, otherwise fall back to original GENOME_FA
	local genome_ref="${GENOME_FA_LINK:-$GENOME_FA}"
	local prefix="$(basename "${genome_ref%.*}")"
	local clean="${FASTA_INDEX_DIR}/${prefix}.jacusa.clean.fna"
	if awk 'BEGIN{IGNORECASE=1} /^>/ {next} { s=$0; gsub(/[ACGTN]/,"",s); if(length(s)>0){ bad=1; exit } } END{exit(bad?1:0)}' "$genome_ref"; then
		JACUSA_REF="${FASTA_INDEX_DIR}/${prefix}.fna"
		[[ -f "$JACUSA_REF" ]] || ln -sf "$genome_ref" "$JACUSA_REF"
		[[ -f "${JACUSA_REF}.fai" ]] || run_cmd docker_run samtools samtools faidx "$JACUSA_REF"
	else
		awk 'BEGIN{IGNORECASE=1} /^>/ {print; next} { t=toupper($0); gsub(/[^ACGTN]/,"N",t); print t }' "$genome_ref" > "$clean"
		run_cmd docker_run samtools samtools faidx "$clean"
		JACUSA_REF="$clean"
	fi
}

ensure_clean_ref() {
  # Create cleaned reference for REDItools (only if ambiguous nucleotides exist)
  ensure_dir "$FASTA_INDEX_DIR"
  
  # Use the symlink in FASTA_INDEX_DIR, not the original GENOME_FA
  local genome_ref="${GENOME_FA_LINK:-$FASTA_INDEX_DIR/genomic.fa}"
  local prefix=$(basename "$genome_ref" .fa)
  local clean_ref="$FASTA_INDEX_DIR/${prefix}.clean.fa"
  
  # Check if reference contains ambiguous codes (check entire file)
  if grep -v "^>" "$genome_ref" | grep -q -i '[WRYKMSBDHV]'; then
	log "REDItools: reference contains ambiguous IUPAC codes, creating cleaned version" >&2
	
	if [[ "$FORCE" -eq 1 ]] || [[ ! -f "$clean_ref" ]]; then
	  awk '
		BEGIN {IGNORECASE=1}
		/^>/ {print; next}
		{
		  t = toupper($0)
		  gsub(/[WRYKMSBDHV]/, "N", t)
		  print t
		}
	  ' "$genome_ref" > "$clean_ref"
	  
	  run_cmd docker_run samtools samtools faidx "$clean_ref" >&2
	  log "REDItools: cleaned reference created at $clean_ref" >&2
	else
	  log "REDItools: using existing cleaned reference" >&2
	fi
	
	echo "$clean_ref"
  else
	log "REDItools: reference is clean (no ambiguous codes), using original" >&2
	echo "$genome_ref"	# Return the symlink path, not GENOME_FA
  fi
}

# ----------------------------- Pipeline Steps (GTF-Optional) ----------------
_fastp_pair(){
	local in1="$1" in2="$2" out1="$3" out2="$4" rep="$5"
	run_cmd docker_run fastp fastp -i "$in1" -I "$in2" -o "$out1" -O "$out2" --thread "${THREADS}" \
		--html "${rep}.html" --json "${rep}.json"
}

_fastp_single(){
	local in1="$1" out1="$2" rep="$3"
	run_cmd docker_run fastp fastp -i "$in1" -o "$out1" --thread "${THREADS}" \
		--html "${rep}.html" --json "${rep}.json"
}

_fastp_on_rna(){
	shopt -s nullglob
	declare -A processed_samples

	# Search for BOTH _1 and _R1 patterns to support different naming conventions
	local all_r1_files=()
	
	# Collect files matching either pattern
	for pattern in "${RAW_DIR}/${RAW_PREFIX}"*_1"${RAW_SUFFIX}" "${RAW_DIR}/${RAW_PREFIX}"*_R1*"${RAW_SUFFIX}"; do
		for file in $pattern; do
			[[ -f "$file" ]] && all_r1_files+=("$file")
		done
	done

	log "Found ${#all_r1_files[@]} R1/1 files to process"

	for r1_file in "${all_r1_files[@]}"; do
		# Get clean sample name
		local sample; sample="$(sample_from_prefix_suffix "$r1_file")"
		
		# Skip if already processed
		[[ "${processed_samples[$sample]:-}" == "1" ]] && continue
		processed_samples[$sample]="1"
		
		# Create sample directory
		local sdir; sdir="$(samples_dir)/$sample"
		local fastp_dir="$sdir/fastp"
		ensure_dir "$fastp_dir"
		
		# Set up output files
		local out1="${fastp_dir}/${FILTERED_PREFIX}$(basename "$r1_file")"
		local report="${fastp_dir}/${FILTERED_PREFIX}${sample}"
		
		# Use mate_of to automatically find the mate (works with both _1/_2 and _R1/_R2)
		local r2_file; r2_file="$(mate_of "$r1_file")"
		
		if [[ -n "$r2_file" && -f "$r2_file" ]]; then
			# Paired-end processing
			local out2="${fastp_dir}/${FILTERED_PREFIX}$(basename "$r2_file")"
			if [[ "$RESUME" -eq 1 && -s "$out1" && -s "$out2" && "$FORCE" -eq 0 ]]; then
				log "fastp(RNA): resume $sample (PE)"; continue
			fi
			log "fastp(RNA): processing $sample (PE: $(basename "$r1_file") + $(basename "$r2_file"))"
			_fastp_pair "$r1_file" "$r2_file" "$out1" "$out2" "$report"
		else
			# Single-end fallback
			if [[ "$RESUME" -eq 1 && -s "$out1" && "$FORCE" -eq 0 ]]; then
				log "fastp(RNA): resume $sample (SE)"; continue
			fi
			warn "No mate found for $r1_file → running single-end"
			log "fastp(RNA): processing $sample (SE: $(basename "$r1_file"))"
			_fastp_single "$r1_file" "$out1" "$report"
		fi
	done
	
	shopt -u nullglob
}

step_fastp(){
	ensure_dir "$(samples_dir)"
	_fastp_on_rna
	
	# DNA QC is optional
	if has_dna_inputs; then
		ensure_dir "${OUT_DIR}/DNA_fastp"
		shopt -s nullglob
		local fq
		for fq in "$RAW_DNA_DIR"/${RAW_DNA_PATTERN}; do
			[[ -f "$fq" ]] || continue
			local layout; layout="$(decide_layout_for_file "$fq")"
			local sample; sample="$(sample_from_prefix_suffix "$fq")"
			local sdir="${OUT_DIR}/DNA_fastp/${sample}"; ensure_dir "$sdir"
			local out1="${sdir}/${FILTERED_PREFIX}$(basename "$fq")"
			
			if [[ "$layout" == "PE" && "$fq" =~ (_R1|\.R1|_1|\.1) ]]; then
				local mate; mate="$(mate_of "$fq")"; [[ -n "$mate" ]] || { warn "Missing DNA R2 for $fq"; continue; }
				local out2="${sdir}/${FILTERED_PREFIX}$(basename "$mate")"
				[[ "$RESUME" -eq 1 && -s "$out1" && -s "$out2" && "$FORCE" -eq 0 ]] || \
					_fastp_pair "$fq" "$mate" "$out1" "$out2" "${sdir}/${FILTERED_PREFIX}${sample}"
			else
				[[ "$RESUME" -eq 1 && -s "$out1" && "$FORCE" -eq 0 ]] || \
					_fastp_single "$fq" "$out1" "${sdir}/${FILTERED_PREFIX}${sample}"
			fi
		done
		shopt -u nullglob
	fi
}

step_star_index(){ ensure_star_index; }
step_star_align(){
	log "STAR-align: ensuring STAR index (will skip if already present)"
	step_star_index

	if ! compgen -G "${RAW_DIR}/${RAW_PREFIX}*${RAW_SUFFIX}" >/dev/null; then
		die "STAR: no raw RNA files matching '${RAW_PREFIX}*${RAW_SUFFIX}' in ${RAW_DIR}"
	fi

	if ! compgen -G "$(samples_dir)/*/fastp/${FILTERED_PREFIX}*" >/dev/null; then
		log "STAR: no filtered reads -> running fastp"
		step_fastp
	fi

	shopt -s nullglob

	for sample_dir in "$(samples_dir)"/*; do
		[[ -d "$sample_dir" ]] || continue

		local sample; sample="$(basename "$sample_dir")"

		if [[ "$sample" =~ _R[12]_ ]]; then
			warn "STAR: Skipping directory with _R1/_R2 in name: $sample"
			warn "This suggests a problem with sample directory naming"
			continue
		fi

		local star_output_dir; star_output_dir="$(star_dir "$sample")"
		ensure_dir "$star_output_dir"
		
		# Enhanced skip check - verify both BAM and log files exist and are non-empty
		local bam="${star_output_dir}/${sample}${BAM_SUFFIX_DEFAULT}"
		local log_final="${star_output_dir}/${sample}_Log.final.out"
		local bai="${bam}.bai"

		if [[ "$RESUME" -eq 1 && -s "$bam" && -s "$log_final" && -s "$bai" && "$FORCE" -eq 0 ]]; then
			log "STAR: resume $sample (BAM, index, and log complete)"; continue
		fi

		# CRITICAL FIX: Initialize both variables with empty strings
		local r1_file=""
		local r2_file=""
		
		# Look for BOTH _1 and _R1 patterns (your files use _1)
		shopt -s nullglob
		local r1_files=("$sample_dir"/fastp/${FILTERED_PREFIX}*_1*.fastq* "$sample_dir"/fastp/${FILTERED_PREFIX}*_R1*.fastq*)
		shopt -u nullglob

		if [[ ${#r1_files[@]} -gt 0 ]]; then
			r1_file="${r1_files[0]}"
			r2_file="$(mate_of "$r1_file")"
		else
			# Fallback: find any fastq file in the fastp directory
			shopt -s nullglob
			local all_files=("$sample_dir"/fastp/${FILTERED_PREFIX}*.fastq*)
			shopt -u nullglob
			if [[ ${#all_files[@]} -gt 0 ]]; then
				r1_file="${all_files[0]}"
				r2_file="$(mate_of "$r1_file")"	 # Try to find mate even for fallback
			fi
		fi

		[[ -z "$r1_file" ]] && { warn "STAR: No filtered files found for $sample"; continue; }

		local read_cmd=()
		[[ "$r1_file" == *.gz ]] && read_cmd=(--readFilesCommand zcat)

		local star_opts=(
			--genomeDir "$STAR_INDEX_DIR"
			--runThreadN "$THREADS"
			--twopassMode Basic
			--outReadsUnmapped Fastx
			--outSAMtype BAM SortedByCoordinate
			--outSAMattributes MD NM
			--outFileNamePrefix "${star_output_dir}/${sample}_"
			--limitBAMsortRAM 5000000000
		)

		if has_gtf; then
			star_opts+=(--quantMode GeneCounts)
			log "STAR: aligning $sample with gene counting enabled"
		else
			log "STAR: aligning $sample without gene counting (no GTF)"
		fi

		if [[ -n "$r2_file" && -f "$r2_file" ]]; then
			log "STAR: processing $sample (PE: $(basename "$r1_file") + $(basename "$r2_file"))"
			run_cmd docker_run star STAR "${star_opts[@]}" "${read_cmd[@]}" --readFilesIn "$r1_file" "$r2_file"
		else
			log "STAR: processing $sample (SE: $(basename "$r1_file"))"
			run_cmd docker_run star STAR "${star_opts[@]}" "${read_cmd[@]}" --readFilesIn "$r1_file"
		fi

		run_cmd docker_run samtools samtools index "$bam"
	done
	shopt -u nullglob
}

step_detect_strand(){
	# Skip if already explicitly set (not auto)
	if [[ -n "${STRAND_LIB:-}" ]] && [[ "${STRAND_LIB}" != "auto" ]]; then
		log "STRAND_LIB already set to ${STRAND_LIB}, skipping auto-detection"
		return 0
	fi
	
	log "Auto-detecting strand orientation..."
	ensure_fasta_index	   
	# Find first available BAM
	shopt -s nullglob
	local first_bam
	first_bam=$(find "$(samples_dir)" -name "*${BAM_SUFFIX}" -type f 2>/dev/null | head -n 1)
	shopt -u nullglob
	
	if [[ -z "$first_bam" || ! -s "$first_bam" ]]; then
		warn "No BAM files found. Run 'star-align' step first."
		warn "Using default STRAND_LIB=2"
		STRAND_LIB=2
		return 0
	fi
	
	# Create BED file from GTF
	local bed_file
	#bed_file=$(ensure_genome_bed) || bed_file=""
	bed_file=$(ensure_genome_bed 2>&1 | tail -n 1)
	
	if [[ -z "$bed_file" || ! -s "$bed_file" ]]; then
		warn "Cannot create BED file from GTF, using default STRAND_LIB=2"
		STRAND_LIB=2
		return 0
	fi
	
	# Run detection
	local detected_strand
	detected_strand=$(infer_strandedness "$first_bam" "$bed_file" 2>/dev/null) || detected_strand=""
	
	if [[ -n "$detected_strand" ]]; then
		STRAND_LIB="$detected_strand"
		log "Successfully detected STRAND_LIB=${STRAND_LIB}"
		
		# Write to a cache file so subsequent runs can use it
		local cache_file="${OUT_DIR}/strand_orientation.cache"
		echo "STRAND_LIB=${STRAND_LIB}" > "$cache_file"
		log "Cached result to: $cache_file"
	else
		warn "Auto-detection failed, using default STRAND_LIB=2"
		STRAND_LIB=2
	fi
}

step_dna_align() {
  ensure_fasta_index
  ensure_dir "${OUT_DIR}/DNA_seq_bwa_aln"
  
  # Use the symlinked FASTA from FASTA_INDEX_DIR instead of original GENOME_FA
  local GENOME_FA_LOCAL="${GENOME_FA_LINK:-${FASTA_INDEX_DIR}/genomic.fa}"

  export GENOME_FA_LOCAL
  
  ensure_bwa_index || return 1

  local outd="${OUT_DIR}/DNA_seq_bwa_aln"
  local R1="${outd}/DNA_R1.merged.fastq"
  local R2="${outd}/DNA_R2.merged.fastq"
  local SORTED="${outd}/dna.sorted.bam"
  local FINAL="${outd}/DNAseq.bam"
  local BWA_LOG="${outd}/bwa.stderr.log"

  # Auto-merge DNA files if merged FASTQs are missing but RAW_DNA_DIR has reads
  if [[ ! -s "$R1" ]]; then
	if has_dna_inputs; then
	  log "DNA alignment: merged FASTQs not found → auto-merging from RAW_DNA_DIR"
	  ensure_dna_merged || { error "DNA alignment: auto-merge failed"; return 1; }
	fi
  fi

  # Input validation checks
  if [[ ! -s "$R1" ]]; then 
	error "DNA alignment: missing '$R1'"
	return 1
  fi
  
  if [[ ! -s "$GENOME_FA_LOCAL" ]]; then 
	error "DNA alignment: missing GENOME_FA_LOCAL '$GENOME_FA_LOCAL'"
	return 1
  fi
  
  if [[ ! -f "${BWA_INDEX_PREFIX}.bwt" ]]; then 
	error "DNA alignment: missing BWA index '${BWA_INDEX_PREFIX}.bwt'"
	return 1
  fi

  # Resume guard
  if [[ "${RESUME:-1}" -eq 1 && -s "$FINAL" && "${FORCE:-0}" -eq 0 ]]; then
	log "DNA alignment: resume — ${FINAL}"
	return 0
  fi

  # Log configuration
  log "DNA alignment: BWA_INDEX_PREFIX='${BWA_INDEX_PREFIX}'"
  log "DNA alignment: GENOME_FA_LOCAL='${GENOME_FA_LOCAL}'"
  log "DNA alignment: R1='${R1}' R2='${R2:-}'"
  log "DNA alignment: streaming bwa→samtools sort …"

  # Stream: bwa mem → samtools sort
  # Capture bwa stderr to ${BWA_LOG} for troubleshooting
  run_cmd bash -lc "
	set -euo pipefail
	: > '${BWA_LOG}'
	docker run --rm -u \$(id -u):\$(id -g) \\
	  -v '${PIPELINE_HOME}:${PIPELINE_HOME}' -w '${PIPELINE_HOME}' ${CONTAINERS[bwa]} \\
	  bash -lc \"set -euo pipefail; \\
		if [[ -s '${R2}' ]]; then \\
		  bwa mem -t ${THREADS:-8} '${BWA_INDEX_PREFIX}' '${R1}' '${R2}'; \\
		else \\
		  bwa mem -t ${THREADS:-8} '${BWA_INDEX_PREFIX}' '${R1}'; \\
		fi\" 2> '${BWA_LOG}' \\
	| docker run -i --rm -u \$(id -u):\$(id -g) \\
	  -v '${PIPELINE_HOME}:${PIPELINE_HOME}' -w '${PIPELINE_HOME}' ${CONTAINERS[samtools]} \\
	  samtools sort -O BAM -o '${SORTED}' -
  "

  # Sanity check the sorted BAM
  if ! docker_run samtools samtools quickcheck -v "${SORTED}"; then
	error "DNA alignment: '${SORTED}' failed samtools quickcheck. See ${BWA_LOG} for bwa errors."
	return 1
  fi

  # Add MD/NM tags + index using the symlinked FASTA
  log "DNA alignment: calmd → ${FINAL}"
  run_cmd docker_run samtools bash -lc "
	set -euo pipefail
	samtools calmd -b '${SORTED}' '${GENOME_FA_LOCAL}' > '${FINAL}'
	
	# Validate BAM integrity before indexing
	if ! samtools quickcheck -v '${FINAL}'; then
	  echo 'ERROR: BAM file ${FINAL} is corrupted' >&2
	  exit 1
	fi
	
	samtools index '${FINAL}'
  "

  # Verify index was created successfully
  if [[ ! -s "${FINAL}.bai" ]]; then
	error "DNA alignment: BAM index creation failed for ${FINAL}"
	return 1
  fi
  log "DNA alignment: index verified → ${FINAL}.bai"

  # Final check + summary
  docker_run samtools samtools view -H "${FINAL}" | head -n 5
  docker_run samtools samtools flagstat "${FINAL}" | head -n 12
  log "DNA alignment: done → ${FINAL}"
}

# Helper function to filter annotated output by DNA genotype
filter_dna_mismatches() {
	local input_file="$1"
	local output_file="$2"
	
	awk -F'\t' 'BEGIN {OFS="\t"}
	NR==1 {
		# Find gAllSubs column
		for(i=1; i<=NF; i++) {
			if($i == "gAllSubs") gallsubs_col = i
		}
		if(!gallsubs_col) {
			print "ERROR: gAllSubs column not found in annotated file" > "/dev/stderr"
			exit 1
		}
		print
		next
	}
	{
		# Keep only rows where gAllSubs is "-" (no DNA variant)
		if($gallsubs_col == "-") {
			print
		}
	}' "$input_file" > "$output_file"
}


step_reditools(){
	ensure_fasta_index
	shopt -s nullglob

	# -------------------- Clean reference & params --------------------
	local REF_TO_USE; REF_TO_USE=$(ensure_clean_ref)
	local SLIB="${STRAND_LIB:-2}"
	local RT_WINDOW="${REDI_WINDOW:-100000}"
	local REPEAT_MINLEN="${REPEATS_MINLEN:-5}"
	local MIN_FREQ="${REDI_MIN_FREQ:-0.01}"
	log "REDItools: min frequency cutoff = ${MIN_FREQ}"

	# Files/dirs
	REPEATS_BED="${FASTA_INDEX_DIR}/repeats.min${REPEAT_MINLEN}.bed"
	export REPEATS_BED

	# -------------------- Extract CDS regions --------------------
	local CDS_BED="${FASTA_INDEX_DIR}/cds_regions.bed"
	local ENABLE_CDS_FILTER=1  # Set to 0 to disable CDS filtering
	
	if [[ "$ENABLE_CDS_FILTER" -eq 1 ]]; then
		if has_gtf; then
			if [[ "$FORCE" -eq 1 || ! -s "$CDS_BED" ]]; then
				log "REDItools: extracting CDS regions from GTF"
				docker_run bedtools bash -lc "
					awk -F'\t' '\$3 == \"CDS\" {
						print \$1\"\t\"\$4-1\"\t\"\$5\"\t.\t0\t\"\$7
					}' '${GENOME_GTF_LINK:-$GENOME_GTF}' | \
					sort -k1,1 -k2,2n > '$CDS_BED'
				"
				if [[ -s "$CDS_BED" ]]; then
					local cds_count=$(wc -l < "$CDS_BED")
					log "REDItools: extracted $cds_count CDS regions"
				else
					warn "REDItools: CDS extraction produced empty file"
					ENABLE_CDS_FILTER=0
				fi
			else
				log "REDItools: using cached CDS regions -> $CDS_BED"
			fi
		else
			warn "REDItools: GTF not available, CDS filtering disabled"
			ENABLE_CDS_FILTER=0
		fi
	fi

	# -------------------- Repeats BED (always) --------------------
	if [[ "$FORCE" -eq 1 || ! -s "$REPEATS_BED" ]]; then
		ensure_dir "$FASTA_INDEX_DIR"
		log "REDItools: generating repeats (min-length=${REPEAT_MINLEN}) -> $REPEATS_BED"
		run_cmd docker_run reditools python -m reditools find-repeats \
			--min-length "$REPEAT_MINLEN" --output "$REPEATS_BED" "$REF_TO_USE"
	else
		log "REDItools: using existing repeats BED -> $REPEATS_BED"
	fi

	# -------------------- (A) Optional DNA branch --------------------
	local DNA_BAM="${OUT_DIR}/DNA_seq_bwa_aln/DNAseq.bam"
	local DO_DNA=0
	if [[ -s "$DNA_BAM" ]]; then
		DO_DNA=1
	elif has_dna_inputs; then
		log "REDItools(DNA): DNA reads found → running dna-align"
		step_dna_align || warn "DNA alignment failed or skipped"
		[[ -s "$DNA_BAM" ]] && DO_DNA=1
	fi

	local DNA_REDI_DIR="" DNA_MERGED=""
	if (( DO_DNA )); then
		DNA_REDI_DIR="${OUT_DIR}/DNA_reditools"; ensure_dir "$DNA_REDI_DIR"
		local dna_out="${DNA_REDI_DIR}/DNAseq_redi_out.tsv"

		if [[ ! -s "${DNA_BAM}.bai" ]]; then
			log "REDItools(DNA): indexing $DNA_BAM"
			run_cmd docker_run samtools samtools index "$DNA_BAM"
		else
			log "REDItools(DNA): BAM index exists, skipping"
		fi

		if [[ "$FORCE" -eq 1 || ! -s "$dna_out" ]]; then
			log "REDItools(DNA): analyzing $DNA_BAM"
			run_cmd docker_run reditools python -m reditools analyze "$DNA_BAM" \
				--reference "$REF_TO_USE" --threads "$THREADS" \
				--window "$RT_WINDOW" --min-read-depth 10 --variants all \
				--dna -o "$dna_out"
		else
			log "REDItools(DNA): $dna_out exists, skipping analyze"
		fi

		local arr=()
		for f in "$DNA_REDI_DIR"/*_redi_out.tsv; do
			[[ -s "$f" ]] && arr+=("$f")
		done
		if (( ${#arr[@]} > 0 )); then
			DNA_MERGED="${DNA_REDI_DIR}/DNA_redi_merged.tsv"
			if [[ "$FORCE" -eq 1 || ! -s "$DNA_MERGED" ]]; then
				log "REDItools(DNA): merging ${#arr[@]} file(s) -> $DNA_MERGED"
				{ head -n1 "${arr[0]}"; for f in "${arr[@]}"; do tail -n +2 "$f"; done; } > "$DNA_MERGED"
			else
				log "REDItools(DNA): merged file exists, skipping"
			fi
		fi
	else
		log "REDItools(DNA): no DNA → skipping DNA outputs"
	fi

	# -------------------- (B) RNA per-sample --------------------
	local sdir
	for sdir in "$(samples_dir)"/*; do
		[[ -d "$sdir" ]] || continue
		local sample; sample="$(basename "$sdir")"
		local bam; bam="$(bam_for_sample "$sample")"
		[[ -s "$bam" ]] || { warn "REDItools(RNA): no BAM for $sample"; continue; }

		local outd="${sdir}/reditools"; ensure_dir "$outd"

		# ============================================================
		# (1) BASIC OUTPUT WORKFLOW
		# ============================================================
		local rna_out_basic="${outd}/${sample}_redi_out_basic.tsv"
		local rna_idx_basic="${outd}/${sample}_redi_index_basic.tsv"

		# Run REDItools analyze on original BAM
		if [[ "$FORCE" -eq 1 || ! -s "$rna_out_basic" ]]; then
			log "REDItools(RNA): analyzing $sample on ORIGINAL BAM"
			run_cmd docker_run reditools python -m reditools analyze "$bam" \
				--reference "$REF_TO_USE" --strand "$SLIB" --threads "$THREADS" \
				--window "$RT_WINDOW" --min-read-depth 10 --variants all \
				-o "$rna_out_basic"
		else
			log "REDItools(RNA): basic analyze output exists for $sample, skipping"
		fi

		# Determine which file to cluster
		local rna_basic_to_cluster="$rna_out_basic"
		
		# ANNOTATE with DNA and filter DNA mismatches (if DNA available)
		if [[ -n "$DNA_MERGED" && -s "$DNA_MERGED" ]]; then
			local rna_basic_annotated_tmp="${outd}/${sample}_redi_out_basic_annotated.tmp.tsv"
			local rna_basic_annotated="${outd}/${sample}_redi_out_basic_annotated.tsv"
			
			if [[ "$FORCE" -eq 1 || ! -s "$rna_basic_annotated" ]]; then
				if [[ -s "$rna_out_basic" ]]; then
					log "REDItools(RNA): annotating basic output with DNA for $sample"
					run_cmd docker_run reditools bash -c \
						"python -m reditools annotate '$rna_out_basic' '$DNA_MERGED' > '$rna_basic_annotated_tmp'"
					
					if [[ -s "$rna_basic_annotated_tmp" ]]; then
						log "REDItools(RNA): filtering DNA mismatches from basic annotated output for $sample"
						filter_dna_mismatches "$rna_basic_annotated_tmp" "$rna_basic_annotated"
						rm -f "$rna_basic_annotated_tmp"
						
						if [[ -s "$rna_basic_annotated" ]]; then
							local kept_sites=$(tail -n +2 "$rna_basic_annotated" | wc -l)
							log "REDItools(RNA): DNA annotation completed for basic output ($sample) - $kept_sites sites"
						fi
					else
						warn "REDItools(RNA): DNA annotation produced empty output for basic ($sample)"
					fi
				fi
			else
				log "REDItools(RNA): annotated basic output exists for $sample, skipping"
			fi
			
			rna_basic_to_cluster="$rna_basic_annotated"
		fi

		# CLUSTER and FILTER the basic output (combined step)
		local rna_out_clustered="${outd}/${sample}_redi_out_basic_clustered.tsv"
		if [[ -n "$DNA_MERGED" && -s "$DNA_MERGED" ]]; then
			rna_out_clustered="${outd}/${sample}_redi_out_basic_annotated_clustered.tsv"
		fi
		
		if [[ "$FORCE" -eq 1 || ! -s "$rna_out_clustered" ]]; then
			if [[ -s "$rna_basic_to_cluster" ]]; then
				log "REDItools(RNA): clustering and filtering by AllSubs for $sample (basic)"
				
				local tmp_bed="${outd}/${sample}_basic.tmp.bed"
				local tmp_clustered="${outd}/${sample}_basic.clustered.bed"
				local tmp_clusterids="${outd}/${sample}_basic.clusterids.txt"
				local tmp_with_clusters="${outd}/${sample}_basic.tmp_clustered.tsv"
				
				# Create BED and run clustering
				tail -n +2 "$rna_basic_to_cluster" | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3, $4}' | \
					sort -k1,1 -k2,2n > "$tmp_bed"
				
				run_cmd docker_run bedtools bedtools cluster -d 100 -i "$tmp_bed" > "$tmp_clustered"
				
				cut -f6 "$tmp_clustered" > "$tmp_clusterids"
				
				# Add ClusterID to header, then paste data
				{
					head -n1 "$rna_basic_to_cluster" | sed 's/$/\tClusterID/'
					paste <(tail -n +2 "$rna_basic_to_cluster") "$tmp_clusterids"
				} > "$tmp_with_clusters"

				# Filter clusters by AllSubs type
				awk -F'\t' 'BEGIN {OFS="\t"}
				NR==1 {
					for(i=1; i<=NF; i++) {
						if($i == "ClusterID") cluster_col = i
						if($i == "AllSubs") allsubs_col = i
					}
					if(!cluster_col || !allsubs_col) {
						print "ERROR: ClusterID or AllSubs column not found" > "/dev/stderr"
						exit 1
					}
					print
					next
				}
				{
					rows[NR] = $0
					cluster = $cluster_col
					allsubs = $allsubs_col
					
					if(!(cluster in cluster_allsubs)) {
						cluster_allsubs[cluster] = allsubs
					} else if(cluster_allsubs[cluster] != allsubs) {
						cluster_allsubs[cluster] = "MULTIPLE"
					}
				}
				END {
					for(i=2; i<=NR; i++) {
						split(rows[i], fields, "\t")
						cluster = fields[cluster_col]
						if(cluster_allsubs[cluster] != "MULTIPLE") {
							print rows[i]
						}
					}
				}' "$tmp_with_clusters" > "$rna_out_clustered"
				
				rm -f "$tmp_bed" "$tmp_clustered" "$tmp_clusterids" "$tmp_with_clusters"
				
				if [[ -s "$rna_out_clustered" ]]; then
					local line_count=$(tail -n +2 "$rna_out_clustered" | wc -l)
					log "REDItools(RNA): clustering completed for $sample (basic) - $line_count sites"
				else
					warn "REDItools(RNA): clustering produced empty output for $sample (basic)"
				fi
			fi
		else
			log "REDItools(RNA): clustered output exists for $sample (basic), skipping"
		fi

		# INDEX basic output
		if [[ "$FORCE" -eq 1 || ! -s "$rna_idx_basic" ]]; then
			if [[ -s "$rna_out_basic" ]]; then
				log "REDItools(RNA): indexing basic output for $sample"
				run_cmd docker_run reditools python -m reditools index "$rna_out_basic" --strand "$SLIB" -o "$rna_idx_basic"
			fi
		fi

		# ============================================================
		# (2) FILTER REPEATS FROM BAM
		# ============================================================
		local bam_norep="${outd}/${sample}.norepeats.bam"
		local bam_norep_sorted="${outd}/${sample}.norepeats.sorted.bam"

		if [[ "$FORCE" -eq 1 || ! -s "$bam_norep_sorted" ]]; then
			log "REDItools(RNA): filtering out repeats from BAM for $sample"
			run_cmd docker_run bedtools bash -lc \
				"bedtools intersect -abam '${bam}' -b '${REPEATS_BED}' -v > '${bam_norep}'"
			run_cmd docker_run samtools samtools sort -o "${bam_norep_sorted}" "${bam_norep}"
			run_cmd docker_run samtools samtools index "${bam_norep_sorted}"
		fi

		# ============================================================
		# (3) FILTRATION WORKFLOW
		# ============================================================
		local rna_out_fil="${outd}/${sample}_redi_out_filtration.tsv"
		local rna_fil_idx="${outd}/${sample}_redi_index_filtration.tsv"
		local tmp_filtered="${outd}/${sample}_redi_out_filtration.tmp.tsv"

		# Run REDItools analyze with filters
		if [[ "$FORCE" -eq 1 || ! -s "$rna_out_fil" ]]; then
			log "REDItools(RNA): analyzing $sample with filtration"
			run_cmd docker_run reditools python -m reditools analyze "$bam_norep_sorted" \
				--reference "$REF_TO_USE" --threads "$THREADS" --window "$RT_WINDOW" \
				--min-read-depth 10 --strand "$SLIB" --variants all -o "$rna_out_fil" \
				--min-base-position 5 --max-base-position 5 --omopolymeric-span 5 \
				--strand-correction --strand-confidence-threshold 0.7 \
				--min-base-quality 30 --min-edits 2 --min-read-depth 10 \
				--max-editing-nucleotides 1

			log "REDItools(RNA): applying frequency threshold ${MIN_FREQ}"
			run_cmd docker_run reditools bash -c \
				"awk -v minfreq=${MIN_FREQ} 'BEGIN{FS=OFS=\"\t\"}
				NR==1{for(i=1;i<=NF;i++) if(\$i==\"Frequency\") col=i; print; next}
				\$col >= minfreq' \"$rna_out_fil\" > \"$tmp_filtered\" && \
				mv \"$tmp_filtered\" \"$rna_out_fil\""
		fi

		# Determine which file to cluster
		local rna_fil_to_cluster="$rna_out_fil"
		
		# ANNOTATE with DNA and filter DNA mismatches (if DNA available)
		if [[ -n "$DNA_MERGED" && -s "$DNA_MERGED" ]]; then
			local rna_fil_annotated_tmp="${outd}/${sample}_redi_out_filtration_annotated.tmp.tsv"
			local rna_fil_annotated="${outd}/${sample}_redi_out_filtration_annotated.tsv"
			
			if [[ "$FORCE" -eq 1 || ! -s "$rna_fil_annotated" ]]; then
				if [[ -s "$rna_out_fil" ]]; then
					log "REDItools(RNA): annotating filtration output with DNA for $sample"
					run_cmd docker_run reditools bash -c \
						"python -m reditools annotate '$rna_out_fil' '$DNA_MERGED' > '$rna_fil_annotated_tmp'"
					
					if [[ -s "$rna_fil_annotated_tmp" ]]; then
						log "REDItools(RNA): filtering DNA mismatches from filtration annotated output for $sample"
						filter_dna_mismatches "$rna_fil_annotated_tmp" "$rna_fil_annotated"
						rm -f "$rna_fil_annotated_tmp"
						
						if [[ -s "$rna_fil_annotated" ]]; then
							local kept_sites=$(tail -n +2 "$rna_fil_annotated" | wc -l)
							log "REDItools(RNA): DNA annotation completed for filtration output ($sample) - $kept_sites sites"
						fi
					else
						warn "REDItools(RNA): DNA annotation produced empty output for filtration ($sample)"
					fi
				fi
			else
				log "REDItools(RNA): annotated filtration output exists for $sample, skipping"
			fi
			
			rna_fil_to_cluster="$rna_fil_annotated"
		fi

		# CLUSTER and FILTER filtration output (combined step)
		local rna_out_fil_clustered="${outd}/${sample}_redi_out_filtration_clustered.tsv"
		if [[ -n "$DNA_MERGED" && -s "$DNA_MERGED" ]]; then
			rna_out_fil_clustered="${outd}/${sample}_redi_out_filtration_annotated_clustered.tsv"
		fi
		
		if [[ "$FORCE" -eq 1 || ! -s "$rna_out_fil_clustered" ]]; then
			if [[ -s "$rna_fil_to_cluster" ]]; then
				log "REDItools(RNA): clustering and filtering by AllSubs for $sample (filtration)"
				
				local tmp_bed="${outd}/${sample}_filtration.tmp.bed"
				local tmp_clustered="${outd}/${sample}_filtration.clustered.bed"
				local tmp_clusterids="${outd}/${sample}_filtration.clusterids.txt"
				local tmp_with_clusters="${outd}/${sample}_filtration.tmp_clustered.tsv"
				
				# Create BED and run clustering
				tail -n +2 "$rna_fil_to_cluster" | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3, $4}' | \
					sort -k1,1 -k2,2n > "$tmp_bed"
				
				run_cmd docker_run bedtools bedtools cluster -d 100 -i "$tmp_bed" > "$tmp_clustered"
				
				cut -f6 "$tmp_clustered" > "$tmp_clusterids"
				
				# Add ClusterID to header, then paste data
				{
					head -n1 "$rna_fil_to_cluster" | sed 's/$/\tClusterID/'
					paste <(tail -n +2 "$rna_fil_to_cluster") "$tmp_clusterids"
				} > "$tmp_with_clusters"

				# Filter clusters by AllSubs type
				awk -F'\t' 'BEGIN {OFS="\t"}
				NR==1 {
					for(i=1; i<=NF; i++) {
						if($i == "ClusterID") cluster_col = i
						if($i == "AllSubs") allsubs_col = i
					}
					if(!cluster_col || !allsubs_col) {
						print "ERROR: ClusterID or AllSubs column not found" > "/dev/stderr"
						exit 1
					}
					print
					next
				}
				{
					rows[NR] = $0
					cluster = $cluster_col
					allsubs = $allsubs_col
					
					if(!(cluster in cluster_allsubs)) {
						cluster_allsubs[cluster] = allsubs
					} else if(cluster_allsubs[cluster] != allsubs) {
						cluster_allsubs[cluster] = "MULTIPLE"
					}
				}
				END {
					for(i=2; i<=NR; i++) {
						split(rows[i], fields, "\t")
						cluster = fields[cluster_col]
						if(cluster_allsubs[cluster] != "MULTIPLE") {
							print rows[i]
						}
					}
				}' "$tmp_with_clusters" > "$rna_out_fil_clustered"
				
				rm -f "$tmp_bed" "$tmp_clustered" "$tmp_clusterids" "$tmp_with_clusters"
				
				if [[ -s "$rna_out_fil_clustered" ]]; then
					local line_count=$(tail -n +2 "$rna_out_fil_clustered" | wc -l)
					log "REDItools(RNA): clustering completed for $sample (filtration) - $line_count sites"
				else
					warn "REDItools(RNA): clustering produced empty output for $sample (filtration)"
				fi
			fi
		else
			log "REDItools(RNA): clustered output exists for $sample (filtration), skipping"
		fi

		# INDEX filtration output
		if [[ "$FORCE" -eq 1 || ! -s "$rna_fil_idx" ]]; then
			if [[ -s "$rna_out_fil" ]]; then
				log "REDItools(RNA): indexing filtration output for $sample"
				run_cmd docker_run reditools python -m reditools index "$rna_out_fil" --strand "$SLIB" -o "$rna_fil_idx"
			fi
		fi

		# ============================================================
		# (4) CDS FILTERING WORKFLOW
		# ============================================================
		if [[ "$ENABLE_CDS_FILTER" -eq 1 && -s "$CDS_BED" ]]; then
			local rna_cds="${outd}/${sample}_redi_out_filtration.cds.tsv"
			local rna_cds_idx="${outd}/${sample}_redi_index_filtration.cds.tsv"

			# Filter filtration output for CDS regions
			if [[ "$FORCE" -eq 1 || ! -s "$rna_cds" ]]; then
				if [[ -s "$rna_out_fil" ]]; then
					log "REDItools(CDS): filtering $sample for CDS regions"

					tail -n +2 "$rna_out_fil" | awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' > "${outd}/${sample}.tmp.bed"

					run_cmd docker_run bedtools bedtools intersect -a "${outd}/${sample}.tmp.bed" -b "$CDS_BED" -wa | \
						sort -k1,1 -k2,2n | uniq > "${outd}/${sample}.cds_sites.bed"

					awk -F'\t' '
					BEGIN {OFS="\t"}
					FNR==NR {
						if(FNR > 1) {
							sites[$1":"$3]=1
						}
						next
					}
					FNR==1 {
						print
						next
					}
					{
						key=$1":"$2
						if (key in sites) {
							print
						}
					}
					' "${outd}/${sample}.cds_sites.bed" "$rna_out_fil" > "$rna_cds"

					if [[ -s "$rna_cds" ]]; then
						local header_check=$(head -n1 "$rna_cds" | grep -c "Strand" || echo 0)
						if [[ "$header_check" -eq 0 ]]; then
							error "REDItools(CDS): CDS output missing Strand column header"
							rm -f "$rna_cds"
						else
							local orig_sites=$(tail -n +2 "$rna_out_fil" | wc -l)
							local cds_sites=$(tail -n +2 "$rna_cds" | wc -l)
							log "REDItools(CDS): kept $cds_sites / $orig_sites sites for $sample"
						fi
					fi

					rm -f "${outd}/${sample}.tmp.bed" "${outd}/${sample}.cds_sites.bed"
				fi
			fi

			# Determine which CDS file to cluster
			local rna_cds_to_cluster="$rna_cds"
			
			# ANNOTATE CDS with DNA and filter DNA mismatches (if DNA available)
			if [[ -n "$DNA_MERGED" && -s "$DNA_MERGED" ]]; then
				local rna_cds_annotated_tmp="${outd}/${sample}_redi_out_filtration.cds_annotated.tmp.tsv"
				local rna_cds_annotated="${outd}/${sample}_redi_out_filtration.cds_annotated.tsv"
				
				if [[ "$FORCE" -eq 1 || ! -s "$rna_cds_annotated" ]]; then
					if [[ -s "$rna_cds" ]]; then
						log "REDItools(CDS): annotating CDS output with DNA for $sample"
						run_cmd docker_run reditools bash -c \
							"python -m reditools annotate '$rna_cds' '$DNA_MERGED' > '$rna_cds_annotated_tmp'"
						
						if [[ -s "$rna_cds_annotated_tmp" ]]; then
							log "REDItools(CDS): filtering DNA mismatches from CDS annotated output for $sample"
							filter_dna_mismatches "$rna_cds_annotated_tmp" "$rna_cds_annotated"
							rm -f "$rna_cds_annotated_tmp"
							
							if [[ -s "$rna_cds_annotated" ]]; then
								local kept_sites=$(tail -n +2 "$rna_cds_annotated" | wc -l)
								log "REDItools(CDS): DNA annotation completed for CDS output ($sample) - $kept_sites sites"
							fi
						else
							warn "REDItools(CDS): DNA annotation produced empty output for CDS ($sample)"
						fi
					fi
				else
					log "REDItools(CDS): annotated CDS output exists for $sample, skipping"
				fi
				
				rna_cds_to_cluster="$rna_cds_annotated"
			fi

			# CLUSTER and FILTER CDS output (combined step)
			local rna_cds_clustered="${outd}/${sample}_redi_out_filtration.cds_clustered.tsv"
			if [[ -n "$DNA_MERGED" && -s "$DNA_MERGED" ]]; then
				rna_cds_clustered="${outd}/${sample}_redi_out_filtration.cds_annotated_clustered.tsv"
			fi
			
			if [[ "$FORCE" -eq 1 || ! -s "$rna_cds_clustered" ]]; then
				if [[ -s "$rna_cds_to_cluster" ]]; then
					log "REDItools(CDS): clustering and filtering by AllSubs for $sample (CDS)"
					
					local tmp_bed="${outd}/${sample}_cds.tmp.bed"
					local tmp_clustered="${outd}/${sample}_cds.clustered.bed"
					local tmp_clusterids="${outd}/${sample}_cds.clusterids.txt"
					local tmp_with_clusters="${outd}/${sample}_cds.tmp_clustered.tsv"
					
					# Create BED and run clustering
					tail -n +2 "$rna_cds_to_cluster" | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3, $4}' | \
						sort -k1,1 -k2,2n > "$tmp_bed"
					
					run_cmd docker_run bedtools bedtools cluster -d 100 -i "$tmp_bed" > "$tmp_clustered"
					
					cut -f6 "$tmp_clustered" > "$tmp_clusterids"
					
					# Add ClusterID to header, then paste data
					{
						head -n1 "$rna_cds_to_cluster" | sed 's/$/\tClusterID/'
						paste <(tail -n +2 "$rna_cds_to_cluster") "$tmp_clusterids"
					} > "$tmp_with_clusters"
					# Filter clusters by AllSubs type
					awk -F'\t' 'BEGIN {OFS="\t"}
					NR==1 {
						for(i=1; i<=NF; i++) {
							if($i == "ClusterID") cluster_col = i
							if($i == "AllSubs") allsubs_col = i
						}
						if(!cluster_col || !allsubs_col) {
							print "ERROR: ClusterID or AllSubs column not found" > "/dev/stderr"
							exit 1
						}
						print
						next
					}
					{
						rows[NR] = $0
						cluster = $cluster_col
						allsubs = $allsubs_col
						
						if(!(cluster in cluster_allsubs)) {
							cluster_allsubs[cluster] = allsubs
						} else if(cluster_allsubs[cluster] != allsubs) {
							cluster_allsubs[cluster] = "MULTIPLE"
						}
					}
					END {
						for(i=2; i<=NR; i++) {
							split(rows[i], fields, "\t")
							cluster = fields[cluster_col]
							if(cluster_allsubs[cluster] != "MULTIPLE") {
								print rows[i]
							}
						}
					}' "$tmp_with_clusters" > "$rna_cds_clustered"
					
					rm -f "$tmp_bed" "$tmp_clustered" "$tmp_clusterids" "$tmp_with_clusters"
					
					if [[ -s "$rna_cds_clustered" ]]; then
						local line_count=$(tail -n +2 "$rna_cds_clustered" | wc -l)
						log "REDItools(CDS): clustering completed for $sample (CDS) - $line_count sites"
					else
						warn "REDItools(CDS): clustering produced empty output for $sample (CDS)"
					fi
				fi
			else
				log "REDItools(CDS): clustered output exists for $sample (CDS), skipping"
			fi

			# INDEX CDS output
			if [[ "$FORCE" -eq 1 || ! -s "$rna_cds_idx" ]]; then
				if [[ -s "$rna_cds" ]]; then
					log "REDItools(CDS): indexing CDS output for $sample"
					run_cmd docker_run reditools python -m reditools index "$rna_cds" --strand "$SLIB" -o "$rna_cds_idx"
				fi
			fi
		fi

	done
	shopt -u nullglob
}

step_rei(){
	warn "REI: requires local installation -> skipping Docker containerization"
	return
}

jacusa_protocol_from_strand(){
	case "${1^^}" in
		2|FIRST|RF|RF-FIRSTSTRAND)	 echo "RF-FIRSTSTRAND" ;;
		1|SECOND|FR|FR-SECONDSTRAND) echo "FR-SECONDSTRAND" ;;
		0|UNSTRANDED|U|UNSTRAND)	 echo "UNSTRANDED" ;;
		*) warn "Unknown STRAND_LIB='$1' → default RF-FIRSTSTRAND"; echo "RF-FIRSTSTRAND" ;;
	esac
}

step_jacusa2(){
	ensure_jacusa_ref
	local proto; proto="$(jacusa_protocol_from_strand "$STRAND_LIB")"
	local opts=( -p "${THREADS}" -R "${JACUSA_REF}" -P "${proto}" )
	
	# Collect RNA BAMs
	local rna_bams=()
	local sdir
	for sdir in "$(samples_dir)"/*; do
		[[ -d "$sdir" ]] || continue
		local s; s="$(basename "$sdir")"
		local bam; bam="$(bam_for_sample "$s")"
		[[ -f "$bam" ]] && rna_bams+=("$bam")
	done
	
	(( ${#rna_bams[@]} > 0 )) || { warn "JACUSA2: no RNA BAMs"; return; }
	
	# Check if DNA BAM exists
	local DNA_BAM="${OUT_DIR}/DNA_seq_bwa_aln/DNAseq.bam"
	local pooled
	local has_dna=0
	
	if [[ -s "$DNA_BAM" ]]; then
		has_dna=1
		# DNA available: use call-2 (RNA vs DNA comparison)
		pooled="${OUT_DIR}/jacusa2_pooled_call2.txt"
		
		# Skip if output already exists
		if [[ -s "$pooled" && "$FORCE" -eq 0 && "$RESUME" -eq 1 ]]; then
			log "JACUSA2: call-2 output exists, skipping"
		else
			log "JACUSA2: DNA BAM found -> running call-2 (RNA vs DNA)"
			
			if [[ -n "$JACUSA_JAR" ]]; then
				run_cmd docker_run jacusa java -jar "$JACUSA_JAR" call-2 "${opts[@]}" \
					-r "$pooled" \
					"$(IFS=,; echo "${rna_bams[*]}")" \
					"$DNA_BAM"
			else
				run_cmd docker_run jacusa JACUSA2 call-2 "${opts[@]}" \
					-r "$pooled" \
					"$(IFS=,; echo "${rna_bams[*]}")" \
					"$DNA_BAM"
			fi
		fi
	else
		# No DNA: use call-1 (RNA-only variant calling)
		pooled="${OUT_DIR}/jacusa2_pooled_call1.txt"
		
		# Skip if output already exists
		if [[ -s "$pooled" && "$FORCE" -eq 0 && "$RESUME" -eq 1 ]]; then
			log "JACUSA2: call-1 output exists, skipping"
		else
			log "JACUSA2: No DNA BAM -> running call-1 (RNA-only)"
			
			if [[ -n "$JACUSA_JAR" ]]; then
				run_cmd docker_run jacusa java -jar "$JACUSA_JAR" call-1 "${opts[@]}" \
					-r "$pooled" \
					"$(IFS=,; echo "${rna_bams[*]}")"
			else
				run_cmd docker_run jacusa JACUSA2 call-1 "${opts[@]}" \
					-r "$pooled" \
					"$(IFS=,; echo "${rna_bams[*]}")"
			fi
		fi
	fi
	
	# Filter DNA SNPs if DNA data is available
	if [[ "$has_dna" -eq 1 && -s "$pooled" ]]; then
		local filtered="${OUT_DIR}/jacusa2_pooled_call2.filtered.txt"
		
		if [[ "$FORCE" -eq 1 || ! -s "$filtered" ]]; then
			log "JACUSA2: Filtering sites with DNA mismatches"
			
			awk -F'\t' 'BEGIN{OFS="\t"}
				NR==1 {
					# Find column indices (header parsing)
					for(i=1; i<=NF; i++) {
						lc=tolower($i); 
						gsub(/^[ \t]+|[ \t]+$/, "", lc);
						h[lc]=i;
					}
					# Find base columns for condition 2 (DNA)
					# JACUSA2 call-2 format: bases1 (RNA), bases2 (DNA)
					bases2_idx = h["bases2"];
					if(!bases2_idx) bases2_idx = h["base2"];
					if(!bases2_idx) bases2_idx = h["base22"];
					
					print; 
					next;
				}
				{
					# Check if DNA (condition 2) has any mismatch
					# bases2 format typically: A:10,C:0,G:0,T:0,N:0,DEL:0
					# A mismatch means non-reference base has count > 0
					
					if(!bases2_idx) {
						# Column not found, keep all sites
						print;
						next;
					}
					
					dna_bases = $bases2_idx;
					
					# Parse base counts (format: A:10,C:0,G:0,T:0,...)
					ref_base = $4;	# Reference base typically in column 4
					has_dna_mismatch = 0;
					
					# Split by comma
					split(dna_bases, base_counts, ",");
					for(i in base_counts) {
						# Split base:count
						split(base_counts[i], pair, ":");
						base = pair[1];
						count = pair[2] + 0;
						
						# Skip reference base, DEL, and N
						if(base == ref_base || base == "DEL" || base == "N") continue;
						
						# If non-reference base has count > 0, it'\''s a DNA mismatch
						if(count > 0) {
							has_dna_mismatch = 1;
							break;
						}
					}
					
					# Keep only sites WITHOUT DNA mismatch
					if(!has_dna_mismatch) {
						print;
					}
				}' "$pooled" > "$filtered"
			
			if [[ -s "$filtered" ]]; then
				local orig_count filtered_count
				orig_count=$(tail -n +2 "$pooled" | wc -l)
				filtered_count=$(tail -n +2 "$filtered" | wc -l)
				log "JACUSA2: Filtered $((orig_count - filtered_count)) DNA SNPs (kept ${filtered_count}/${orig_count} sites)"
			else
				warn "JACUSA2: Filtered output is empty"
			fi
		else
			log "JACUSA2: Filtered output exists, skipping"
		fi
	fi
	
	log "JACUSA2: Completed"
}
step_jacusa2_OLD(){
	local pooled="${OUT_DIR}/jacusa2_pooled_call1.txt"
	
	# Skip if pooled output already exists
	if [[ -s "$pooled" ]]; then
		log "JACUSA2: pooled output exists, skipping"
		return 0
	fi
	ensure_jacusa_ref
	local proto; proto="$(jacusa_protocol_from_strand "$STRAND_LIB")"
	local opts=( -p "${THREADS}" -R "${JACUSA_REF}" -P "${proto}" )
	local bams=()
	local sdir
	
	for sdir in "$(samples_dir)"/*; do
		[[ -d "$sdir" ]] || continue
		local s; s="$(basename "$sdir")"
		local bam; bam="$(bam_for_sample "$s")"
		[[ -f "$bam" ]] && bams+=("$bam")
	done
	
	(( ${#bams[@]} > 0 )) || { warn "JACUSA2: no BAMs"; return; }

	#local pooled="${OUT_DIR}/jacusa2_pooled_call1.txt"

	if [[ -n "$JACUSA_JAR" ]]; then
		# Use explicit JAR if the user provided one
		run_cmd docker_run jacusa java -jar "$JACUSA_JAR" call-1 "${opts[@]}" -r "$pooled" "$(IFS=,; echo "${bams[*]}")"
	else
		# Prefer the BioContainers launcher on PATH
		run_cmd docker_run jacusa JACUSA2 call-1 "${opts[@]}" -r "$pooled" "$(IFS=,; echo "${bams[*]}")"
	fi
}

step_hyperedit(){
	warn "HyperEdit: requires Python2 and specific setup -> skipping Docker containerization"
	return
}

ensure_sprint_image(){
	local tag="${SPRINT_IMAGE}"
	if ! docker image inspect "$tag" >/dev/null 2>&1; then
		log "SPRINT: image '$tag' not found → building from local Dockerfile"
		# Build context = directory containing this script (expects a file named 'Dockerfile')
		local script_dir; script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
		run_cmd docker build -t "$tag" "$script_dir"
	fi
	# Keep the container map in sync in case SPRINT_IMAGE was overridden
	CONTAINERS["sprint"]="$tag"
}

ensure_sprint_prepare(){
	ensure_fasta_index # creates ${FASTA_INDEX_DIR}/genomic.fa (+.fai)
	ensure_bwa_index # CRITICAL: ensure BWA index exists first
	ensure_sprint_image

	local prep_dir="${SPRINT_INDEX_DIR}"
	local mark="${prep_dir}/.prepare.done"
	local gtf_mark="${prep_dir}/.prepare_with_gtf"
	
	# Check if we need to re-run due to GTF change
	local needs_rerun=0
	if has_gtf && [[ ! -f "$gtf_mark" ]]; then
		log "SPRINT prepare: GTF now available, need to re-run with GTF annotation"
		needs_rerun=1
	elif ! has_gtf && [[ -f "$gtf_mark" ]]; then
		log "SPRINT prepare: GTF no longer available, need to re-run without GTF"
		needs_rerun=1
	fi
	
	# Check if prepare already completed with same GTF status
	if [[ "$FORCE" -eq 0 && -f "$mark" && "$needs_rerun" -eq 0 ]]; then
		log "SPRINT prepare: Already completed (mark found: $mark)"
		return 0
	fi

	ensure_dir "$prep_dir"

	# Symlink the reference FASTA into the SPRINT index directory
	ln -sf "${FASTA_INDEX_DIR}/genomic.fa" "${prep_dir}/genomic.fa"
	
	# Copy FASTA index (.fai) so SPRINT doesn't regenerate it
	if [[ -s "${FASTA_INDEX_DIR}/genomic.fa.fai" ]]; then
		cp "${FASTA_INDEX_DIR}/genomic.fa.fai" "${prep_dir}/genomic.fa.fai"
		log "SPRINT prepare: copied FASTA index → ${prep_dir}/genomic.fa.fai"
	else
		warn "SPRINT prepare: no .fai found at ${FASTA_INDEX_DIR}/genomic.fa.fai"
	fi

	# Run sprint prepare with GTF if available
	log "SPRINT prepare: Running full preparation (including BWA index and additional processing)..."
	
	# Build the sprint prepare command with optional GTF
	local sprint_prep_cmd="sprint prepare"
	if has_gtf && [[ -s "${FASTA_INDEX_DIR}/genomic.gtf" ]]; then
		sprint_prep_cmd+=" -t '${FASTA_INDEX_DIR}/genomic.gtf'"
		log "SPRINT prepare: running WITH GTF annotation"
	else
		log "SPRINT prepare: running WITHOUT GTF annotation"
	fi
	sprint_prep_cmd+=" '${prep_dir}/genomic.fa' '${SPRINT_BWA_PATH}'"
	
	run_cmd docker_run sprint bash -lc "
		set -euo pipefail
		cd '${prep_dir}'
		printf '[SPRINT prepare cwd] %s\n' \"\$(pwd)\"
		${sprint_prep_cmd}
	"

	# Mark completion and GTF status
	touch "$mark"
	if has_gtf; then
		touch "$gtf_mark"
	else
		rm -f "$gtf_mark"
	fi
}

ensure_sprint_prepare_OLD(){
	ensure_fasta_index		# creates ${FASTA_INDEX_DIR}/genomic.fa (+.fai)
	ensure_bwa_index		# CRITICAL: ensure BWA index exists first
	ensure_sprint_image

	local prep_dir="${SPRINT_INDEX_DIR}"
	local mark="${prep_dir}/.prepare.done"
	ensure_dir "$prep_dir"

	# Symlink the reference FASTA into the SPRINT index directory
	ln -sf "${FASTA_INDEX_DIR}/genomic.fa" "${prep_dir}/genomic.fa"
	
	# NEW: Copy FASTA index (.fai) so SPRINT doesn't regenerate it
	if [[ -s "${FASTA_INDEX_DIR}/genomic.fa.fai" ]]; then
		cp "${FASTA_INDEX_DIR}/genomic.fa.fai" "${prep_dir}/genomic.fa.fai"
		log "SPRINT prepare: copied FASTA index → ${prep_dir}/genomic.fa.fai"
	else
		warn "SPRINT prepare: no .fai found at ${FASTA_INDEX_DIR}/genomic.fa.fai"
	fi
	
	# CRITICAL FIX: Link ALL BWA index files to SPRINT prep directory
	log "SPRINT prepare: Linking BWA index files from ${BWA_INDEX_PREFIX}.*"
	for ext in amb ann bwt pac sa; do
		local src="${BWA_INDEX_PREFIX}.${ext}"
		local dst="${prep_dir}/genomic.fa.${ext}"
		if [[ -f "$src" ]]; then
			ln -sf "$src" "$dst"
			log "SPRINT prepare: Linked ${ext} index"
		else
			warn "SPRINT prepare: Missing BWA index file: $src"
		fi
	done

	# Skip actual sprint prepare if BWA index already linked and not forcing
	if [[ "$FORCE" -eq 0 && -s "${prep_dir}/genomic.fa.bwt" ]]; then
		log "SPRINT prepare: BWA index already present in ${prep_dir}, skipping prepare step"
		touch "$mark"
		return 0
	fi
	
	# Run sprint prepare only if BWA index doesn't exist
	log "SPRINT prepare: Building BWA index..."
	run_cmd docker_run sprint bash -lc "
		set -euo pipefail
		cd '${prep_dir}'
		printf '[SPRINT prepare cwd] %s\n' \"\$(pwd)\"
		sprint prepare '${prep_dir}/genomic.fa' '${SPRINT_BWA_PATH}'
	"
	
	touch "$mark"
}


SPRINT_CLEANUP_FASTQ="${SPRINT_CLEANUP_FASTQ:-1}"  # Set to 1 to auto-delete decompressed FASTQs
step_sprint(){
	# 1) REQUIRED: sprint prepare (once per project)
	ensure_sprint_prepare

	# === NEW: discover & clean REDItools repeats BED for optional SPRINT -rp ===
	local DEFAULT_RP="" CLEAN_BED="" REDI_RP=""
	shopt -s nullglob
	local rp_candidates=("${FASTA_INDEX_DIR}"/repeats*.bed)
	shopt -u nullglob
	if [[ ${#rp_candidates[@]} -gt 0 ]]; then
		IFS=$'\n' rp_candidates=($(printf "%s\n" "${rp_candidates[@]}" | sort))
		unset IFS
		REDI_RP="${rp_candidates[0]}"
		CLEAN_BED="${REDI_RP%.bed}.clean.bed"

		if [[ "$FORCE" -eq 1 || ! -s "$CLEAN_BED" || "$REDI_RP" -nt "$CLEAN_BED" ]]; then
			log "SPRINT: normalizing repeats to BED6 → $(basename "$CLEAN_BED")"
			run_cmd docker_run bedtools bash -lc "
				set -euo pipefail
				awk -v OFS='\t' '
					NF>=3 && \$1!~/^#/ {
						chr=\$1; s=\$2+0; e=\$3+0;
						if(e<=s) next;
						print chr, s, e
					}' '$REDI_RP' \
				| sort -k1,1 -k2,2n \
				| bedtools merge -i - \
				| awk -v OFS='\t' '{print \$1,\$2,\$3,\"repeat\",0,\".\"}' \
				> '${CLEAN_BED}.tmp'

				awk -v OFS='\t' 'NF==6 && \$3>\$2' '${CLEAN_BED}.tmp' > '${CLEAN_BED}'
				rm -f '${CLEAN_BED}.tmp'
			"
		else
			log "SPRINT: using cached cleaned repeats → $(basename "$CLEAN_BED")"
		fi

		[[ -s "$CLEAN_BED" ]] && DEFAULT_RP="$CLEAN_BED"
	else
		log "SPRINT: no REDItools repeat BED found under ${FASTA_INDEX_DIR} (optional)"
	fi
	# === END NEW BLOCK ===

	# 2) Per-sample: sprint main
	shopt -s nullglob
	local sdir
	for sdir in "$(samples_dir)"/*; do
		[[ -d "$sdir" ]] || continue
		local sample; sample="$(basename "$sdir")"
		local bam; bam="$(bam_for_sample "$sample")"
		[[ -s "$bam" ]] || { warn "SPRINT: no BAM for $sample"; continue; }

		local outd="${sdir}/sprint"; ensure_dir "$outd"
		local res="${outd}/SPRINT_identified_regular.res"

		if [[ "$RESUME" -eq 1 && -s "$res" && "$FORCE" -eq 0 ]]; then
			local line_count
			line_count=$(wc -l < "$res" 2>/dev/null || echo "0")
			if (( line_count > 10 )); then
				log "SPRINT: resume $sample ($line_count lines in result)"; continue
			else
				warn "SPRINT: Output exists but appears empty ($line_count lines), rerunning"
			fi
		fi

		local ss_opt=""
		case "${STRAND_LIB^^}" in
			2|FIRST|RF|RF-FIRSTSTRAND)
				ss_opt="-ss 1"
				log "SPRINT: Using fr-firststrand orientation (read1=sense)"
				;;
			1|SECOND|FR|FR-SECONDSTRAND)
				ss_opt="-ss 0"
				log "SPRINT: Using fr-secondstrand orientation (read1=antisense)"
				;;
			0|UNSTRANDED|U|UNSTRAND)
				ss_opt=""
				log "SPRINT: Unstranded library - omitting -ss flag (SPRINT default)"
				;;
			*)
				warn "SPRINT: Invalid STRAND_LIB='$STRAND_LIB', defaulting to -ss 1"
				ss_opt="-ss 1"
				;;
		esac

		# Get FASTQ files for this sample
		local fastp_dir="${sdir}/fastp"
		local r1_file="" r2_file=""

		shopt -s nullglob
		local r1_files=("$fastp_dir"/${FILTERED_PREFIX}*_1*.fastq* "$fastp_dir"/${FILTERED_PREFIX}*_R1*.fastq*)
		shopt -u nullglob

		if [[ ${#r1_files[@]} -gt 0 ]]; then
			r1_file="${r1_files[0]}"
			r2_file="$(mate_of "$r1_file")"
		else
			warn "SPRINT: No R1 FASTQ found for $sample, skipping"
			continue
		fi

		# === NEW: Handle compressed FASTQ files ===
		local r1_to_use="$r1_file"
		local r2_to_use="$r2_file"
		local need_cleanup=0

		# Check and decompress R1
		if [[ "$r1_file" =~ \.gz$ ]]; then
			log "SPRINT: R1 is compressed, decompressing to sprint directory"
			local r1_basename; r1_basename="$(basename "$r1_file" .gz)"
			r1_to_use="${outd}/${r1_basename}"
			
			# Decompress only if not already done or FORCE is set
			if [[ "$FORCE" -eq 1 || ! -s "$r1_to_use" || "$r1_file" -nt "$r1_to_use" ]]; then
				# Use pigz if available (faster), else gunzip
				if command -v pigz >/dev/null 2>&1; then
					run_cmd pigz -dc -p "${THREADS}" "$r1_file" > "$r1_to_use"
				else
					run_cmd gunzip -c "$r1_file" > "$r1_to_use"
				fi
			else
				log "SPRINT: Using cached decompressed R1"
			fi
			need_cleanup=1
		fi

		# Check and decompress R2 if exists
		if [[ -n "$r2_file" && -s "$r2_file" ]]; then
			if [[ "$r2_file" =~ \.gz$ ]]; then
				log "SPRINT: R2 is compressed, decompressing to sprint directory"
				local r2_basename; r2_basename="$(basename "$r2_file" .gz)"
				r2_to_use="${outd}/${r2_basename}"
				
				if [[ "$FORCE" -eq 1 || ! -s "$r2_to_use" || "$r2_file" -nt "$r2_to_use" ]]; then
					if command -v pigz >/dev/null 2>&1; then
						run_cmd pigz -dc -p "${THREADS}" "$r2_file" > "$r2_to_use"
					else
						run_cmd gunzip -c "$r2_file" > "$r2_to_use"
					fi
				else
					log "SPRINT: Using cached decompressed R2"
				fi
			fi
		fi
		# === END DECOMPRESSION BLOCK ===

		# Repeat file precedence: user-provided SPRINT_RP > cleaned REDItools BED > none
		local rp_opt=""
		if [[ -n "${SPRINT_RP:-}" && -s "$SPRINT_RP" ]]; then
			rp_opt="-rp $SPRINT_RP"
			log "SPRINT: using repeat file (from config): $SPRINT_RP"
		elif [[ -n "$DEFAULT_RP" && -s "$DEFAULT_RP" ]]; then
			rp_opt="-rp $DEFAULT_RP"
			log "SPRINT: using repeat file (REDItools-clean): $DEFAULT_RP"
		else
			log "SPRINT: no repeat file (optional)"
		fi

		# Build SPRINT main command with decompressed files
		local sprint_cmd="sprint main"
		sprint_cmd+=" -1 '$r1_to_use'"
		if [[ -n "$r2_to_use" && -s "$r2_to_use" ]]; then
			sprint_cmd+=" -2 '$r2_to_use'"
			log "SPRINT: processing $sample (PE mode)"
		else
			log "SPRINT: processing $sample (SE mode)"
		fi
		[[ -n "$rp_opt" ]] && sprint_cmd+=" $rp_opt"
		[[ -n "$ss_opt" ]] && sprint_cmd+=" $ss_opt"
		sprint_cmd+=" -cd ${SPRINT_CD:-200}"
		sprint_cmd+=" -p ${THREADS}"
		[[ -n "${SPRINT_ARGS:-}" ]] && sprint_cmd+=" ${SPRINT_ARGS}"
		sprint_cmd+=" '${SPRINT_INDEX_DIR}/genomic.fa' '$outd' '${SPRINT_BWA_PATH}' '${SPRINT_SAMTOOLS_PATH}'"

		log "SPRINT: running sprint main for ${sample}"
		run_cmd docker_run sprint bash -lc "$sprint_cmd"

		# === NEW: Optional cleanup of decompressed files ===
		# # Uncomment the following block to remove decompressed FASTQs after SPRINT completes
		# if [[ "$need_cleanup" -eq 1 && "${SPRINT_CLEANUP_FASTQ:-1}" -eq 1 ]]; then
			# log "SPRINT: Cleaning up decompressed FASTQ files"
			# [[ "$r1_to_use" != "$r1_file" && -f "$r1_to_use" ]] && rm -f "$r1_to_use"
			# [[ "$r2_to_use" != "$r2_file" && -f "$r2_to_use" ]] && rm -f "$r2_to_use"
		# fi
		# #	 === END CLEANUP BLOCK ===

		# Validate SPRINT output
		if [[ ! -s "$res" ]]; then
			error "SPRINT: Output file is empty or missing: $res"
			error "SPRINT: Check ${outd}/ for intermediate files and logs"
			continue
		fi
		local content_lines
		content_lines=$(tail -n +2 "$res" 2>/dev/null | wc -l || echo "0")
		if (( content_lines == 0 )); then
			warn "SPRINT: No editing sites detected for $sample (empty result)"
			warn "SPRINT: Possible causes: 1) Low coverage, 2) Wrong strand setting, 3) Over-masking by repeats, 4) No detectable A→G editing"
		else
			log "SPRINT: Detected $content_lines editing sites for $sample"
		fi
	done
	shopt -u nullglob
}
step_sprint_V1(){
	# 1) REQUIRED: sprint prepare (once per project)
	ensure_sprint_prepare

	# === NEW: discover & clean REDItools repeats BED for optional SPRINT -rp ===
	# We look for any repeats*.bed in FASTA_INDEX_DIR (e.g., repeats.min5.bed from REDItools find-repeats)
	local DEFAULT_RP="" CLEAN_BED="" REDI_RP=""
	shopt -s nullglob
	local rp_candidates=("${FASTA_INDEX_DIR}"/repeats*.bed)
	shopt -u nullglob
	if [[ ${#rp_candidates[@]} -gt 0 ]]; then
		# Prefer the most specific name if multiple exist (simple: take first sorted)
		IFS=$'\n' rp_candidates=($(printf "%s\n" "${rp_candidates[@]}" | sort))
		unset IFS
		REDI_RP="${rp_candidates[0]}"
		CLEAN_BED="${REDI_RP%.bed}.clean.bed"

		# Clean only if missing, forced, or source is newer
		if [[ "$FORCE" -eq 1 || ! -s "$CLEAN_BED" || "$REDI_RP" -nt "$CLEAN_BED" ]]; then
			log "SPRINT: normalizing repeats to BED6 → $(basename "$CLEAN_BED")"
			run_cmd docker_run bedtools bash -lc "
				set -euo pipefail
				# 1) Read REDItools repeats (>=3 cols), coerce to BED3, drop invalids, coerce to BED6
				awk -v OFS='\t' '
					NF>=3 && \$1!~/^#/ {
						chr=\$1; s=\$2+0; e=\$3+0;
						if(e<=s) next;				   # drop zero/neg intervals
						print chr, s, e
					}' '$REDI_RP' \
				| sort -k1,1 -k2,2n \
				| bedtools merge -i - \
				| awk -v OFS='\t' '{print \$1,\$2,\$3,\"repeat\",0,\".\"}' \
				> '${CLEAN_BED}.tmp'

				# 2) Final sanity: exactly 6 tab-separated fields, length>0, no blanks
				awk -v OFS='\t' 'NF==6 && \$3>\$2' '${CLEAN_BED}.tmp' > '${CLEAN_BED}'
				rm -f '${CLEAN_BED}.tmp'
			"
		else
			log "SPRINT: using cached cleaned repeats → $(basename "$CLEAN_BED")"
		fi


		[[ -s "$CLEAN_BED" ]] && DEFAULT_RP="$CLEAN_BED"
	else
		log "SPRINT: no REDItools repeat BED found under ${FASTA_INDEX_DIR} (optional)"
	fi
	# === END NEW BLOCK ===

	# 2) Per-sample: sprint main
	shopt -s nullglob
	local sdir
	for sdir in "$(samples_dir)"/*; do
		[[ -d "$sdir" ]] || continue
		local sample; sample="$(basename "$sdir")"
		local bam; bam="$(bam_for_sample "$sample")"
		[[ -s "$bam" ]] || { warn "SPRINT: no BAM for $sample"; continue; }

		local outd="${sdir}/sprint"; ensure_dir "$outd"
		local res="${outd}/SPRINT_identified_regular.res"

		if [[ "$RESUME" -eq 1 && -s "$res" && "$FORCE" -eq 0 ]]; then
			local line_count
			line_count=$(wc -l < "$res" 2>/dev/null || echo "0")
			if (( line_count > 10 )); then
				log "SPRINT: resume $sample ($line_count lines in result)"; continue
			else
				warn "SPRINT: Output exists but appears empty ($line_count lines), rerunning"
			fi
		fi

		local ss_opt=""
		case "${STRAND_LIB^^}" in
			2|FIRST|RF|RF-FIRSTSTRAND)
				ss_opt="-ss 1"
				log "SPRINT: Using fr-firststrand orientation (read1=sense)"
				;;
			1|SECOND|FR|FR-SECONDSTRAND)
				ss_opt="-ss 0"
				log "SPRINT: Using fr-secondstrand orientation (read1=antisense)"
				;;
			0|UNSTRANDED|U|UNSTRAND)
				ss_opt=""
				log "SPRINT: Unstranded library - omitting -ss flag (SPRINT default)"
				;;
			*)
				warn "SPRINT: Invalid STRAND_LIB='$STRAND_LIB', defaulting to -ss 1"
				ss_opt="-ss 1"
				;;
		esac

		# Get FASTQ files for this sample
		local fastp_dir="${sdir}/fastp"
		local r1_file="" r2_file=""

		shopt -s nullglob
		local r1_files=("$fastp_dir"/${FILTERED_PREFIX}*_1*.fastq* "$fastp_dir"/${FILTERED_PREFIX}*_R1*.fastq*)
		shopt -u nullglob

		if [[ ${#r1_files[@]} -gt 0 ]]; then
			r1_file="${r1_files[0]}"
			r2_file="$(mate_of "$r1_file")"
		else
			warn "SPRINT: No R1 FASTQ found for $sample, skipping"
			continue
		fi

		# Repeat file precedence: user-provided SPRINT_RP > cleaned REDItools BED > none
		local rp_opt=""
		if [[ -n "${SPRINT_RP:-}" && -s "$SPRINT_RP" ]]; then
			rp_opt="-rp $SPRINT_RP"
			log "SPRINT: using repeat file (from config): $SPRINT_RP"
		elif [[ -n "$DEFAULT_RP" && -s "$DEFAULT_RP" ]]; then
			rp_opt="-rp $DEFAULT_RP"
			log "SPRINT: using repeat file (REDItools-clean): $DEFAULT_RP"
		else
			log "SPRINT: no repeat file (optional)"
		fi

		# Build SPRINT main command
		local sprint_cmd="sprint main"
		sprint_cmd+=" -1 '$r1_file'"
		if [[ -n "$r2_file" && -s "$r2_file" ]]; then
			sprint_cmd+=" -2 '$r2_file'"
			log "SPRINT: processing $sample (PE mode)"
		else
			log "SPRINT: processing $sample (SE mode)"
		fi
		[[ -n "$rp_opt" ]] && sprint_cmd+=" $rp_opt"
		[[ -n "$ss_opt" ]] && sprint_cmd+=" $ss_opt"
		sprint_cmd+=" -cd ${SPRINT_CD:-200}"
		sprint_cmd+=" -p ${THREADS}"
		[[ -n "${SPRINT_ARGS:-}" ]] && sprint_cmd+=" ${SPRINT_ARGS}"
		sprint_cmd+=" '${SPRINT_INDEX_DIR}/genomic.fa' '$outd' '${SPRINT_BWA_PATH}' '${SPRINT_SAMTOOLS_PATH}'"

		log "SPRINT: running sprint main for ${sample}"
		run_cmd docker_run sprint bash -lc "$sprint_cmd"

		# Validate SPRINT output
		if [[ ! -s "$res" ]]; then
			error "SPRINT: Output file is empty or missing: $res"
			error "SPRINT: Check ${outd}/ for intermediate files and logs"
			continue
		fi
		local content_lines
		content_lines=$(tail -n +2 "$res" 2>/dev/null | wc -l || echo "0")
		if (( content_lines == 0 )); then
			warn "SPRINT: No editing sites detected for $sample (empty result)"
			warn "SPRINT: Possible causes: 1) Low coverage, 2) Wrong strand setting, 3) Over-masking by repeats, 4) No detectable A→G editing"
		else
			log "SPRINT: Detected $content_lines editing sites for $sample"
		fi
	done
	shopt -u nullglob
}
step_sprint_OLD(){
	# 1) REQUIRED: sprint prepare (once per project)
	ensure_sprint_prepare
	
	# 2) Per-sample: sprint main
	shopt -s nullglob
	local sdir
	for sdir in "$(samples_dir)"/*; do
		[[ -d "$sdir" ]] || continue
		local sample; sample="$(basename "$sdir")"
		local bam; bam="$(bam_for_sample "$sample")"
		[[ -s "$bam" ]] || { warn "SPRINT: no BAM for $sample"; continue; }
		
		local outd="${sdir}/sprint"; ensure_dir "$outd"
		local res="${outd}/SPRINT_identified_regular.res"
		
		# if [[ "$RESUME" -eq 1 && -s "$res" && "$FORCE" -eq 0 ]]; then
			# log "SPRINT: resume $sample"; continue
		# fi
		if [[ "$RESUME" -eq 1 && -s "$res" && "$FORCE" -eq 0 ]]; then
			local line_count
			line_count=$(wc -l < "$res" 2>/dev/null || echo "0")
			if (( line_count > 10 )); then	# More than just headers
				log "SPRINT: resume $sample ($line_count lines in result)"; continue
			else
				warn "SPRINT: Output exists but appears empty ($line_count lines), rerunning"
			fi
		fi
		
		local ss_opt=""
		case "${STRAND_LIB^^}" in
			2|FIRST|RF|RF-FIRSTSTRAND) 
				ss_opt="-ss 1"
				log "SPRINT: Using fr-firststrand orientation (read1=sense)"
				;;
			1|SECOND|FR|FR-SECONDSTRAND) 
				ss_opt="-ss 0"
				log "SPRINT: Using fr-secondstrand orientation (read1=antisense)"
				;;
			0|UNSTRANDED|U|UNSTRAND) 
				ss_opt=""  # Omit flag for unstranded data
				log "SPRINT: Unstranded library - omitting -ss flag (SPRINT will use default)"
				;;
			*) 
				warn "SPRINT: Invalid STRAND_LIB='$STRAND_LIB', defaulting to -ss 1"
				ss_opt="-ss 1"
				;;
		esac
		
		# Get FASTQ files for this sample
		local fastp_dir="${sdir}/fastp"
		local r1_file="" r2_file=""
		
		# Find R1 file
		shopt -s nullglob
		local r1_files=("$fastp_dir"/${FILTERED_PREFIX}*_1*.fastq* "$fastp_dir"/${FILTERED_PREFIX}*_R1*.fastq*)
		shopt -u nullglob
		
		if [[ ${#r1_files[@]} -gt 0 ]]; then
			r1_file="${r1_files[0]}"
			r2_file="$(mate_of "$r1_file")"
		else
			warn "SPRINT: No R1 FASTQ found for $sample, skipping"
			continue
		fi
		
		# ===== MODIFICATION 2: ONLY use global SPRINT_RP if provided =====
		local rp_opt=""
		if [[ -n "${SPRINT_RP:-}" && -s "$SPRINT_RP" ]]; then
			rp_opt="-rp $SPRINT_RP"
			log "SPRINT: using repeat file: $SPRINT_RP"
		else
			log "SPRINT: no repeat file (optional)"
		fi
		
		# Build SPRINT main command
		local sprint_cmd="sprint main"
		
		# Add R1 (required)
		sprint_cmd+=" -1 '$r1_file'"
		
		# Add R2 if paired-end
		if [[ -n "$r2_file" && -s "$r2_file" ]]; then
			sprint_cmd+=" -2 '$r2_file'"
			log "SPRINT: processing $sample (PE mode)"
		else
			log "SPRINT: processing $sample (SE mode)"
		fi
		
		# Add repeat file if available
		[[ -n "$rp_opt" ]] && sprint_cmd+=" $rp_opt"
		
		# Add strand-specific flag if defined
		[[ -n "$ss_opt" ]] && sprint_cmd+=" $ss_opt"
		
		# Add cluster distance parameter (use SPRINT default 200 if not specified)
		sprint_cmd+=" -cd ${SPRINT_CD:-200}"
		
		# Add threading
		sprint_cmd+=" -p ${THREADS}"
		
		# Add any additional user-specified arguments
		[[ -n "${SPRINT_ARGS:-}" ]] && sprint_cmd+=" ${SPRINT_ARGS}"
		
		# Add required positional arguments: reference_genome output_path bwa_path samtools_path
		sprint_cmd+=" '${SPRINT_INDEX_DIR}/genomic.fa' '$outd' '${SPRINT_BWA_PATH}' '${SPRINT_SAMTOOLS_PATH}'"
		
		log "SPRINT: running sprint main for ${sample}"
		run_cmd docker_run sprint bash -lc "$sprint_cmd"
		
		# Validate SPRINT output
		if [[ ! -s "$res" ]]; then
			error "SPRINT: Output file is empty or missing: $res"
			error "SPRINT: Check ${outd}/ for intermediate files and logs"
			continue
		fi
		local content_lines
		content_lines=$(tail -n +2 "$res" 2>/dev/null | wc -l || echo "0")
		if (( content_lines == 0 )); then
			warn "SPRINT: No editing sites detected for $sample (empty result)"
			warn "SPRINT: Possible causes: 1) Low coverage, 2) Wrong strand setting, 3) No A-to-G editing present"
		else
			log "SPRINT: Detected $content_lines editing sites for $sample"
		fi
	done
	shopt -u nullglob
}

step_lodei(){
	if [[ -z "${GROUPS_CSV:-}" ]]; then
		warn "LODEI: requires groups CSV (two groups). Skipping."
		return
	fi
	local sdir
	for sdir in "$(samples_dir)"/*; do
		[[ -d "$sdir" ]] || continue
		ensure_dir "${sdir}/lodei"
		[[ -s "${sdir}/lodei/README.txt" ]] || echo "LODEI integration hook – implement here." > "${sdir}/lodei/README.txt"
	done
}

# ----------------------------- Step Registry ---------------------------------
declare -A STEP_FUNCS=()
step_register(){ STEP_FUNCS["$1"]="$2"; }
step_exists(){ [[ -n "${STEP_FUNCS[$1]:-}" ]]; }
step_run(){ local s="$1"; step_exists "$s" || die "Unknown step: $s"; log "=== STEP: $s ==="; "${STEP_FUNCS[$s]}"; }

# ----------------------------- CLI & Config ----------------------------------
print_usage() {
	cat <<'EOF'

rnae – RNA Editing Pipeline (Docker-enabled, GTF-optional)

DESCRIPTION:
	A comprehensive RNA editing analysis pipeline using Docker containers for 
	reproducible results. Supports paired-end and single-end RNA-seq data with
	optional GTF annotation for enhanced functionality.

USAGE:
	./rnae_docker.sh [OPTIONS]

REQUIRED OPTIONS:
	--pipeline-home DIR		   Project home directory containing data and outputs
	--genome-fa FILE		   Reference genome FASTA file

CONFIGURATION:
	--config FILE			   Configuration file (recommended approach)
	-h, --help				   Show this help message and exit

INPUT/OUTPUT:
	--out-dir DIR			   Output directory (default: PIPELINE_HOME/rnae_out)
	--raw-dir DIR			   Directory with raw RNA FASTQ files
	--raw-prefix PREFIX		   Prefix pattern for raw FASTQ files (optional)
	--raw-suffix SUFFIX		   Suffix pattern for raw FASTQ files (default: .fastq)
	--layout [SE|PE|auto]	   Library layout (default: auto)

DNA DATA (Optional):
	--raw-dna-dir DIR		   Directory with raw DNA FASTQ files
	--raw-dna-pattern PATTERN  Pattern for DNA FASTQ files (default: *.fastq*)

REFERENCE GENOME:
	--genome-fa FILE		   Reference genome FASTA file (REQUIRED)
	--genome-gtf FILE		   Genome annotation GTF file (optional but recommended)
	--genome-gff FILE		   Genome annotation GFF file (alternative to GTF)

INDEX DIRECTORIES:
	--star-index-dir DIR	   Directory for STAR index files
	--bwa-index-dir DIR		   Directory for BWA index files
	--bwa-index-prefix PREFIX  Prefix for BWA index files

PIPELINE EXECUTION:
	--steps LIST			   Comma-separated pipeline steps to run:
								• fastp		  - Quality control and filtering
								• star-index  - Build STAR genome index
								• star-align  - RNA-seq alignment with STAR
								• dna-align	  - DNA alignment with BWA (if DNA provided)
								• reditools	  - RNA editing detection with REDItools
								• rei		  - RNA Editing Index calculation
								• jacusa2	  - Alternative RNA editing caller
								• sprint	  - SPRINT analysis integration
								• hyperedit	  - HyperEdit analysis
								• all		  - Run all steps (default)

	--threads N				   Number of CPU threads (default: 15)
	--resume				   Resume from previous run (default: enabled)
	--no-resume				   Start fresh, ignore previous outputs
	--force					   Overwrite existing outputs
	--dry-run				   Show commands without execution (testing mode)
	--verbose				   Enable verbose logging

READ LENGTH CONFIGURATION:
	--read-length N			   Read length for STAR index optimization (default: 75)
	--auto-detect-read-length  Auto-detect read length from first FASTQ file
	--no-auto-detect-read-length  Disable auto-detection

RNA EDITING ANALYSIS:
	--annotate				   Enable RNA editing site annotation (default: enabled)
	--no-annotate			   Disable annotation step
	--virus-pattern PATTERN	   Regex pattern to identify viral sequences
	--groups-csv FILE		   Sample grouping file for differential analysis

ANNOTATION SETTINGS:
	--annotator-bin FILE	   Path to annotation script (default: enhanced_gene_annotator.py)
	--annotator-threads N	   Threads for annotation process
	--annotator-batch-size N   Batch size for annotation processing
	--no-annotator-calc-rei	   Disable REI calculation during annotation

TOOL PATHS:
	--rei-bin FILE			   Path to RNAEditingIndex executable
	--jacusa-jar FILE		   Path to JACUSA2 jar file

FEATURES:
	✓ Docker containerized tools for reproducibility
	✓ GTF annotation optional (reduced functionality without GTF)
	✓ Automatic read length detection
	✓ Resume capability for interrupted runs
	✓ Memory-optimized STAR indexing
	✓ Support for both paired-end and single-end data
	✓ Comprehensive logging and error handling

EXAMPLES:

	# Basic run with config file:
	./rnae_docker.sh --config config.env --steps all --threads 16

	# Quality control and alignment only:
	./rnae_docker.sh --config config.env --steps fastp,star --threads 8

	# Custom read length for specific dataset:
	./rnae_docker.sh --config config.env --read-length 150 --steps star-index,star-align

	# Auto-detect read length and run RNA editing analysis:
	./rnae_docker.sh --config config.env --auto-detect-read-length --steps reditools

	# Minimal run without GTF (reduced functionality):
	./rnae_docker.sh --pipeline-home /data/project --genome-fa genome.fa --steps fastp,star,reditools

	# Dry run to test configuration:
	./rnae_docker.sh --config config.env --dry-run --verbose

DIRECTORY STRUCTURE:
	PIPELINE_HOME/
	├── Raw_Data/			   # Input RNA FASTQ files
	├── Genome/				   # Reference genome files
	│	└── indexes/		   # Generated index files
	├── samples/			   # Per-sample analysis results
	│	└── SAMPLE_NAME/
	│		├── fastp/		   # Quality control results
	│		├── star/		   # Alignment results
	│		└── reditools/	   # RNA editing results
	├── logs/				   # Pipeline execution logs
	└── config.env			   # Configuration file

REQUIREMENTS:
	• Docker installed and accessible
	• Input RNA FASTQ files (gzipped or uncompressed)
	• Reference genome FASTA file
	• Sufficient disk space for indexes and outputs

For more information, see the pipeline documentation or contact KranzlerLab.

EOF
}

init_layout(){
	mkdir -p "${OUT_DIR}/"{samples,logs,resources} \
			 "$FASTA_INDEX_DIR" "$STAR_INDEX_DIR" "$BWA_INDEX_DIR"
	LOG_FILE="${OUT_DIR}/logs/pipeline_$(date +%Y%m%d_%H%M%S).log"
	exec > >(tee -a "$LOG_FILE") 2>&1
}

expand_steps(){
	local in="${1,,}"
	if [[ "$in" == "all" ]]; then
		STEPS="fastp,star-index,star-align,detect-strand,reditools,rei,jacusa2,sprint,hyperedit"
		return
	fi
	
	local out=""
	IFS=',' read -r -a arr <<< "$in"
	local x
	for x in "${arr[@]}"; do
		x="$(echo "$x" | xargs)"
		case "$x" in
			all) out+="fastp,star-index,star-align,detect-strand,reditools,rei,jacusa2,sprint,hyperedit";;
			qc|fastp) out+=",fastp";;
			star) out+=",star-index,star-align";;
			dna) out+=",dna-align";;
			redi|reditools) out+=",reditools";;
			he|hyper|hyperedit) out+=",hyperedit";;
			post) out+="detect-strand,reditools,rei,jacusa2,sprint,hyperedit";;
			*) out+=",$x";;
		esac
	done
	
	local seen="," uniq=""
	IFS=',' read -r -a arr2 <<< "${out#,}"
	for x in "${arr2[@]}"; do
		[[ "$seen" == *",$x,"* ]] || { uniq+=",${x}"; seen+="$x,"; }
	done
	STEPS="${uniq#,}"
}

prescan_config(){
	local i=0; local args=("$@")
	while [[ $i -lt ${#args[@]} ]]; do
		local k="${args[$i]}"; local v="${args[$((i+1))]:-}"
		case "$k" in
			--config)
				[[ -f "$v" ]] || die "Config file not found: $v"
				# shellcheck source=/dev/null
				. "$v"
				i=$((i+2))
				;;
			*) i=$((i+1));;
		esac
	done
}

# ----------------------------- Main ------------------------------------------
main(){
	# Register steps
	step_register fastp		 step_fastp
	step_register star-index step_star_index
	step_register star-align step_star_align
	step_register detect-strand step_detect_strand
	step_register dna-align	 step_dna_align
	step_register reditools	 step_reditools
	step_register rei		 step_rei
	step_register jacusa2	 step_jacusa2
	step_register lodei		 step_lodei
	step_register sprint	 step_sprint
	step_register hyperedit	 step_hyperedit
	
	if [[ $# -eq 0 ]]; then print_usage; exit 0; fi
	
	# Pre-scan for --config and source it
	prescan_config "$@"
	
	# Parse CLI arguments
	ARGS=("$@")
	i=0
	while [[ $i -lt $# ]]; do
		k="${ARGS[$i]}"; v="${ARGS[$((i+1))]:-}"
		case "$k" in
			-h|--help) print_usage; exit 0;;
			--config) i=$((i+2));; # already handled
			--pipeline-home) PIPELINE_HOME="$v"; i=$((i+2));;
			--out-dir) OUT_DIR="$v"; i=$((i+2));;
			--raw-dir) RAW_DIR="$v"; i=$((i+2));;
			--raw-prefix) RAW_PREFIX="$v"; i=$((i+2));;
			--raw-suffix) RAW_SUFFIX="$v"; i=$((i+2));;
			--raw-dna-dir) RAW_DNA_DIR="$v"; i=$((i+2));;
			--raw-dna-pattern) RAW_DNA_PATTERN="$v"; i=$((i+2));;
			--layout) LAYOUT="$v"; i=$((i+2));;
			--genome-fa) GENOME_FA="$v"; i=$((i+2));;
			--genome-gtf) GENOME_GTF="$v"; i=$((i+2));;
			--genome-gff) GENOME_GFF="$v"; i=$((i+2));;
			--star-index-dir) STAR_INDEX_DIR="$v"; i=$((i+2));;
			--bwa-index-dir) BWA_INDEX_DIR="$v"; i=$((i+2));;
			--bwa-index-prefix) BWA_INDEX_PREFIX="$v"; i=$((i+2));;
			--threads) THREADS="$v"; i=$((i+2));;
			--steps) STEPS_REQ="$v"; i=$((i+2));;
			--resume) RESUME=1; i=$((i+1));;
			--no-resume) RESUME=0; i=$((i+1));;
			--force) FORCE=1; i=$((i+1));;
			--dry-run) DRY_RUN=1; i=$((i+1));;
			--verbose) VERBOSE=1; i=$((i+1));;
			--groups-csv) GROUPS_CSV="$v"; i=$((i+2));;
			--virus-pattern) VIRUS_PATTERN="$v"; i=$((i+2));;
			--read-length) READ_LENGTH="$v"; i=$((i+2));;
			--auto-detect-read-length) READ_LENGTH_AUTO_DETECT=1; i=$((i+1));;
			--no-auto-detect-read-length) READ_LENGTH_AUTO_DETECT=0; i=$((i+1));;
			--annotate) ANNOTATE_REDI=1; i=$((i+1));;
			--no-annotate) ANNOTATE_REDI=0; i=$((i+1));;
			--annotator-bin) ANNOTATOR_BIN="$v"; i=$((i+2));;
			--annotator-threads) ANNOTATOR_THREADS="$v"; i=$((i+2));;
			--annotator-batch-size) ANNOTATOR_BATCH_SIZE="$v"; i=$((i+2));;
			--no-annotator-calc-rei) ANNOTATOR_CALC_REI=0; i=$((i+1));;
			--strand-lib) STRAND_LIB="$v"; i=$((i+2));;
			--redi-index-strand) REDI_INDEX_STRAND="$v"; i=$((i+2));;
			--redi-window) REDI_WINDOW="$v"; i=$((i+2));;
			--he-runner) HE_RUNNER="$v"; i=$((i+2));;
			--no-he-analyze) HE_ANALYZE_ENABLE=0; i=$((i+1));;
			--rei-bin) REI_BIN="$v"; i=$((i+2));;
			--jacusa-jar) JACUSA_JAR="$v"; i=$((i+2));;
			--sprint-image) SPRINT_IMAGE="$v"; i=$((i+2));;
			--sprint-repeat) SPRINT_RP="$v"; i=$((i+2));;
			--sprint-cd) SPRINT_CD="$v"; i=$((i+2));;
			--sprint-args) SPRINT_ARGS="$v"; i=$((i+2));;
			--sprint-index-dir) SPRINT_INDEX_DIR="$v"; i=$((i+2));;
			#*) die "Unknown option: $k (use --help)";;
		esac
	done
	
	# Initialize
	init_layout
	expand_steps "${STEPS_REQ:-all}"

	# Auto detect read length if enabled
	#auto_detect_read_length

	# Preflight checks (GTF now optional)
	log "rnae-docker start | pipeline-home=${PIPELINE_HOME} | out=${OUT_DIR} | steps=${STEPS} | threads=${THREADS}"

	#####TEST######################################################### 
	# Load cached strand orientation if it exists and STRAND_LIB not set
	if [[ -z "${STRAND_LIB:-}" ]] || [[ "${STRAND_LIB}" == "auto" ]]; then
		local cache_file="${OUT_DIR}/strand_orientation.cache"
		if [[ -s "$cache_file" ]]; then
			source "$cache_file"
			log "Loaded cached STRAND_LIB=${STRAND_LIB} from previous detection"
		fi
	fi

	# Final fallback
	STRAND_LIB="${STRAND_LIB:-2}"
	log "Using STRAND_LIB=${STRAND_LIB}"
	########################################################################
	
	
	# GTF is now optional - just warn if missing
	if ! has_gtf; then
		warn "No GTF file specified - pipeline will run with reduced functionality"
		warn "STAR: gene counting disabled, splice junction detection reduced"
		warn "REDItools: annotation step will be skipped"
		warn "REDItools: annotation step will be skipped"
	else
		log "GTF file found: ${GENOME_GTF} - running with full functionality"
	fi

	# Check Docker availability
	command -v docker >/dev/null 2>&1 || die "Docker is required but not found in PATH"
	
	# Check RNA inputs
	compgen -G "${RAW_DIR}/${RAW_PREFIX}*${RAW_SUFFIX}" >/dev/null || die "No RNA files found matching: ${RAW_DIR}/${RAW_PREFIX}*${RAW_SUFFIX}"
	
	# Run pipeline steps
	local s
	IFS=',' read -r -a STEPS_ARR <<< "$STEPS"
	for s in "${STEPS_ARR[@]}"; do
		step_run "$s"
	done
	
	log "rnae-docker completed successfully."
}

main "$@"