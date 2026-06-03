#!/bin/zsh

set -u
set -o pipefail
setopt NO_NOMATCH

# =========================
# User settings
# =========================
K=40   # 并行 job 数
MAX_CONCURRENT=40   # 同时最多跑多少个，防止机器打满
OUTPUT_W_CSV=1      # 设为 1 或使用 --w-csv：只输出 corrected W CSV，不输出 ROOT

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
EXE="$SCRIPT_DIR/../build/AnalysisDVCS"
OUTPUT_BASE="$SCRIPT_DIR/../build/dvcs_parallel_output"
typeset -a INPUT_DIRS
INPUT_DIRS=(
  /cache/clas12/rg-k/production/recon/fall2018/torus+1/7546MeV/pass2/v0/dst/recon/*/
)
MAX_TOTAL_FILES=10000000000

# AnalysisDVCS 第3个参数
ARG3=1

# root 输出文件
typeset -a OUTPUT_FILES
OUTPUT_FILES=(
  dfSelected_afterFid_afterCorr.root
  dfSelected_afterFid.root
  dfSelected.root
)

# csv 输出文件（如果某些 job 没有产生，也没关系）
typeset -a CSV_FILES
CSV_FILES=(
  events_per_run_afterFid.csv
)

SUCCESS_GREP='Finished processing all events'

usage() {
  cat <<EOF
Usage: $0 [--w-csv]

Options:
  --w-csv, --output-w-csv   Run AnalysisDVCS in corrected-W CSV-only mode.
  --help, -h               Show this help message.
EOF
}

while (( $# > 0 )); do
  case "$1" in
    --w-csv|--output-w-csv)
      OUTPUT_W_CSV=1
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown option: $1"
      usage
      exit 1
      ;;
  esac
  shift
done

if (( OUTPUT_W_CSV )); then
  OUTPUT_FILES=()
  CSV_FILES=(electron_w_afterCorr.csv)
  SUCCESS_GREP='Wrote final corrected electron W CSV'
fi

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
WORKDIR="$OUTPUT_BASE/jobs_${TIMESTAMP}"
LISTDIR="$WORKDIR/lists"
MASTER_LOG="$WORKDIR/master.log"
RUN_CONFIG_LOG="$WORKDIR/run_config.txt"
SUMMARY_TXT="$WORKDIR/summary.txt"
SUMMARY_CSV="$WORKDIR/job_summary.csv"

# =========================
# Helpers
# =========================
log_msg() {
  local msg="$1"
  echo "$msg" | tee -a "$MASTER_LOG"
}

write_run_config() {
  cat > "$RUN_CONFIG_LOG" <<EOF
timestamp=$TIMESTAMP
exe=$EXE
input_dirs=${(j:,:)INPUT_DIRS}
output_base=$OUTPUT_BASE
workdir=$WORKDIR
K=$K
max_concurrent=$MAX_CONCURRENT
arg3=$ARG3
output_w_csv=$OUTPUT_W_CSV
success_grep=$SUCCESS_GREP

root_outputs=${(j:,:)OUTPUT_FILES}
csv_outputs=${(j:,:)CSV_FILES}
EOF
}

merge_csv_keep_header_once() {
  local outfile="$1"
  shift
  local -a infiles=("$@")

  (( ${#infiles[@]} > 0 )) || return 1

  : > "$outfile"

  local first=1
  local f
  for f in "${infiles[@]}"; do
    [[ -s "$f" ]] || continue
    if (( first )); then
      cat "$f" >> "$outfile"
      first=0
    else
      tail -n +2 "$f" >> "$outfile"
    fi
  done

  [[ -s "$outfile" ]]
}

sum_run_events_csv() {
  local infile="$1"
  local outfile="$2"

  [[ -s "$infile" ]] || return 1

  awk -F',' '
    NR==1 { next }
    NF>=2 {
      gsub(/^[ \t]+|[ \t]+$/, "", $1)
      gsub(/^[ \t]+|[ \t]+$/, "", $2)
      if ($1 != "" && $2 != "") sum[$1] += $2
    }
    END {
      print "run,events_afterFid"
      for (r in sum) print r "," sum[r]
    }
  ' "$infile" | sort -t',' -k1,1n > "$outfile"
}

# =========================
# Checks
# =========================
if [[ ! -x "$EXE" ]]; then
  echo "[ERROR] Executable not found or not executable: $EXE"
  exit 1
fi

for input_dir in "${INPUT_DIRS[@]}"; do
  if [[ ! -d "$input_dir" ]]; then
    echo "[ERROR] Input directory not found: $input_dir"
    exit 1
  fi
done

if (( ! OUTPUT_W_CSV )) && ! command -v hadd >/dev/null 2>&1; then
  echo "[ERROR] hadd not found. Please source ROOT first."
  exit 1
fi

mkdir -p "$OUTPUT_BASE"
OUTPUT_BASE="$(cd "$OUTPUT_BASE" && pwd)"

WORKDIR="$OUTPUT_BASE/jobs_${TIMESTAMP}"
LISTDIR="$WORKDIR/lists"
MASTER_LOG="$WORKDIR/master.log"
RUN_CONFIG_LOG="$WORKDIR/run_config.txt"
SUMMARY_TXT="$WORKDIR/summary.txt"
SUMMARY_CSV="$WORKDIR/job_summary.csv"

mkdir -p "$WORKDIR" "$LISTDIR"
: > "$MASTER_LOG"

write_run_config

log_msg "[INFO] EXE            = $EXE"
log_msg "[INFO] INPUT_DIRS     = ${(j:,:)INPUT_DIRS}"
log_msg "[INFO] OUTPUT_BASE    = $OUTPUT_BASE"
log_msg "[INFO] WORKDIR        = $WORKDIR"
log_msg "[INFO] K              = $K"
log_msg "[INFO] MAX_CONCURRENT = $MAX_CONCURRENT"
log_msg "[INFO] ARG3           = $ARG3"
log_msg "[INFO] OUTPUT_W_CSV   = $OUTPUT_W_CSV"

# =========================
# Collect input files
# =========================
ALLLIST="$LISTDIR/all_files.txt"
for input_dir in "${INPUT_DIRS[@]}"; do
  find "$input_dir" -type f -name "*.hipo"
done | sort | head -n $MAX_TOTAL_FILES > "$ALLLIST"

NF=$(wc -l < "$ALLLIST")
if (( NF == 0 )); then
  log_msg "[ERROR] No .hipo files found in: ${(j:,:)INPUT_DIRS}"
  exit 1
fi

log_msg "[INFO] Found $NF hipo files."

if (( K > NF )); then
  K=$NF
  log_msg "[INFO] Reduce K to $K because total file count is only $NF."
fi

# =========================
# Split files round-robin
# =========================
for (( i=0; i<K; i++ )); do
  kid=$(printf "%03d" $i)
  : > "$LISTDIR/files_${kid}.txt"
done

j=0
while IFS= read -r f; do
  kid=$(printf "%03d" $j)
  echo "$f" >> "$LISTDIR/files_${kid}.txt"
  ((j++))
  (( j >= K )) && j=0
done < "$ALLLIST"

# =========================
# Prepare job directories
# =========================
typeset -a JOBS_TO_RUN
JOBS_TO_RUN=()

for (( i=0; i<K; i++ )); do
  kid=$(printf "%03d" $i)
  LIST="$LISTDIR/files_${kid}.txt"
  [[ -s "$LIST" ]] || continue

  JOBDIR="$WORKDIR/job_${kid}"
  INDIR="$JOBDIR/input"
  OUTDIR="$JOBDIR/out"
  LOG="$JOBDIR/job_${kid}.log"
  FILEMANIFEST="$JOBDIR/input_files.txt"

  mkdir -p "$INDIR" "$OUTDIR"
  rm -f "$INDIR"/*.hipo 2>/dev/null

  cp "$LIST" "$FILEMANIFEST"

  while IFS= read -r ff; do
    ln -sf "$ff" "$INDIR/"
  done < "$LIST"

  JOBS_TO_RUN+=("$kid")
done

log_msg "[INFO] Will run ${#JOBS_TO_RUN[@]} jobs."

# =========================
# Run jobs with concurrency limit
# =========================
typeset -A PID_TO_JOB
running=0

for kid in "${JOBS_TO_RUN[@]}"; do
  while (( running >= MAX_CONCURRENT )); do
    for pid in ${(k)PID_TO_JOB}; do
      if ! kill -0 "$pid" 2>/dev/null; then
        wait "$pid" 2>/dev/null
        log_msg "[INFO] Job ${PID_TO_JOB[$pid]} finished."
        unset PID_TO_JOB[$pid]
        ((running--))
      fi
    done
    sleep 1
  done

  JOBDIR="$WORKDIR/job_${kid}"
  INDIR="$JOBDIR/input"
  OUTDIR="$JOBDIR/out"
  LOG="$JOBDIR/job_${kid}.log"
  LIST="$LISTDIR/files_${kid}.txt"

  NFILES=$(wc -l < "$LIST")
  if (( NFILES <= 0 )); then
    log_msg "[WARN] Job $kid has no input files, skip."
    continue
  fi

  (
    cd "$OUTDIR" || exit 1
    echo "[INFO] Running job $kid"
    echo "[INFO] Working directory: $OUTDIR"
    echo "[INFO] Input directory: $INDIR"
    echo "[INFO] Number of files passed to AnalysisDVCS: $NFILES"
    if (( OUTPUT_W_CSV )); then
      echo "[INFO] CSV-only W output mode enabled."
      echo "[INFO] Full command: DISANA_OUTPUT_W_CSV=1 $EXE $INDIR $NFILES $ARG3"
    else
      echo "[INFO] Full command: $EXE $INDIR $NFILES $ARG3"
    fi
    echo "[INFO] Input file list:"
    cat "$JOBDIR/input_files.txt"
    if (( OUTPUT_W_CSV )); then
      DISANA_OUTPUT_W_CSV=1 "$EXE" "$INDIR" "$NFILES" "$ARG3"
    else
      "$EXE" "$INDIR" "$NFILES" "$ARG3"
    fi
  ) &> "$LOG" &

  pid=$!
  PID_TO_JOB[$pid]="$kid"
  ((running++))
  log_msg "[INFO] Started job $kid with PID $pid"
done

while (( running > 0 )); do
  for pid in ${(k)PID_TO_JOB}; do
    if ! kill -0 "$pid" 2>/dev/null; then
      wait "$pid" 2>/dev/null
      log_msg "[INFO] Job ${PID_TO_JOB[$pid]} finished."
      unset PID_TO_JOB[$pid]
      ((running--))
    fi
  done
  sleep 1
done

log_msg "[INFO] All jobs finished."

# =========================
# Check outputs
# =========================
typeset -a OKJ
typeset -a BADJ
OKJ=()
BADJ=()

echo "job_id,job_status,n_input_files,root_outputs_found,csv_outputs_found,log_file" > "$SUMMARY_CSV"

for kid in "${JOBS_TO_RUN[@]}"; do
  JOBDIR="$WORKDIR/job_${kid}"
  OUTDIR="$JOBDIR/out"
  LOG="$JOBDIR/job_${kid}.log"
  LIST="$LISTDIR/files_${kid}.txt"
  NFILES=$(wc -l < "$LIST")

  ok_this_job=0
  root_found=0
  csv_found=0

  if [[ -f "$LOG" ]]; then
    for f in "${OUTPUT_FILES[@]}"; do
      if [[ -s "$OUTDIR/$f" ]]; then
        ok_this_job=1
        ((root_found++))
      fi
    done

    for f in "${CSV_FILES[@]}"; do
      [[ -s "$OUTDIR/$f" ]] && ((csv_found++))
    done
  fi

  if (( OUTPUT_W_CSV )); then
    (( csv_found > 0 )) && ok_this_job=1
  else
    (( root_found > 0 )) && ok_this_job=1
  fi

  if (( ok_this_job )); then
    OKJ+=("$kid")
    job_status="OK"
  else
    BADJ+=("$kid")
    job_status="FAILED"
  fi

  echo "$kid,$job_status,$NFILES,$root_found,$csv_found,$LOG" >> "$SUMMARY_CSV"
done

{
  echo "Successful jobs: ${#OKJ[@]}"
  echo "Failed jobs    : ${#BADJ[@]}"
  if (( ${#BADJ[@]} > 0 )); then
    echo "Failed job IDs : ${BADJ[*]}"
  fi
} > "$SUMMARY_TXT"

log_msg "[INFO] Successful jobs: ${#OKJ[@]}"
log_msg "[INFO] Failed jobs    : ${#BADJ[@]}"
if (( ${#BADJ[@]} > 0 )); then
  log_msg "[WARN] Failed job IDs: ${BADJ[*]}"
fi

# =========================
# Merge ROOT outputs
# =========================
MERGED_DIR="$OUTPUT_BASE/merged_${TIMESTAMP}"
mkdir -p "$MERGED_DIR"

for f in "${OUTPUT_FILES[@]}"; do
  typeset -a merge_inputs
  merge_inputs=()

  for kid in "${OKJ[@]}"; do
    OUTDIR="$WORKDIR/job_${kid}/out"
    [[ -s "$OUTDIR/$f" ]] && merge_inputs+=("$OUTDIR/$f")
  done

  if (( ${#merge_inputs[@]} > 0 )); then
    log_msg "[INFO] Merging ROOT $f from ${#merge_inputs[@]} jobs ..."
    hadd -f "$MERGED_DIR/$f" "${merge_inputs[@]}" >> "$MASTER_LOG" 2>&1
  else
    log_msg "[INFO] No job produced ROOT $f, skip."
  fi
done

# =========================
# Merge CSV outputs
# =========================
for f in "${CSV_FILES[@]}"; do
  typeset -a csv_inputs
  csv_inputs=()

  for kid in "${OKJ[@]}"; do
    OUTDIR="$WORKDIR/job_${kid}/out"
    [[ -s "$OUTDIR/$f" ]] && csv_inputs+=("$OUTDIR/$f")
  done

  if (( ${#csv_inputs[@]} > 0 )); then
    raw_out="$MERGED_DIR/${f:r}_merged_raw.csv"
    log_msg "[INFO] Merging CSV $f from ${#csv_inputs[@]} jobs -> $raw_out"
    merge_csv_keep_header_once "$raw_out" "${csv_inputs[@]}"

    # 如果正好是 run,events_afterFid 这种格式，再额外按 run 求和
    if [[ "$f" == "events_afterFid.csv" ]]; then
      sum_out="$MERGED_DIR/${f:r}_merged_summed.csv"
      sum_run_events_csv "$raw_out" "$sum_out"
      log_msg "[INFO] Summed CSV by run -> $sum_out"
    fi
  else
    log_msg "[INFO] No job produced CSV $f, skip."
  fi
done

log_msg "[INFO] Done."
log_msg "[INFO] Master log      : $MASTER_LOG"
log_msg "[INFO] Run config      : $RUN_CONFIG_LOG"
log_msg "[INFO] Summary text    : $SUMMARY_TXT"
log_msg "[INFO] Summary csv     : $SUMMARY_CSV"
log_msg "[INFO] Merged outputs  : $MERGED_DIR"
