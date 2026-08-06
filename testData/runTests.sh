#!/bin/bash
#
# Regression test for deML: demultiplexes testData/ in both BAM and fastq mode
# and compares md5sums of the output against testData/expected.md5.
#
# The md5s are taken over *normalized* output, never over the raw files: BAM and
# fastq.gz are compressed, and the exact bytes depend on the zlib/htslib version
# and on the number of BGZF threads. Normalizing means
#
#   BAM      -> SAM text (via htslib/htsfile), minus the @PG lines, which carry
#               the command line and program version
#   fastq.gz -> the uncompressed fastq
#   summary  -> already plain text, used as is
#
# Usage:
#   testData/runTests.sh            # run the tests, compare against expected.md5
#   testData/runTests.sh --update   # regenerate testData/expected.md5
#
# Run it from the top of the repo, or via "make test" / "make test-update".

set -u

BASEDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )
DEML="${BASEDIR}/src/deML"
HTSFILE="${BASEDIR}/htslib/htsfile"
DATA="${BASEDIR}/testData"
EXPECTED="${DATA}/expected.md5"
OUT="${DATA}/testout"

UPDATE=0
if [ $# -gt 0 ]; then
    case "$1" in
	--update) UPDATE=1 ;;
	*) echo "Usage: $0 [--update]" >&2; exit 2 ;;
    esac
fi

for prog in "$DEML" "$HTSFILE"; do
    if [ ! -x "$prog" ]; then
	echo "TEST ERROR: $prog is missing, run \"make\" first" >&2
	exit 2
    fi
done

# md5sum on Linux, md5 on macOS; both reduced to the bare hash.
if command -v md5sum >/dev/null 2>&1; then
    md5of(){ md5sum | cut -d' ' -f1; }
elif command -v md5 >/dev/null 2>&1; then
    md5of(){ md5 -q; }
else
    echo "TEST ERROR: neither md5sum nor md5 was found" >&2
    exit 2
fi

rm -rf "$OUT"
mkdir -p "$OUT/fq" "$OUT/norm/fq"

echo "[1/4] BAM input, single threaded"
"$DEML" -i "${DATA}/index.txt" \
	-s "${OUT}/bam.summary.txt" \
	-o "${OUT}/demultiplexed.bam" \
	"${DATA}/todemultiplex.bam" > "${OUT}/bam.stdout" 2>&1
if [ $? -ne 0 ]; then
    echo "TEST FAILED: deML exited non-zero on BAM input" >&2
    cat "${OUT}/bam.stdout" >&2
    exit 1
fi

echo "[2/4] BAM input, 4 threads (output must stay in input order)"
"$DEML" -i "${DATA}/index.txt" -t 4 \
	-s "${OUT}/bam.t4.summary.txt" \
	-o "${OUT}/demultiplexed.t4.bam" \
	"${DATA}/todemultiplex.bam" > "${OUT}/bam.t4.stdout" 2>&1
if [ $? -ne 0 ]; then
    echo "TEST FAILED: deML exited non-zero on BAM input with -t 4" >&2
    cat "${OUT}/bam.t4.stdout" >&2
    exit 1
fi

echo "[3/4] fastq input"
"$DEML" -i "${DATA}/index.txt" \
	-s "${OUT}/fq.summary.txt" \
	-f  "${DATA}/todemultiplex.fq1.gz" \
	-r  "${DATA}/todemultiplex.fq2.gz" \
	-if1 "${DATA}/todemultiplex.i1.gz" \
	-if2 "${DATA}/todemultiplex.i2.gz" \
	-o  "${OUT}/fq/demultiplexed" > "${OUT}/fq.stdout" 2>&1
if [ $? -ne 0 ]; then
    echo "TEST FAILED: deML exited non-zero on fastq input" >&2
    cat "${OUT}/fq.stdout" >&2
    exit 1
fi

echo "[4/4] normalizing output and checking md5sums"

# BAM -> SAM text without the @PG lines (command line and version are not part
# of the result being tested).
"$HTSFILE" --view "${OUT}/demultiplexed.bam"    | grep -v '^@PG' > "${OUT}/norm/demultiplexed.sam"
"$HTSFILE" --view "${OUT}/demultiplexed.t4.bam" | grep -v '^@PG' > "${OUT}/norm/demultiplexed.t4.sam"

cp "${OUT}/bam.summary.txt"    "${OUT}/norm/bam.summary.txt"
cp "${OUT}/bam.t4.summary.txt" "${OUT}/norm/bam.t4.summary.txt"
cp "${OUT}/fq.summary.txt"     "${OUT}/norm/fq.summary.txt"

for gz in "${OUT}"/fq/*.fq.gz; do
    base=$( basename "$gz" .gz )
    gunzip -c "$gz" > "${OUT}/norm/fq/${base}"
done

# -t 4 must reproduce the single threaded result exactly; compare here rather
# than through expected.md5 so a mismatch names the actual problem.
if ! cmp -s "${OUT}/norm/demultiplexed.sam" "${OUT}/norm/demultiplexed.t4.sam"; then
    echo "TEST FAILED: -t 4 output differs from -t 1 output" >&2
    exit 1
fi
if ! cmp -s "${OUT}/norm/bam.summary.txt" "${OUT}/norm/bam.t4.summary.txt"; then
    echo "TEST FAILED: -t 4 summary differs from -t 1 summary" >&2
    exit 1
fi
rm -f "${OUT}/norm/demultiplexed.t4.sam" "${OUT}/norm/bam.t4.summary.txt"

# One "md5  relative/path" line per normalized file, in a stable order.
: > "${OUT}/observed.md5"
( cd "${OUT}/norm" && find . -type f | LC_ALL=C sort ) | \
while read -r f; do
    printf '%s  %s\n' "$( md5of < "${OUT}/norm/${f#./}" )" "${f#./}" >> "${OUT}/observed.md5"
done

nfiles=$( wc -l < "${OUT}/observed.md5" )
if [ "$nfiles" -lt 2 ]; then
    echo "TEST FAILED: only $nfiles output file(s) were produced" >&2
    exit 1
fi

if [ "$UPDATE" -eq 1 ]; then
    cp "${OUT}/observed.md5" "$EXPECTED"
    echo "Wrote $nfiles md5sums to $EXPECTED"
    exit 0
fi

if [ ! -f "$EXPECTED" ]; then
    echo "TEST ERROR: $EXPECTED is missing, create it with \"make test-update\"" >&2
    exit 2
fi

if diff -u "$EXPECTED" "${OUT}/observed.md5" > "${OUT}/md5.diff" 2>&1; then
    echo "TEST PASSED: $nfiles output files match $EXPECTED"
    rm -rf "$OUT"
    exit 0
else
    echo "TEST FAILED: md5sums differ (expected vs observed)" >&2
    cat "${OUT}/md5.diff" >&2
    echo "Output kept in $OUT for inspection" >&2
    exit 1
fi
