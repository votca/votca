#!/bin/bash
# Generates the libint2 export tarball that BUILD_OWN_LIBINT builds from.
# Run once per libint version by a maintainer (not per machine, not per build);
# upload the tarball as a VOTCA release asset and put its URL and SHA256 into
# the CI / build instructions (-DVOTCA_LIBINT_EXPORT_URL, -DVOTCA_LIBINT_EXPORT_SHA256).
#
# Usage: make_libint_export.sh [VERSION] [OUTDIR]     (default 2.13.1, .)
# Needs: cmake >= 3.19, a C++ compiler, Boost headers, GMP (with C++ bindings),
#        Eigen3, internet access to github.com. The export step is serial and
#        can take several hours for these settings.
#
# Generator settings (only these need derivatives; 4-centre ERIs stay at
# energies only, XTP uses RI):
#   1-body, 2-centre and 3-centre integrals: 1st derivatives (analytic forces)
#   max AM 6 for orbital basis functions, 7 on the auxiliary centre of 2/3-centre
#   integrals (as in the conda-forge build), for energies and 1st derivatives
set -euo pipefail
# no AppleDouble (._*) files in tarballs made on macOS
export COPYFILE_DISABLE=1
VER=${1:-2.13.1}
mkdir -p "${2:-.}"; OUT=$(cd "${2:-.}" && pwd)
WORK=$(mktemp -d)
trap 'echo "work dir kept: $WORK"' ERR
cd "$WORK"
curl -fsSL -o src.tar.gz "https://github.com/evaleev/libint/archive/refs/tags/v${VER}.tar.gz"
tar xzf src.tar.gz
cmake -S "libint-${VER}" -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DLIBINT2_ENABLE_ONEBODY=1 \
  -DLIBINT2_ENABLE_ERI=0 \
  -DLIBINT2_ENABLE_ERI2=1 \
  -DLIBINT2_ENABLE_ERI3=1 \
  -DLIBINT2_MAX_AM="6;6" \
  -DLIBINT2_ERI2_MAX_AM="7;7" \
  -DLIBINT2_ERI3_MAX_AM="7;7"
cmake --build build --target export
TGZ=$(find build -maxdepth 3 -name "libint-${VER}*.tgz" | head -1)
[ -n "$TGZ" ] || { echo "no export tarball found under $WORK/build"; exit 1; }
NAME="libint-${VER}-votca-deriv1.tgz"
cp "$TGZ" "$OUT/$NAME"
# sanity check: the export's config.h must carry the derivative orders
mkdir chk && tar xzf "$OUT/$NAME" -C chk
grep -hE "define LIBINT_INCLUDE_(ONEBODY|ERI2|ERI3) " chk/*/include/libint2/config.h | tee /dev/stderr \
  | awk '$3 < 1 {bad=1} END {exit bad}' || { echo "export lacks derivatives"; exit 1; }
SHA=$(sha256sum "$OUT/$NAME" 2>/dev/null | cut -d' ' -f1 || shasum -a 256 "$OUT/$NAME" | cut -d' ' -f1)
rm -rf "$WORK"
echo
echo "export: $OUT/$NAME"
echo "  -DBUILD_OWN_LIBINT=ON -DVOTCA_LIBINT_EXPORT_URL=<url of $NAME> -DVOTCA_LIBINT_EXPORT_SHA256=$SHA"
