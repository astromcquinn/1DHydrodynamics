#!/bin/bash
# Build the gPLUTO executable for this problem on a CPU (g++), no MPI.
#   PLUTO_DIR=/path/to/gPLUTO ./build_gpluto.sh  [build_dir]
# gPLUTO: https://gitlab.com/PLUTO-code/gPLUTO  (tested with commit 870aa80, June 2026)
set -e
: "${PLUTO_DIR:?set PLUTO_DIR to the gPLUTO source directory}"
HERE=$(cd "$(dirname "$0")" && pwd)
B=${1:-$HERE/build}
mkdir -p "$B"
cp "$HERE"/definitions.hpp "$HERE"/init.cpp "$HERE"/radiat.cpp "$HERE"/pluto.ini "$B"/
cd "$B"
# pre-seed the architecture so setup.py does not open its interactive menu
printf "ARCH         = linux.g++.defs\n" > makefile
python3 "$PLUTO_DIR"/setup.py --auto-update --no-curses > setup.log
make -j4 > make.log
echo "built $B/pluto"
