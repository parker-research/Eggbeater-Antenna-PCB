#!/usr/bin/env bash
# Regenerate both boards, upgrade to the installed KiCad format, fill zones,
# run DRC and export JLCPCB fabrication files + preview images.
set -euo pipefail
cd "$(dirname "$0")/.."
FREQ="${1:-436.5}"

python3 scripts/generate_boards.py --freq-mhz "$FREQ"

for b in main_board phasing_board; do
  pcb="kicad/$b/$b.kicad_pcb"
  kicad-cli pcb upgrade "$pcb" >/dev/null
  kicad-cli pcb drc --refill-zones --save-board --severity-all --exit-code-violations \
    -o "kicad/$b/drc_report.txt" "$pcb"

  out="fab/$b"
  rm -rf "$out"; mkdir -p "$out/gerbers"
  if [ "$b" = phasing_board ]; then
    layers="F.Cu,In1.Cu,In2.Cu,B.Cu,F.Paste,B.Paste,F.SilkS,B.SilkS,F.Mask,B.Mask,Edge.Cuts"
  else
    layers="F.Cu,B.Cu,F.Paste,B.Paste,F.SilkS,B.SilkS,F.Mask,B.Mask,Edge.Cuts"
  fi
  kicad-cli pcb export gerbers --no-protel-ext -l "$layers" -o "$out/gerbers/" "$pcb" >/dev/null
  kicad-cli pcb export drill --format excellon --excellon-separate-th --generate-map --map-format gerberx2 \
    -o "$out/gerbers/" "$pcb" >/dev/null
  python3 -c "import shutil,sys; shutil.make_archive(sys.argv[1], 'zip', sys.argv[2])" "$out/${b}_gerbers" "$out/gerbers"

  mkdir -p docs/img
  kicad-cli pcb render --side top --quality high -w 1600 -h 1000 -o "docs/img/${b}_top.png" "$pcb" >/dev/null
  kicad-cli pcb render --side bottom --quality high -w 1600 -h 1000 -o "docs/img/${b}_bottom.png" "$pcb" >/dev/null
done
echo "done: fab/*/ *_gerbers.zip, docs/img/*.png"
