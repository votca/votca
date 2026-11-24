"""Input of the qmmm_geometry integration test.

Takes the QM geometry (job frame) and the final orbitals of the
qmmm_relaxed run and writes them as geometry.xyz and guess.orb, the
qmregion geometry/guess inputs of the next step of a geometry
optimization. Run with the qmmm_relaxed run directory and the target
directory as arguments.
"""
import os
import shutil
import sys
import xml.etree.ElementTree as ET

source, target = sys.argv[1], sys.argv[2]
root = ET.parse(os.path.join(source, "qmmm_relaxed.jobs")).getroot()
geometry = root.find(".//geometry")
if geometry is None:
    sys.exit("no geometry block in the qmmm_relaxed output")
lines = []
for atom in geometry.findall("atom"):
    x, y, z = (float(v) for v in atom.text.split())
    lines.append(f"{atom.get('element')} {x:.10f} {y:.10f} {z:.10f}")
with open(os.path.join(target, "geometry.xyz"), "w") as f:
    f.write(f"{len(lines)}\nqmmm_relaxed geometry\n" + "\n".join(lines) + "\n")
orbs = [os.path.join(d, n) for d, _, files in os.walk(source)
        for n in files if n == "qmregion_final.orb"]
if len(orbs) != 1:
    sys.exit(f"expected one qmregion_final.orb, found {len(orbs)}")
shutil.copy(orbs[0], os.path.join(target, "guess.orb"))
