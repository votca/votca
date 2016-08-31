#! /bin/bash
echo "benzene orbital overlap"
xtp_overlap --conjseg benzene.xml --pos1 benzene1.pos --pos2 benzene2.pos --pdb benzene.pdb
echo "perylene orbital overlap"
xtp_overlap --conjseg perylene.xml --pos1 perylene1.pos --pos2 perylene2.pos --pdb perylene.pdb
