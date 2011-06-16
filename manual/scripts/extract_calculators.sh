#!/bin/bash -e

CALCULATOR_SRC="../md2qm/src/libmd2qm/calculators"
out="$PWD/ref_calculators.tex"

if [ ! -d $CALCULATOR_SRC ]; then
  echo "error: calculator directory not found"
fi

cd $CALCULATOR_SRC

echo "running doxygen"

doxygen -g > doxygen.log 2>&1
doxygen >> doxygen.log 2>&1

echo "% Calculators" > $out

for file in latex/*.tex; do
  name=$(sed -ne '/\\subsection{Detailed Description}/,/\\subsection/p' $file | awk '/Callname/{print $2}')
  if [ ! -z "$name" ]; then
    echo "Generating documentation for calculator $name"
    echo "\subsection{$name}" >> $out
    sed -ne '/\\subsection{Detailed Description}/,/\\subsection/p' $file | sed -e '1d' -e '$d' -e '/Callname/d' >> $out
  fi
done



