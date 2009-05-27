#!/bin/bash
rm errors.ibm
rm errors.imc

for i in $(seq 1 250); do
  sed -e '/[uo]/d' ibm.err.$i | sed -e '1d' | awk "{sum+=\$4;n++;}END{print $i*16, sum/n}" >> errors.ibm
sed -e '/[uo]/d' imc.err.$i | sed -e '1d' | awk "{sum+=\$4;n++;}END{print $i*16, sum/n}" >> errors.imc
done
