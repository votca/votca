#! /bin/bash

echo "The input is the name of a pdb file name #1and#2.pdb. Inside this file molecule 1 and molecule 2 are present and each of the atoms in there are 
labelled NAMEAT-#1 or NAMEAT-#2"
echo "NB: at the moment you cannot specify chkpoint change by hand if needed, also _you_ gotta run the g03 jobs!!!"

method="b3lyp/6-311+g*"
namemol=$1
idmol1=$( echo $namemol | sed -e 's/.pdb//' -e 's/and/ /' | awk '{print $1}' )
idmol2=$( echo $namemol | sed -e 's/.pdb//' -e 's/and/ /' | awk '{print $2}' )
 
if [ ! -e $namemol ]
then
	echo $namemol not found
	exit 1
fi

echo $idmol1
echo $idmol2
if [ ! -e $idmol1"and"$idmol2 ] 
then
	mkdir $idmol1"and"$idmol2
fi
mv $namemol $idmol1"and"$idmol2
cd $idmol1"and"$idmol2

echo "MODEL 0" >   ${idmol1}".pdb"
grep "[A-Z]-${idmol1}" $namemol | sed 's/-'$idmol1'//' >>  ${idmol1}".pdb"
echo "END" >> ${idmol1}".pdb"

echo "MODEL 0" > ${idmol2}".pdb"
grep "[A-Z]-${idmol2}" $namemol | sed 's/-'$idmol2'//' >>  ${idmol2}".pdb"
echo "END" >> ${idmol2}".pdb"

babel -i pdb ${idmol1}".pdb" -o xyz ${idmol1} &> babel1.log
babel -i pdb ${idmol2}".pdb" -o xyz ${idmol2} &> babel2.log

echo "%nproc=1
%mem=100Mb
#p pop=full METHOD nosymm

autogen

0 1" > header


mkdir part1 
mkdir part2
mkdir dim



sed "s:METHOD:${method}:" header > part1/part1.com
cat $idmol1 >> part1/part1.com
awk '{printf "%s-Bq %f %f %f \n", $1, $2, $3,$4}' $idmol2 >> part1/part1.com
echo >>  part1/part1.com


sed "s:METHOD:${method}:" header > part2/part2.com
awk '{printf "%s-Bq %f %f %f \n", $1, $2, $3,$4}' $idmol1 >> part2/part2.com
cat $idmol2 >> part2/part2.com
echo >>  part2/part2.com


echo "%chk=dimer.chk" > dim/dim.com
sed "s:METHOD:${method}:" header >> dim/dim.com
cat $idmol1 >> dim/dim.com
cat $idmol2 >> dim/dim.com
echo >>  dim/dim.com
echo "--Link1--
%chk=dimer.chk
%mem=100Mb
#p geom(allcheck) guess(read,only) IOp(3/33=1) ${method} nosymm

" >> dim/dim.com

cd ..

tar -cvf $idmol1"and"$idmol2.tar $idmol1"and"$idmol2 &> /dev/null

rm -rf $idmol1"and"$idmol2
