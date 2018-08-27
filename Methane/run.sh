#!/bin/bash

#convienience function to change xml option
changeoption(){
    sed -i "s&<${1}.*>.*</${1}>&<${1}>${2}</${1}>&" $3
}
#convienience function to delete xml option
deleteoption(){
 sed -i "s&<${1}.*>.*</${1}>&&" $2
}

echo $VOTCASHARE

#runs the mapping from MD coordinates to segments and creates .sql file

xtp_map -t MD_FILES/topol.tpr -c MD_FILES/conf.gro -s system.xml -f state.sql

# you can explore the created .sql file with e.g. sqlitebrowser

# output MD and QM mappings into extract.trajectory_md.pdb and extract.trajectory_qm.pdb files

xtp_dump -e trajectory2pdb -f state.sql

#make OPTIONFILE folder, you can put all options into a single options.xml file but experience has shown, that this is better.

mkdir OPTIONFILES

#copy neighborlistfile from votca share folder to here

cp $VOTCASHARE/xtp/xml/neighborlist.xml OPTIONFILES/

changeoption constant 0.6 OPTIONFILES/neighborlist.xml

#run neighborlist calculator

xtp_run -e neighborlist -o OPTIONFILES/neighborlist.xml -f state.sql


# read in reorganisation energies stored in system.xml to state.sql

cp $VOTCASHARE/xtp/xml/einternal.xml OPTIONFILES/

xtp_run -e einternal -o OPTIONFILES/einternal.xml -f state.sql


#site energies

#setup jobfile xqmultipole has no own jobfile thus wehave to use the jobwriter

cp $VOTCASHARE/xtp/xml/jobwriter.xml OPTIONFILES/

changeoption keys "mps.monomer mps.background" OPTIONFILES/jobwriter.xml 
changeoption states "n e h" OPTIONFILES/jobwriter.xml 

xtp_run -e jobwriter -o OPTIONFILES/jobwriter.xml -f state.sql -s 0
mv jobwriter.mps.background.tab mps.tab
mv jobwriter.mps.monomer.xml xqmultipole.jobs
#Only run the first 3 jobs and set the rest to complete
sed -i "s/AVAILABLE/COMPLETE/g" xqmultipole.jobs
sed -i '0,/COMPLETE/s/COMPLETE/AVAILABLE/' xqmultipole.jobs
sed -i '0,/COMPLETE/s/COMPLETE/AVAILABLE/' xqmultipole.jobs
sed -i '0,/COMPLETE/s/COMPLETE/AVAILABLE/' xqmultipole.jobs


#running xqmultipole

cp $VOTCASHARE/ctp/xml/xqmultipole.xml OPTIONFILES/

changeoption job_file xqmultipole.jobs OPTIONFILES/xqmultipole.xml
changeoption emp_file mps.tab OPTIONFILES/xqmultipole.xml
#switch polarisation off for tutorial
changeoption induce 0 OPTIONFILES/xqmultipole.xml

changeoption pdb_check 0 OPTIONFILES/xqmultipole.xml
deleteoption write_chk OPTIONFILES/xqmultipole.xml
echo "Running xqmultipole, rerouting output to xqmultipole.log"
xtp_parallel -e xqmultipole -f state.sql -o OPTIONFILES/xqmultipole.xml -s 0 -t 1 -c 1000 -j "run" > xqmultipole.log

# xqmultipole has no parser to read the siteenergies into the sql file, write a python script or look at https://github.com/JensWehner/votca-scripts/blob/master/xtp/xtp_parseewald.py


#running eqm
#eqm runs qm calculations for each segment in the sql file, it consists of three stages first writing a jobfile, then running the calculations, if necessary 
# on multiple machines and then reading them into the .sql file
echo "Running eQM"

cp $VOTCASHARE/xtp/xml/eqm.xml OPTIONFILES/
cp $VOTCASHARE/xtp/packages/mbgft.xml OPTIONFILES/
cp $VOTCASHARE/xtp/packages/xtpdft.xml OPTIONFILES/
cp $VOTCASHARE/xtp/packages/esp2multipole.xml OPTIONFILES/

changeoption dftpackage OPTIONFILES/xtpdft.xml OPTIONFILES/eqm.xml
changeoption gwbse_options OPTIONFILES/mbgft.xml OPTIONFILES/eqm.xml

changeoption threads 0 OPTIONFILES/xtpdft.xml
changeoption openmp 0 OPTIONFILES/mbgft.xml
changeoption ranges full OPTIONFILES/mbgft.xml
changeoption openmp 0 OPTIONFILES/esp2multipole.xml

xtp_parallel -e eqm -o OPTIONFILES/eqm.xml -f state.sql -s 0 -j "write"
sed -i "s/AVAILABLE/COMPLETE/g" eqm.jobs
sed -i '0,/COMPLETE/s/COMPLETE/AVAILABLE/' eqm.jobs
sed -i '0,/COMPLETE/s/COMPLETE/AVAILABLE/' eqm.jobs

xtp_parallel -e eqm -o OPTIONFILES/eqm.xml -f state.sql -s 0 -j run -c 1 -t 1

#running iqm
#iqm runs qm calculations for each pair in the sql file, it consists of three stages first writing a jobfile, then running the calculations, if necessary 
# on multiple machines and then reading them into the .sql file
echo "Running iQM"

cp $VOTCASHARE/xtp/xml/iqm.xml OPTIONFILES/
cp $VOTCASHARE/xtp/packages/mbgft.xml OPTIONFILES/mbgft_pair.xml
cp $VOTCASHARE/xtp/packages/xtpdft.xml OPTIONFILES/xtpdft_pair.xml
cp $VOTCASHARE/xtp/packages/bsecoupling.xml OPTIONFILES/

changeoption bsecoupling_options OPTIONFILES/bsecoupling.xml OPTIONFILES/iqm.xml
changeoption dftpackage OPTIONFILES/xtpdft_pair.xml OPTIONFILES/iqm.xml
changeoption gwbse_options OPTIONFILES/mbgft_pair.xml OPTIONFILES/iqm.xml
changeoption read_guess 1 OPTIONFILES/xtpdft_pair.xml
changeoption energy 1e-2 OPTIONFILES/xtpdft_pair.xml

changeoption threads 0 OPTIONFILES/xtpdft_pair.xml
changeoption ranges full OPTIONFILES/mbgft_pair.xml
changeoption openmp 0 OPTIONFILES/mbgft_pair.xml
changeoption openmp 0 OPTIONFILES/bsecoupling.xml
changeoption openmp 0 OPTIONFILES/bsecoupling.xml

changeoption singlet "Methane:s1" OPTIONFILES/iqm.xml
changeoption triplet "Methane:t1" OPTIONFILES/iqm.xml
changeoption electron "Methane:e1" OPTIONFILES/iqm.xml
changeoption hole "Methane:h1" OPTIONFILES/iqm.xml

changeoption tasks "singlets,triplets,iqm" OPTIONFILES/mbgft_pair.xml

xtp_parallel -e iqm -o OPTIONFILES/iqm.xml -f state.sql -s 0 -j "write"
sed -i "s/AVAILABLE/COMPLETE/g" iqm.jobs
sed -i '0,/COMPLETE/s/COMPLETE/AVAILABLE/' iqm.jobs

xtp_parallel -e iqm -o OPTIONFILES/iqm.xml -f state.sql -s 0 -j run -c 1 -t 1

xtp_parallel -e iqm -o OPTIONFILES/iqm.xml -f state.sql -j "read"

#running qmmm 

cp $VOTCASHARE/xtp/xml/qmmm.xml OPTIONFILES/
cp $VOTCASHARE/xtp/packages/mbgft.xml OPTIONFILES/mbgft_qmmm.xml
cp $VOTCASHARE/xtp/packages/xtpdft.xml OPTIONFILES/xtpdft_qmmm.xml


#writing jobfile
cp $VOTCASHARE/xtp/xml/jobwriter.xml OPTIONFILES/jobwriter_qmmm.xml

changeoption keys "mps.monomer mps.background" OPTIONFILES/jobwriter_qmmm.xml 

xtp_run -e jobwriter -o OPTIONFILES/jobwriter_qmmm.xml -f state.sql -s 0
mv jobwriter.mps.background.tab mps.tab
mv jobwriter.mps.monomer.xml qmmm.jobs
sed -i "s/AVAILABLE/COMPLETE/g" qmmm.jobs
sed -i '0,/COMPLETE/s/COMPLETE/AVAILABLE/' qmmm.jobs
sed -i '0,/COMPLETE/s/COMPLETE/AVAILABLE/' qmmm.jobs

#config options
changeoption archiving iterations OPTIONFILES/qmmm.xml
changeoption dftpackage OPTIONFILES/xtpdft_qmmm.xml OPTIONFILES/qmmm.xml
changeoption gwbse_options OPTIONFILES/mbgft_qmmm.xml OPTIONFILES/qmmm.xml
changeoption cutoff1 1 OPTIONFILES/qmmm.xml
changeoption cutoff2 1.2 OPTIONFILES/qmmm.xml

changeoption threads 0 OPTIONFILES/xtpdft_qmmm.xml
changeoption openmp 0 OPTIONFILES/mbgft_qmmm.xml
changeoption ranges full OPTIONFILES/mbgft_qmmm.xml
echo "Running qmmm, rerouting output to qmmm.log"
xtp_parallel -e qmmm -o OPTIONFILES/qmmm.xml -f state.sql -s 0 >qmmm.log


















