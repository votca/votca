#!/bin/bash

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

sed -i "s_<constant>.*</constant>_<constant>0.25</constant>_" OPTIONFILES/neighborlist.xml

#run neighborlist calculator

xtp_run -e neighborlist -o OPTIONFILES/neighborlist.xml -f state.sql


# read in reorganisation energies stored in system.xml to state.sql

cp $VOTCASHARE/xtp/xml/einternal.xml OPTIONFILES/

xtp_run -e einternal -o OPTIONFILES/einternal.xml -f state.sql


#site energies

#setup jobfile xqmultipole has no own jobfile thus wehave to use the jobwriter

cp $VOTCASHARE/xtp/xml/jobwriter.xml OPTIONFILES/

sed -i "s_<keys.*>.*</keys>_<keys>mps.chrg</keys>_" OPTIONFILES/jobwriter.xml 

xtp_run -e jobwriter -o OPTIONFILES/jobwriter.xml -f state.sql -s 0

mv jobwriter.mps.chrg.xml xqmultipole.jobs

#running xqmultipole














