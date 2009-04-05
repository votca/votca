#! /bin/bash

if [ "$1" = "--help" ]; then
  echo Start the script to run ibm, imc, etc.
  echo Usage: ${0##*/} setting_file.xml 
  exit 0
fi

if [ -z "$1" ]; then
  echo Error: Missing setting file > /dev/stderr
  exit 1
fi
if [ -f "./$1" ]; then
  export CSGXMLFILE="$PWD/$1"
else
  echo Error: file could not read > /dev/stderr
  exit 1
fi

#check for CSGSHARE 
if [ -z "${CSGSHARE}" ]; then
  echo Error: CSGSHARE not definded > /dev/stderr
  exit 1
fi
if [ ! -d "$CSGSHARE" ]; then
  echo "CSGSHARE '$CSGSHARE' is not a dir" > /dev/stderr
  exit 1
fi
#export CSGSHARE

#we need csg_property
if [ -z "$(type -p csg_property)" ]; then
  echo Error: csg_property not found, check your PATH > /dev/stderr
  exit 1
fi

CSGSCRIPTDIR="$(csg_property --file $CSGXMLFILE --path cg.inverse.scriptdir --short --print .)"
#scriptdir maybe contains $PWD or something
eval CSGSCRIPTDIR=$CSGSCRIPTDIR
if [ ! -d "$CSGSCRIPTDIR" ]; then
  echo "CSGSCRIPTDIR '$CSGSCRIPTDIR' is not a dir" > /dev/stderr
  exit 1
fi
export CSGSCRIPTDIR

#find source_wrapper.sh
#first local, then local scriptdir. then CSGSHARE
if [ -f "${CSGSCRIPTDIR}/source_wrapper.sh" ]; then
   SOURCE_WRAPPER="${CSGSCRIPTDIR}/source_wrapper.sh"
elif [ -f "${CSGSHARE}/source_wrapper.sh" ]; then
   SOURCE_WRAPPER="${CSGSHARE}/source_wrapper.sh"
else
   echo Could not find source_wrapper.sh > /dev/stderr
   exit 1
fi
export SOURCE_WRAPPER

source $($SOURCE_WRAPPER --direct inverse_functions.sh)

method="$(get_sim_property method)" || exit 1
echo "We are doing Method: $method"

sim_prog="$(get_sim_property program)" || exit 1
echo "We using Sim Program: $sim_prog"
source $($SOURCE_WRAPPER functions $sim_prog) || exit 1

iterations="$(get_sim_property iterations_max)" || exit 1
echo "We are doing $iterations iterations."

filelist="$(get_sim_property filelist)" || exit 1
echo "We extra need $filelist to run the simulation"

p_target="$(get_sim_property p_target)" || exit 1

#main script
[[ -f done ]] && echo Job is already done && exit 0

if [ -d step_00 ]; then
  echo Skiping prepare 
  if [ ! -f step_00/done ]; then
    echo Incomplete step 00
    exit 1
  fi
else
  echo ------------------------
  echo Prepare \(make step_00\)
  echo ------------------------
  mkdir step_00
  cd step_00

  #copy all rdf in step_00
  for_all non-bonded "cp ../rdf_\${type1}_\${type2}_aim.xvg ." || exit 1

  do_external init $method for_all non-bonded || exit 1

  #convert generic table file in gromacs table file (.xvg)
  do_external convert_potential $sim_prog for_all non-bonded "table_\${type1}_\${type2}_new.d" || exit 1
   
  #make confout.gro
  do_external init $sim_prog || exit 1

  for_all non-bonded cp "table_\${type1}_\${type2}_new.d ../table_\${type1}_\${type2}_final.d" || exit 1
  touch done
  cd ..
fi

for ((i=1;i<$iterations+1;i++)); do
   echo ---------------------------------
   echo Doing iteration $i \(make step_$i\)
   echo -------------------------------
   last_dir=$(printf step_%02i $((i-1)) )
   this_dir=$(printf step_%02i $i)
   if [ -d $this_dir ]; then
      if [ -f $this_dir/done ]; then
         echo is already done - skipping
         continue
      else
         echo Incomplete step $i
         exit 1
      fi
   fi
   mkdir $this_dir
   
   #copy all rdf in step_00
   for_all non-bonded "cp rdf_\${type1}_\${type2}_aim.xvg $this_dir" || exit 1
   
   #get need files
   for myfile in $filelist; do
      cp ./$myfile ./$this_dir/ || exit 1
   done
   cd $this_dir
   
   for_all non-bonded "cp ../$last_dir/table_\${type1}_\${type2}_new.d ./table_\${type1}_\${type2}.d" || exit 1
   do_external convert_potential $sim_prog for_all non-bonded "table_\${type1}_\${type2}.d" || exit 1

   #Run simulation maybe change to Espresso or whatever
   do_external prepare $sim_prog $last_dir || exit 1
   do_external run $sim_prog || exit 1

   p_now="$(do_external pressure $sim_prog)" || exit 1
   echo New pressure $p_now
   echo Target pressure was $p_target
   
   #calc pressure correction
   pressure_cor=$($SOURCE_WRAPPER --direct pressure_cor.pl) || exit 1
   run_or_exit $pressure_cor $p_target $p_now pressure_cor.d 

   do_external rdf $sim_prog for_all non-bonded
   
   do_external update $method for_all non-bonded $i

   #copy latest results
   for_all non-bonded "cp table_\${type1}_\${type2}_new.d ../table_\${type1}_\${type2}_final.d"
   touch done
   cd ..
done

touch done
exit 0
