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
   for_all --cmd non-bonded cp ../rdf_\${atom1}_\${atom2}_aim.xvg . || exit 1

   exit
   do_external init $method || exit 1

   #convert generic table file in gromacs table file (.xvg)
   for ((a1=0;a1<${#atoms[*]};a1++)); do
      atom1=${atoms[$a1]}
      for ((a2=$a1;a2<${#atoms[*]};a2++)); do
         atom2=${atoms[$a1]}
         do_external convert_potential $sim_prog table_${atom1}_${atom2}_new.d || exit 1
      done
   done
   
   #make confout.gro
   do_external init $sim_prog || exit 1

   for ((a1=0;a1<${#atoms[*]};a1++)); do
      atom1=${atoms[$a1]}
      for ((a2=$a1;a2<${#atoms[*]};a2++)); do
         atom2=${atoms[$a2]}
         cp table_${atom1}_${atom2}_new.d ../table_${atom1}_${atom2}_final.d || exit 1
      done
   done
   touch done
   cd ..
fi

for ((i=1;i<$iterations+1;i++)); do
   echo ---------------------------------
   echo Doing iteration $i \(make step_$i\)
   echo -------------------------------
   last=$i
   ((last--))
   last_dir=$(printf step_%02i $last)
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
   for ((a1=0;a1<${#atoms[*]};a1++)); do
      atom1=${atoms[$a1]}
      for ((a2=$a1;a2<${#atoms[*]};a2++)); do
         atom2=${atoms[$a2]}
         cp rdf_${atom1}_${atom2}_aim.xvg $this_dir || exit 1
      done
   done
   
   #get need files
   for myfile in $filelist; do
      cp ./$myfile ./$this_dir/ || exit 1
   done
   cd $this_dir
   
   for ((a1=0;a1<${#atoms[*]};a1++)); do
      atom1=${atoms[$a1]}
      for ((a2=$a1;a2<${#atoms[*]};a2++)); do
         atom2=${atoms[$a2]}
         cp ../$last_dir/table_${atom1}_${atom2}_new.d ./table_${atom1}_${atom2}.d || exit 1
         do_external convert_potential $sim_prog table_${atom1}_${atom2}.d || exit 1
      done
   done

   #Run simulation maybe change to Espresso or whatever
   do_external prepare $sim_prog
   do_external run $sim_prog

   do_external pressure $sim_prog
   echo New pressure $p_now
   echo Target pressure was $p_target
   
   #calc pressure correction
   pressure_cor=$($source_wrapper --direct pressure_cor.pl) || exit 1
   run_or_exit $pressure_cor $p_target $p_now pressure_cor.d 

   do_external rdf $sim_prog
   
   scheme_nr=$(( ( i - 1 ) % ${#scheme[@]} ))
   pscheme_nr=$(( ( i - 1 ) % ${#pscheme[@]} ))
   echo Doing Scheme step $scheme_nr : ${scheme[$scheme_nr]}
   echo Doing Pressure Scheme step $pscheme_nr : ${pscheme[$pscheme_nr]}
   for ((a1=0;a1<${#atoms[*]};a1++)); do
      atom1=${atoms[$a1]}
      for ((a2=$a1;a2<${#atoms[*]};a2++)); do
         atom2=${atoms[$a2]}
         #if atom is in this step
         if [ -n "${scheme[$scheme_nr]}" ] && [ -z "${scheme[$scheme_nr]/*${atom1}-${atom2}*}" ]; then
            echo Update ${atom1}-${atom2} - $method
            #update ibm
            do_external update $method $atom1 $atom2
            add_POT=$($source_wrapper --direct add_POT.pl) || exit 1
            run_or_exit --log log_add_POT_${atom1}_${atom2} $add_POT table_${atom1}_${atom2}.d delta_pot_${atom1}_${atom2}.d table_${atom1}_${atom2}_new1.d
         else
            echo Just copying ${atom1}-${atom2} - no $method
            cp table_${atom1}_${atom2}.d table_${atom1}_${atom2}_new1.d
         fi
         if [ -n "${pscheme[$pscheme_nr]}" ] && [ -z "${pscheme[$pscheme_nr]/*${atom1}-${atom2}*}" ]; then
            echo Update presuure ${atom1}-${atom2}
            add_POT=$($source_wrapper --direct add_POT.pl) || exit 1
            run_or_exit --log log_add_POT2_${atom1}_${atom2} $add_POT table_${atom1}_${atom2}_new1.d pressure_cor.d table_${atom1}_${atom2}_new.d
         else
            echo Just copying ${atom1}-${atom2} - no pressure correction
            cp  table_${atom1}_${atom2}_new1.d table_${atom1}_${atom2}_new.d
         fi
         #copy latest results
         cp table_${atom1}_${atom2}_new.d ../table_${atom1}_${atom2}_final.d
      done
   done
   touch done
   cd ..
done

touch done
exit 0
