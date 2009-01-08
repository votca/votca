#! /bin/bash

if [ "$1" = "--help" ]; then
   echo Start the script to run ibm, imc, etc.
   echo Usage: ${0##*/} setting_file
   exit 0
fi

if [ -z "$1" ]; then
   echo Error: Missing setting file > /dev/stderr
   exit 1
else 
   source $1
   shift
fi

#check if needed variables are set
for variable in p_target iterations filelist method sim_prog CSGSHARE scriptdir; do
   if [ -z "$(declare | grep -e "^$variable" )" ]; then
      echo Error $variable not definded, check setting file > /dev/stderr
      exit 1
   fi
done
#check if needed arrays are set
for array in atoms scheme pscheme; do
   if [ -z "$(declare | grep -e "^$array=(.*)$" )" ]; then
      echo Error $array not definded or is not an array, check setting file > /dev/stderr
      exit 1
   fi
done

export scriptdir
export CSGSHARE

#find source_wrapper.sh
#first local, then local scriptdir. then CSGSHARE
if [ -n  "${scriptdir}" ] && [ -f "${scriptdir}/source_wrapper.sh" ]; then
   source_wrapper="${scriptdir}/source_wrapper.sh"
elif [ -n  "${CSGSHARE}" ] && [ -f "${CSGSHARE}/source_wrapper.sh" ]; then
   source_wrapper="${CSGSHARE}/source_wrapper.sh"
else
   echo Could not find source_wrapper.sh > /dev/stderr
   exit 1
fi
export source_wrapper

do_external() {
   local script
   script="$($source_wrapper $1 $2)" || exit 1
   shift 2
   echo Running subscript ${script##*/} "$@"
   source $script "$@"
}

#useful subroutine check if a command was succesful AND log the output
run_or_exit() {
   local prog mylog
   [[ "$1" = "--log" ]] && { mylog="$2"; shift 2; }
   prog=$1
   shift
   [[ -n "$prog" ]] || { echo Error give one argument >&2; exit 1; }
   [[ -z "$mylog" ]] && mylog="log_${prog##*/}"
   echo Running ${prog##*/} $* \&\> $mylog
   $prog $* &> $mylog
   [[ $? -eq 0 ]] || { echo Error at ${prog##*/}; exit 1; }
}

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
   for ((a1=0;a1<${#atoms[*]};a1++)); do
      atom1=${atoms[$a1]}
      for ((j=$a1;a2<${#atoms[*]};a2++)); do
         atom2=${atoms[$a2]}
         cp ../rdf_${atom1}_${atom2}_aim.xvg . || exit 1
      done
   done

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
