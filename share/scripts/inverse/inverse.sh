#! /bin/bash


#FILES we need

#main dir
#conf.gro -- starting configuration
#grompp.mdp -- grompp seeting
#topol.top -- topology file
#rdf_atom1_atom2_aim.xvg -- aim rdfs

#if you have 
#table_atom1_atom2_guess.d -- generic table (2columns r v(r))


#DEFAULTS

#target pressure_cor 
p_target=1
#grep this from g_energy out with:
#p_target=$(awk '/^Pressure/{print $3}' log_g_energy)

#atomsname here
#atoms=( X1 X2 X3 )
atoms=( CG )

#number of iterations
iterations=50
#file to run simulation
filelist="grompp.mdp rdf_*_aim.xvg topol.top table.xvg"
#update scheme
#scheme=( "X1-X1 X2-X2" "X3-X3" "X2-X3" )
scheme=( CG-CG "" )

#additional iteration for pressure
pinterations=30
#pressure update scheme
#pscheme=( "X1-X1 X2-X2" "X3-X3" "X2-X3" )
pscheme=( "" CG-CG )

method="ibm"
sim_prog="gromacs"
export scriptdir="$PWD"
source_wrapper="$scriptdir/source_wrapper.sh"

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
   echo Prepare \(make step_00\)
   mkdir step_00
   cd step_00

   #copy all rdf in step_00
   for ((i=0;i<${#atoms[*]};i++)); do
      atom1=${atoms[$i]}
      for ((j=$i;j<${#atoms[*]};j++)); do
         atom2=${atoms[$j]}
         cp ../rdf_${atom1}_${atom2}_aim.xvg . || exit 1
      done
   done

   do_external init $method || exit 1

   #convert generic table file in gromacs table file (.xvg)
   for ((i=0;i<${#atoms[*]};i++)); do
      atom1=${atoms[$i]}
      for ((j=$i;j<${#atoms[*]};j++)); do
         atom2=${atoms[$j]}
         do_external convert_potential $sim_prog table_${atom1}_${atom2}_new.d || exit 1
      done
   done
   
   #make confout.gro
   do_external init $sim_prog || exit 1

   for ((i=0;i<${#atoms[*]};i++)); do
      atom1=${atoms[$i]}
      for ((j=$i;j<${#atoms[*]};j++)); do
         atom2=${atoms[$j]}
         cp table_${atom1}_${atom2}_new.d ../table_${atom1}_${atom2}_final.d || exit 1
      done
   done
   touch done
   cd ..
fi

for ((i=1;i<$iterations+1;i++)) ; do
   echo Doing iteration $i
   last=$i
   ((last--))
   last_dir=$(printf step_%02i $last)
   this_dir=$(printf step_%02i $i)
   if [ -d $this_dir ]; then
      if [ -f $this_dir/done ]; then
         continue
      else
         echo Incomplete step $i
         exit 1
      fi
   fi
   mkdir $this_dir
   
   #get need files
   for myfile in $filelist; do
      cp ./$myfile ./$this_dir/
   done
   cd $this_dir
   
   for ((i=0;i<${#atoms[*]};i++)); do
      atom1=${atoms[$i]}
      for ((j=$i;j<${#atoms[*]};j++)); do
         atom2=${atoms[$j]}
         cp ../$last_dir/table_${atom1}_${atom2}_new.d ./table_${atom1}_${atom2}.d
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
   run_or_exit ../pressure_cor.pl $p_target $p_now pressure_cor.d

   do_external rdf $sim_prog
   
   scheme_nr=$((i % ${#scheme[@]} - 1 ))
   echo Doing Scheme step $scheme_nr - update ${scheme[$scheme_nr]}
   echo Doing Pressure Scheme step $scheme_nr - update ${pscheme[$scheme_nr]}
   for ((i=0;i<${#atoms[*]};i++)); do
      atom1=${atoms[$i]}
      for ((j=$i;j<${#atoms[*]};j++)); do
         atom2=${atoms[$j]}
         #if atom is in this step
         if [ -n "${scheme[$scheme_nr]}" ] && [ -z "${scheme[$scheme_nr]/*${atom1}-${atom2}*}" ]; then
            echo Update ${atom1}-${atom2} - $method
            #update ibm
            do_external update $method $atom1 $atom2
            run_or_exit --log log_add_POT_${atom1}_${atom2} ../add_POT.pl table_${atom1}_${atom2}.d delta_pot_${atom1}_${atom2}.d table_${atom1}_${atom2}_new1.d
            run_or_exit --log log_add_POT2_${atom1}_${atom2} ../add_POT.pl table_${atom1}_${atom2}_new1.d pressure_cor.d table_${atom1}_${atom2}_new.d
         else
            echo Just copying ${atom1}-${atom2} - no $method
            cp table_${atom1}_${atom2}.d table_${atom1}_${atom2}_new.d
         fi
         if [ -n "${pscheme[$scheme_nr]}" ] && [ -z "${pscheme[$scheme_nr]/*${atom1}-${atom2}*}" ]; then
            echo Update presuure ${atom1}-${atom2}
            run_or_exit --log log_add_POT2_${atom1}_${atom2} ../add_POT.pl table_${atom1}_${atom2}_new1.d pressure_cor.d table_${atom1}_${atom2}_new.d
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
exit 0;
