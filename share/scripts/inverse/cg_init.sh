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
scheme=( CG CG )

#additional iteration for pressure
pinterations=30
#pressure update scheme
#pscheme=( "X1-X1 X2-X2" "X3-X3" "X2-X3" )
pscheme=( CG CG )

#useful subroutine check if a command was succesful AND log the output
run_or_exit() {
   local prog mylog
   [[ "$1" = "--log" ]] && { mylog="$2"; shift 2; }
   prog=$1
   shift
   [[ -n "$prog" ]] || { echo Error give one argument >&2; exit 1; }
   [[ -z "$mylog" ]] && mylog="log_${prog##*/}"
   echo Running $prog $* \&\> $mylog
   $prog $* &> $mylog
   [[ $? -eq 0 ]] || { echo Error at $prog; exit 1; }
}

get_from_mdp() {
   [[ -n "$1" ]] || { echo What?; exit 1;}
   sed -n -e "s#$1[[:space:]]*=[[:space:]]*\(.*\)\$#\1#p" grompp.mdp | sed -e 's#;.*##'
}

#main script
[[ -f done ]] && echo Job is already done && exit 0

if [ -d step_00 ]; then
   echo Skiping prepare 
   if [ -f step_00/done ]; then
      continue
   else
     echo Incomplete step 00
      exit 1
   fi
else
   echo Prepare
   mkdir step_00
   cd step_00
   cp ../conf.gro confout.gro
   for rdf in $(ls ../rdf_*_aim.xvg); do
      cp $rdf .
   done
   
   for ((i=0;i<${#atoms[*]};i++)); do
      atom1=${atoms[$i]}
      for ((j=$i;j<${#atoms[*]};j++)); do
         echo $j
         atom2=${atoms[$j]}
         cp ../rdf_${atom1}_${atom2}_aim.xvg . || exit 1
         if [ -f ../table_${atom1}_${atom2}_guess.d ]; then
            echo Using given table for $atom1-$atom2
            cp ../table_${atom1}_${atom2}_guess.d table_${atom1}_${atom2}_new.d || exit 1
         else
            # RDF_to_POT.pl just does log g(r) = extrapolation
            echo Using intial guess from RDF for ${atom1}-${atom2}
            run_or_exit  --log log_RDF_to_POT_${atom1}_${atom2} ../RDF_to_POT.pl rdf_${atom1}_${atom2}_aim.xvg table_${atom1}_${atom2}_new.d
         fi
         #convert generic table file in gromacs table file (.xvg)
         run_or_exit --log log_table_to_xvg_${atom1}_${atom2} ../table_to_xvg.pl table_${atom1}_${atom2}_new.d table_${atom1}_${atom2}_new.xvg
         cp table_${atom1}_${atom2}_new.d ../table_${atom1}_${atom2}_final.d
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
   
   #move output to input
   cp ../$last_dir/confout.gro ./conf.gro
   for ((i=0;i<${#atoms[*]};i++)); do
      atom1=${atoms[$i]}
      for ((j=$i;j<${#atoms[*]};j++)); do
         atom2=${atoms[$j]}
         cp ../$last_dir/table_${atom1}_${atom2}_new.d ./table_${atom1}_${atom2}.d
         cp ../$last_dir/table_${atom1}_${atom2}_new.xvg ./table_${atom1}_${atom2}.xvg
      done
   done

   #Run simulation maybe change to Espresso or whatever
   echo -e "a ${atoms[*]}\nq" | make_ndx -f conf.gro &> log_make_ndx
   run_or_exit grompp -v -n index.ndx
   run_or_exit mdrun -n index.ndx

   nsteps=$(get_from_mdp nsteps)
   dt=$(get_from_mdp dt)
   #20 % is warmup
   equi=$(awk "BEGIN{print 0.2*$nsteps*$dt}")
   echo equi = $equi
   
   echo Running g_energy
   echo "Pressure" | g_energy -b $equi &> log_g_energy
   [[ $? -ne 0 ]] && echo Error at running g_energy && exit 1
   p_now=$(awk '/^Pressure/{print $3}' log_g_energy)
   echo New pressure $p_now
   echo Target pressure was $p_target
   #calc pressure correction
   run_or_exit ../pressure_cor.pl $p_target $p_now pressure_cor.d
   

   #inv boltzmann scheme
   scheme_nr=$((i % ${#scheme[@]} - 1 ))
   echo Doing Scheme step $scheme_nr - update ${scheme[$scheme_nr]}
   for ((i=0;i<${#atoms[*]};i++)); do
      atom1=${atoms[$i]}
      for ((j=$i;j<${#atoms[*]};j++)); do
         atom2=${atoms[$j]}
         echo Running g_rdf for ${atom1}-${atom2}
         echo -e "${atom1}\n${atom2}" | g_rdf -b $equi -n index.ndx -bin 0.01 -o rdf_${atom1}_${atom2}.xvg &> log_g_rdf_${atom1}_${atom2}
         [[ $? -ne 0 ]] && echo Error at g_rdf ${atom1}-${atom2} && exit 1
         #if atom is in this step
         if [ -z "${scheme[$scheme_nr]/*${atom1}-${atom2}*}" ]; then
            echo Update ${atom1}-${atom2}
            run_or_exit --log log_update_POT_${atom1}_${atom2} ../update_POT.pl rdf_${atom1}_${atom2}_aim.xvg rdf_${atom1}_${atom2}.xvg delta_pot_${atom1}_${atom2}.d
            run_or_exit --log log_add_POT_${atom1}_${atom2} ../add_POT.pl table_${atom1}_${atom2}.d delta_pot_${atom1}_${atom2}.d table_${atom1}_${atom2}_new1.d
            run_or_exit --log log_add_POT2_${atom1}_${atom2} ../add_POT.pl table_${atom1}_${atom2}_new1.d pressure_cor.d table_${atom1}_${atom2}_new.d
         else
            echo Just copying ${atom1}-${atom2}
            cp table_${atom1}_${atom2}.d table_${atom1}_${atom2}_new.d
         fi
         run_or_exit --log log_table_to_xvg_${atom1}_${atom2} ../table_to_xvg.pl table_${atom1}_${atom2}_new.d table_${atom1}_${atom2}_new.xvg
         
         #copy latest results
         cp rdf_${atom1}_${atom2}.xvg ../rdf_${atom1}_${atom2}_final.xvg
         cp table_${atom1}_${atom2}_new.d ../table_${atom1}_${atom2}_final.d
         cp table_${atom1}_${atom2}_new.xvg ../table_${atom1}_${atom2}_final.xvg
      done
   done
   touch done
   cd ..
done

if [ -d step_p00 ]; then
   echo Skiping Pressure prepare
   if [ -f step_p00/done ]; then
      continue
   else
     echo Incomplete step p00
      exit 1
   fi
else
   echo Pressure Prepare
   mkdir step_p00
   for ((i=0;i<${#atoms[*]};i++)); do
      atom1=${atoms[$i]}
      for ((j=$i;j<${#atoms[*]};j++)); do
         atom2=${atoms[$j]}
         cp ./$this_dir/table_${atom1}_${atom2}_new.d ./step_p00
         cp ./$this_dir/table_${atom1}_${atom2}_new.xvg ./step_p00
      done
   done
   cp ./$this_dir/confout.gro ./step_p00
   touch step_p00/done
fi

for ((i=1;i<$piterations+1;i++)) ; do
   echo Doing Pressure iteration $i
   last=$i
   ((last--))
   last_dir=$(printf step_p%02i $last)
   this_dir=$(printf step_p%02i $i)
   if [ -d $this_dir ]; then
      echo Skiping Pressure iteration $i
      if [ -f $this_dir/done ]; then
         continue
      else
         echo Incomplete step $i
         exit 1
      fi
   fi
   date
   mkdir $this_dir
   for myfile in $filelist; do
      cp ./$myfile ./$this_dir/
   done
   cd $this_dir

   #move output to input
   cp ../$last_dir/confout.gro ./conf.gro
   for ((i=0;i<${#atoms[*]};i++)); do
      atom1=${atoms[$i]}
      for ((j=$i;j<${#atoms[*]};j++)); do
         atom2=${atoms[$j]}
         cp ../$last_dir/table_${atom1}_${atom2}_new.d ./table_${atom1}_${atom2}.d
         cp ../$last_dir/table_${atom1}_${atom2}_new.xvg ./table_${atom1}_${atom2}.xvg
      done
   done

   #Run simulation maybe change to Espresso or whatever
   echo -e "a ${atoms[*]}\nq" | make_ndx -f conf.gro &> log_make_ndx
   run_or_exit grompp -v -n index.ndx
   run_or_exit mdrun -n index.ndx

   nsteps=$(get_from_mdp nsteps)
   dt=$(get_from_mdp dt)
   #20 % is warmup
   equi=$(awk "BEGIN{print 0.2*$nsteps*$dt}")
   echo equi = $equi
   
   echo Running g_energy
   echo "Pressure" | g_energy -b $equi &> log_g_energy
   [[ $? -ne 0 ]] && echo Error at running g_energy && exit 1
   p_now=$(awk '/^Pressure/{print $3}' log_g_energy)
   echo New pressure $p_now
   echo Target pressure was $p_target
   #calc pressure correction
   run_or_exit ../pressure_cor.pl $p_target $p_now pressure_cor.d
   
   scheme_nr=$((i % ${#pscheme[@]} - 1))
   echo Doing Scheme step $scheme_nr - update ${pscheme[$scheme_nr]}
   for ((i=0;i<${#atoms[*]};i++)); do
      atom1=${atoms[$i]}
      for ((j=$i;j<${#atoms[*]};j++)); do
         atom2=${atoms[$j]}
         echo Running g_rdf for ${atom1}-${atom2}
         echo -e "${atom1}\n${atom2}" | g_rdf -b $equi -n index.ndx -bin 0.01 -o rdf_${atom1}_${atom2}.xvg &> log_g_rdf_${atom1}_${atom2}
         [[ $? -ne 0 ]] && echo Error at g_rdf ${atom1}-${atom2} && exit 1
         #if atom is in this step
         if [ -z "${scheme[$scheme_nr]/*${atom1}-${atom2}*}" ]; then
            echo Update ${atom1}-${atom2}
            run_or_exit --log log_add_POT_${atom1}_${atom2} ../add_POT.pl table_${atom1}_${atom2}.d pressure_cor.d table_${atom1}_${atom2}_new.d
         else
            echo Just copying ${atom1}-${atom2}
            cp table_${atom1}_${atom2}.d table_${atom1}_${atom2}_new.d
         fi
         run_or_exit --log log_table_to_xvg_${atom1}_${atom2} ../table_to_xvg.pl table_${atom1}_${atom2}_new.d table_${atom1}_${atom2}_new.xvg
         
         #copy latest results
         cp rdf_${atom1}_${atom2}.xvg ../rdf_${atom1}_${atom2}_final.xvg
         cp table_${atom1}_${atom2}_new.d ../table_${atom1}_${atom2}_final.d
         cp table_${atom1}_${atom2}_new.xvg ../table_${atom1}_${atom2}_final.xvg
      done
   done
   touch done
   cd ..
done

touch done
exit 0;
