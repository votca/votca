#! /bin/bash

run_or_exit(){
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

[[ -f done ]] && echo Job is already done && exit 0

echo Initial files
number=${PWD##*/}
number=${number%%_*}
((number--))
number=$(printf "%02i" $number)
mydirname=$(cd ..;ls -d $number*)
cp ../$mydirname/confout.gro ./conf_ex.gro
for rdf in $(cd ../$mydirname; ls rdf_*.xvg); do
   aim_name=${rdf//OW/CG}
   cp ../$mydirname/$rdf ./${aim_name%.xvg}_aim.xvg
done
n_mol=$(perl -ne 'print $1 if /^SOL\s+(\S+?)$/' ../$mydirname/topol.top)
echo Found $n_mol SOL Molecules
n_salt=$(perl -ne 'print $1 if /^Cl\s+(\S+?)$/' ../$mydirname/topol.top)
echo Found $n_salt Na-Cl Molecules

p_target=$(awk '/^Pressure/{print $3}' ../$mydirname/log_g_energy)
echo found p_target $p_target
p_target=1
echo found p_target $p_target

if [ -d step_00 ]; then
   echo Skiping prepare
else
   echo Prepare
   mkdir step_00
   cd step_00
   cp ../conf_ex.gro .
   for rdf in $(ls ../rdf_*_aim.xvg); do
      cp $rdf .
   done
   total_mol=$(awk "BEGIN{print 2*$n_salt+$n_mol}")
   echo New total molecules $total_mol
   run_or_exit ../EX_to_CG.pl conf_ex.gro confout.gro
   sed -i "2s/.*/$total_mol/" confout.gro
   for atom in CG Na Cl; do
      if [ -f ../table_CG_${atom}_guess.d ]; then
         echo Using given table for CG-${atom}
         cp ../table_CG_${atom}_guess.d table_CG_${atom}_new.d
      else
         echo Using intial guess from RDF for CG-${atom}
         run_or_exit  --log log_RDF_to_POT_CG_${atom} ../RDF_to_POT.pl rdf_CG_${atom}_aim.xvg table_CG_${atom}_new.d
      fi
      n_lines=$(wc -l table_CG_${atom}_new.d | awk '{print 5*($1-1)}')
      echo Spline lines are $n_lines for CG-${atom}
      spline -n $n_lines table_CG_${atom}_new.d > smooth_table_CG_${atom}_new.d
      run_or_exit --log log_table_to_xvg_CG_${atom} ../table_to_xvg.pl smooth_table_CG_${atom}_new.d table_CG_${atom}_new.xvg
      cp smooth_table_CG_${atom}_new.d ../table_CG_${atom}_final.d
   done
   cd ..
fi

iterations=50
filelist="grompp.mdp rdf_*_aim.xvg cg_model.itp salt.itp topol.top table.xvg"
#5x cg then 5x Na and Cl
scheme=( CG CG CG CG CG "Na Cl" "Na Cl" "Na Cl" "Na Cl" "Na Cl" )
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
   date
   mkdir $this_dir
   for myfile in $filelist; do
      cp ./$myfile ./$this_dir/
   done
   cd $this_dir
   cp ../$last_dir/confout.gro ./conf.gro
   for atom in CG Na Cl; do
      cp ../$last_dir/table_CG_${atom}_new.d ./table_CG_${atom}.d
      cp ../$last_dir/table_CG_${atom}_new.xvg ./table_CG_${atom}.xvg
   done

   sed -i "s/N_MOL/$n_mol/" topol.top
   sed -i "s/N_SALT/$n_salt/" topol.top

   echo -e "a CG\nq" | make_ndx -f conf.gro &> log_make_ndx
   run_or_exit grompp -v -n index.ndx -maxwarn 1
   run_or_exit mdrun -v

   nsteps=$(perl -ne 'print $1 if /nsteps\s*=\s*(.*?)$/' grompp.mdp)
   dt=$(perl -ne 'print $1 if /dt\s*=\s*(.*?)$/' grompp.mdp)
   equi=$(awk "BEGIN{print 0.2*$nsteps*$dt}")
   echo equi = $equi
   
   echo Running g_energy
   echo "Pressure" | g_energy -b $equi &> log_g_energy
   [[ $? -ne 0 ]] && echo Error at running g_energy && exit 1
   p_now=$(awk '/^Pressure/{print $3}' log_g_energy)
   echo found p_now $p_now
   echo p_target was $p_target
   run_or_exit ../pressure_cor.pl $p_target $p_now pressure_cor.d
   
   scheme_nr=$((i % ${#scheme[@]} - 1 ))
   echo Doing Scheme step $scheme_nr - update ${scheme[$scheme_nr]}
   for atom in CG Na Cl; do
      echo Running g_rdf for CG-$atom
      echo -e "CG\n${atom}" | g_rdf -b $equi -n index.ndx -bin 0.01 -o rdf_CG_${atom}.xvg &> log_g_rdf_CG_${atom}
      [[ $? -ne 0 ]] && echo Error at g_rdf CG-${atom} && exit 1
      #if atom is in this step
      if [ -z "${scheme[$scheme_nr]/*$atom*}" ]; then
         echo Update CG-$atom
         run_or_exit --log log_update_POT_CG_${atom} ../update_POT.pl rdf_CG_${atom}_aim.xvg rdf_CG_${atom}.xvg delta_pot_CG_${atom}.d
         run_or_exit --log log_add_POT_CG_${atom} ../add_POT.pl table_CG_${atom}.d delta_pot_CG_${atom}.d table_CG_${atom}_new1.d
         run_or_exit --log log_add_POT2_CG_${atom} ../add_POT.pl table_CG_${atom}_new1.d pressure_cor.d table_CG_${atom}_new.d
      else
         echo Just copying CG-$atom
         cp table_CG_${atom}.d table_CG_${atom}_new.d
      fi
      n_lines=$(wc -l table_CG_${atom}_new.d | awk '{print 5*($1-1)}')
      echo Spline lines are $n_lines for CG-$atom
      spline -n $n_lines table_CG_${atom}_new.d > smooth_table_CG_${atom}_new.d
      run_or_exit --log log_table_to_xvg_CG_${atom} ../table_to_xvg.pl smooth_table_CG_${atom}_new.d table_CG_${atom}_new.xvg
      cp smooth_table_CG_${atom}_new.d ../table_CG_${atom}_final.d
      cp table_CG_${atom}_new.xvg ../table_CG_${atom}_final.xvg
   done
   touch done
   cd ..
done

piterations=30
if [ -d step_p00 ]; then
   echo Skiping Pressure prepare
else
   echo Pressure Prepare
   mkdir step_p00
   cp ./$this_dir/confout.gro ./step_p00
   for atom in CG Na Cl; do
      cp ./$this_dir/table_CG_${atom}_new.d ./step_p00
      cp ./$this_dir/table_CG_${atom}_new.xvg ./step_p00
   done
fi

#5x cg then 5x Na and Cl
pscheme=( CG CG CG CG CG "Na Cl" "Na Cl" "Na Cl" "Na Cl" "Na Cl" )
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
   cp ../$last_dir/confout.gro ./conf.gro
   for atom in CG Na Cl; do
      cp ../$last_dir/table_CG_${atom}_new.d ./table_CG_${atom}.d
      cp ../$last_dir/table_CG_${atom}_new.xvg ./table_CG_${atom}.xvg
   done

   sed -i "s/N_MOL/$n_mol/" topol.top
   sed -i "s/N_SALT/$n_salt/" topol.top

   echo -e "a CG\nq" | make_ndx -f conf.gro &> log_make_ndx
   run_or_exit grompp -v -n index.ndx -maxwarn 1
   run_or_exit mdrun -v

   nsteps=$(perl -ne 'print $1 if /nsteps\s*=\s*(.*?)$/' grompp.mdp)
   dt=$(perl -ne 'print $1 if /dt\s*=\s*(.*?)$/' grompp.mdp)
   equi=$(awk "BEGIN{print 0.2*$nsteps*$dt}")
   echo equi = $equi
   
   echo Running g_energy
   echo "Pressure" | g_energy -b $equi &> log_g_energy
   [[ $? -ne 0 ]] && echo Error at running g_energy && exit 1
   p_now=$(awk '/^Pressure/{print $3}' log_g_energy)
   echo found p_now $p_now
   echo p_target was $p_target
   run_or_exit ../pressure_cor.pl $p_target $p_now pressure_cor.d
   
   scheme_nr=$((i % ${#pscheme[@]} - 1))
   echo Doing Scheme step $scheme_nr - update ${pscheme[$scheme_nr]}
   for atom in CG Na Cl; do
      echo Running g_rdf for CG-$atom
      echo -e "CG\n${atom}" | g_rdf -b $equi -n index.ndx -bin 0.01 -o rdf_CG_${atom}.xvg &> log_g_rdf_CG_${atom}
      [[ $? -ne 0 ]] && echo Error at g_rdf CG-${atom} && exit 1
      #if atom is in this step
      if [ -z "${pscheme[$scheme_nr]/*$atom*}" ]; then
         echo Update CG-$atom
         run_or_exit --log log_add_POT_CG_${atom} ../add_POT.pl table_CG_${atom}.d pressure_cor.d table_CG_${atom}_new.d
      else
         echo Just copying CG-$atom
         cp table_CG_${atom}.d table_CG_${atom}_new.d
      fi
      n_lines=$(wc -l table_CG_${atom}_new.d | awk '{print 5*($1-1)}')
      echo Spline lines are $n_lines for CG-$atom
      spline -n $n_lines table_CG_${atom}_new.d > smooth_table_CG_${atom}_new.d
      run_or_exit --log log_table_to_xvg_CG_${atom} ../table_to_xvg.pl smooth_table_CG_${atom}_new.d table_CG_${atom}_new.xvg
      cp smooth_table_CG_${atom}_new.d ../table_CG_${atom}_final.d
   done
   
   touch done
   cd ..
done

touch done
exit 0;
