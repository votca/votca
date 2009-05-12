#! /bin/bash

if [ "$1" = "--help" ]; then
  echo Start the script to run ibm, imc, etc.
  echo Usage: ${0##*/} setting_file.xml 
  exit 0
fi

#do all start up checks
source "${0%.sh}_start.sh"  "$@" || exit 1

#----------------End of pre checking--------------------------------
echo For a more verbose log see: $CSGLOG
#log is created
echo "Sim started $(date)" > $CSGLOG || exit 1

method="$(get_sim_property method)" 
log "We are doing Method: $method"

sim_prog="$(get_sim_property program)"
log "We using Sim Program: $sim_prog"
source $($SOURCE_WRAPPER functions $sim_prog) || die "$SOURCE_WRAPPER functions $sim_prog failed" 

iterations="$(get_sim_property iterations_max)" 
log "We are doing $iterations iterations."

filelist="$(get_sim_property filelist)"  
log "We extra need $filelist to run the simulation"

run_or_exit $SOURCE_WRAPPER --status

#main script
[[ ! -f done ]] || { msg "Job is already done"; exit 0; }

if [ -d step_00 ]; then
  msg "Skiping prepare"
  [[ -f step_00/done ]] || die "Incomplete step 00"
else
  msg ------------------------
  msg Prepare \(make step_00\)
  msg ------------------------
  mkdir step_00 || die "mkdir step_00 failed"

  #copy+resample all rdf in step_00
  do_external resample calc for_all non-bonded step_00
  cd step_00 || die "cd step_00 failed"

  do_external init $method for_all non-bonded 

  #get confout.gro
  do_external init $sim_prog 

  for_all non-bonded cp '$($csg_get name).pot.new ..' 
  touch done
  msg "step_00 done"
  cd ..
fi

for ((i=1;i<$iterations+1;i++)); do
  msg -------------------------------
  msg Doing iteration $i \(make step_$i\)
  msg -------------------------------
  last_dir=$(printf step_%02i $((i-1)) )
  this_dir=$(printf step_%02i $i)
  if [ -d $this_dir ]; then
    if [ -f $this_dir/done ]; then
      msg "step $i is already done - skipping"
      continue
    else
      die "Incomplete step $i"
    fi
  fi
  log "Step $i started at $(date)"
  mkdir $this_dir || die "mkdir $this_dir failed"
  
  #copy+resample all rdf in this_dir 
  do_external resample calc for_all non-bonded $this_dir
  
  #get need files
  for myfile in $filelist; do
    run_or_exit cp ./$myfile ./$this_dir/  
  done
  cd $this_dir || die "cd $this_dir failed"

  #get new pot from last step and make it current potential 
  for_all non-bonded "cp ../$last_dir/\$(\$csg_get name).pot.new ./\$(\$csg_get name).pot.cur" 

  #convert potential in format for sim_prog
  do_external convert_potential $sim_prog for_all non-bonded

  #Run simulation maybe change to Espresso or whatever
  do_external prepare $sim_prog "../$last_dir" 

  msg "Simualtion runs"
  do_external run $sim_prog 

  msg "Make update $method" 
  do_external update $method $i

  msg "Post update"
  do_external post update $i

  msg "Adding up potential"
  do_external add_pot $method

  msg "Post add"
  do_external post add $i

  #copy latest results
  for_all non-bonded 'cp $($csg_get name).pot.new ..'

  touch done
  msg "step $i done"
  cd ..
done

touch done
log "All done at $(date)"
exit 0

