#! /bin/bash

#defaults
usage="Usage: ${0##*/} setting_file.xml"
do_iterations=""
clean="no"

show_help () {
  cat << eof
Start the script to run ibm, imc, etc.
$usage

OPTIONS:
-N, --do-iterations N         only do N iterationso
    --clean                   clean out the PWD, dangerous


USES: csg_get_property date \$SOURCE_WRAPPER msg mkdir for_all do_external mark_done cp die is_done log run_or_exit csg_get_interaction_property date \$CSGLOG
NEEDS: cg.inverse.method cg.inverse.program cg.inverse.iterations_max cg.inverse.filelist name
eof
}

### begin parsing options
shopt -s extglob
while [ "${1#-}" != "$1" ]; do
 if [ "${1#--}" = "$1" ] && [ -n "${1:2}" ]; then
    #short opt with arguments here: fc
    if [ "${1#-[fc]}" != "${1}" ]; then
       set -- "${1:0:2}" "${1:2}" "${@:2}"
    else
       set -- "${1:0:2}" "-${1:2}" "${@:2}"
    fi
 fi
 case $1 in 
   --do-iterations)
    do_iterations=$2
    shift 2 ;;
   -[0-9]*)
    do_iterations=${1#-}
    shift ;;
   --clean)
    clean="yes"
    shift ;;
   -h | --help)
    show_help
    exit 0;;
  *)
   echo "Unknown option '$1'"
   exit 1;;
 esac
done

### end parsing options 

#do all start up checks option stuff
source "${0%.sh}_start.sh"  "$@" || exit 1
#shift away xml file
shift 1

#----------------End of pre checking--------------------------------
if [ -f "$CSGLOG" ]; then
  log "\n\n#################################"
  log "# Appending to existing logfile #"
  log "#################################\n\n"
  log "Sim started $(date)"
  echo "Appending to existing logfile $CSGLOG"
else
  echo For a more verbose log see: $CSGLOG
  #log is created in the next line
  echo "Sim started $(date)" > $CSGLOG || exit 1
fi

method="$(csg_get_property cg.inverse.method)" 
log "We are doing Method: $method"

sim_prog="$(csg_get_property cg.inverse.program)"
log "We using Sim Program: $sim_prog"
source $($SOURCE_WRAPPER functions $sim_prog) || die "$SOURCE_WRAPPER functions $sim_prog failed" 

iterations="$(csg_get_property cg.inverse.iterations_max)" 
log "We are doing $iterations iterations."

filelist="$(csg_get_property cg.inverse.filelist)"  
log "We extra need $filelist to run the simulation"

run_or_exit $SOURCE_WRAPPER --status
run_or_exit $SOURCE_WRAPPER --check

#main script
[[ ! -f done ]] || { msg "Job is already done"; exit 0; }

this_dir="$(get_stepname 0)"
if [ -d "$this_dir" ]; then
  msg "Skiping prepare"
  [[ -f $this_dir/done ]] || die "Incomplete step 00"
else
  msg ------------------------
  msg "Prepare (make $this_dir)"
  msg ------------------------
  mkdir -p $this_dir || die "mkdir -p $this_dir failed"

  cd $this_dir || die "cd $this_dir failed"

  #copy+resample all rdf in $this_dir
  for_all non-bonded do_external resample calc ..

  for_all "non-bonded" do_external init $method  

  #get confout.gro
  do_external init $sim_prog 

  for_all non-bonded cp '$(csg_get_interaction_property name).pot.new ..' 
  touch done
  msg "$this_dir done"
  cd ..
fi

for ((i=1;i<$iterations+1;i++)); do
  last_dir=$(get_stepname $((i-1)) )
  this_dir=$(get_stepname $i)
  msg -------------------------------
  msg "Doing iteration $i (make $this_dir)"
  msg -------------------------------
  if [ -d $this_dir ]; then
    if [ -f $this_dir/done ]; then
      msg "step $i is already done - skipping"
      continue
    else
      msg "Incomplete step $i"
      [[ -f "${this_dir}/${CSGRESTART}" ]] || die "No restart file found"
    fi
  else
    log "Step $i started at $(date)"
    mkdir -p $this_dir || die "mkdir -p $this_dir failed"
  fi

  cd $this_dir || die "cd $this_dir failed"

  if is_done "Initialize"; then
    msg "Initialization already done"
  else
    #copy+resample all rdf in this_dir 
    for_all non-bonded do_external resample calc ..

    #get need files
    for myfile in $filelist; do
      run_or_exit cp ../$myfile .  
    done

    #get new pot from last step and make it current potential 
    for_all non-bonded "cp ../$last_dir/\$(csg_get_interaction_property name).pot.new ./\$(csg_get_interaction_property name).pot.cur" 

    #convert potential in format for sim_prog
    for_all non-bonded do_external convert_potential $sim_prog

    #Run simulation maybe change to Espresso or whatever
    do_external prepare $sim_prog "../$last_dir" 
    mark_done "Initialize"
  fi

  if is_done "Simulation"; then
    msg "Simulation is already done"
  else
    msg "Simulation runs"
    do_external run $sim_prog 
    mark_done "Simulation"
  fi

  msg "Make update $method" 
  do_external update $method $i

  msg "Post update"
  do_external post update $i

  msg "Adding up potential"
  do_external add_pot $method

  msg "Post add"
  do_external post add $i

  #copy latest results
  for_all non-bonded 'cp $(csg_get_interaction_property name).pot.new ..'

  touch done
  msg "step $i done"
  cd ..

  if [ -n "$do_iterations" ]; then
    ((do_iterations--))
    if [ $do_iterations -lt 1 ] ; then
      log "Stopping at step $i, user requested to take some rest after this amount of iterations"
      exit 0
    fi
  fi
done

touch done
log "All done at $(date)"
exit 0

