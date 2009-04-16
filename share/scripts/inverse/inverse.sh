#! /bin/bash

if [ "$1" = "--help" ]; then
  echo Start the script to run ibm, imc, etc.
  echo Usage: ${0##*/} setting_file.xml 
  exit 0
fi

#for now we replace it later
die(){ echo "$*" > /dev/stderr; exit 1; }

[[ -n "$1" ]] || die_ "Error: Missing xml file"

if [ -f "./$1" ]; then
  export CSGXMLFILE="$PWD/$1"
else
  die "Error: file could not read"
fi

#check for CSGSHARE 
[[ -n "${CSGSHARE}" ]] || die "Error: CSGSHARE not definded"
[[ -d "$CSGSHARE" ]] || die "CSGSHARE '$CSGSHARE' is not a dir"

#we need csg_property
[[ -n "$(type -p csg_property)" ]] || die "Error: csg_property not found, check your PATH"

CSGSCRIPTDIR="$(csg_property --file $CSGXMLFILE --path cg.inverse.scriptdir --short --print .)" || \
  die "csg_property --file $CSGXMLFILE --path cg.inverse.scriptdir --short --print . failed" 
#scriptdir maybe contains $PWD or something

eval CSGSCRIPTDIR=$CSGSCRIPTDIR
[[ -d "$CSGSCRIPTDIR" ]] || die "CSGSCRIPTDIR '$CSGSCRIPTDIR' is not a dir"
export CSGSCRIPTDIR

#find source_wrapper.pl
#first local, then local scriptdir. then CSGSHARE
if [ -f "${CSGSCRIPTDIR}/source_wrapper.pl" ]; then
   SOURCE_WRAPPER="${CSGSCRIPTDIR}/source_wrapper.pl"
elif [ -f "${CSGSHARE}/source_wrapper.pl" ]; then
   SOURCE_WRAPPER="${CSGSHARE}/source_wrapper.pl"
else
   die "Could not find source_wrapper.pl"
fi
export SOURCE_WRAPPER

CSGLOG="$(csg_property --file $CSGXMLFILE --path cg.inverse.log_file --short --print .)" || die "Could not get logfile"
[[ -n "${CSGLOG}" ]] || die "logfile is empty"
CSGLOG="$PWD/$CSGLOG"
export CSGLOG

function_file=$($SOURCE_WRAPPER functions common) || die "$SOURCE_WRAPPER functions common failed"
#die() is overwritten here
source ${function_file} || exit 1
unset function_file

#----------------End of pre checking--------------------------------
echo For a more verbose log see: $CSGLOG
#log is created
echo "Sim started $(date)" > $CSGLOG || exit 1

method="$(get_sim_property method)" 
log "We are doing Method: $method"

sim_prog="$(get_sim_property program)"
echo $sim_prog
log "We using Sim Program: $sim_prog"
source $($SOURCE_WRAPPER functions $sim_prog) || die "$SOURCE_WRAPPER functions $sim_prog failed" 

iterations="$(get_sim_property iterations_max)" 
log "We are doing $iterations iterations."

filelist="$(get_sim_property filelist)"  
log "We extra need $filelist to run the simulation"

#main script
[[ ! -f done ]] || { msg "Job is already done"; exit 0; }

if [ -d step_00 ]; then
  msg "Skiping prepare"
  [[ -f step_00/done ]] || die "Incomplete step 00"
else
  msg ------------------------
  msg Prepare \(make step_00\)
  msg ------------------------
  mkdir step_00

  # TODO: check weather dir could be created!!!
  
  #copy+resample all rdf in step_00
  cd step_00 || die "cd step_00 failed"
  for_all non-bonded 'cp ../$($csg_get name).dist.tgt .'
  #for_all non-bonded $($SOURCE_WRAPPER --direct resample_to_calc.sh) step_00
  # cd step_00 || die "cd step_00 failed"

  do_external init $method for_all non-bonded 

  #get confout.gro
  do_external init $sim_prog 

  for_all non-bonded cp '$($csg_get name).pot.new ..' 
  touch done
  msg "step_00 done"
  cd ..
fi

for ((i=1;i<$iterations+1;i++)); do
  msg ---------------------------------
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
  mkdir $this_dir
  
  #copy+resample all rdf in step_00
  for_all non-bonded "cp \$(\$csg_get name).dist.tgt $this_dir" 
  # for_all non-bonded $($SOURCE_WRAPPER --direct resample_to_calc.sh) $this_dir
  
  #get need files
  for myfile in $filelist; do
    run_or_exit cp ./$myfile ./$this_dir/  
  done
  cd $this_dir || die "cd $this_dir failed"
 
  for_all non-bonded "cp ../$last_dir/\$(\$csg_get name).pot.new ./\$(\$csg_get name).pot.cur" 
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
  do_external post add

  #copy latest results
  for_all non-bonded 'cp $($csg_get name).pot.new ..'

  touch done
  msg "step $i done"
  cd ..
done

touch done
exit 0
