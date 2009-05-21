#!/bin/bash
name="CG-CG"
method=$1
nblocks=16

echo "calculate error-bars"
export CSGINVERSE=$CSGSHARE/scripts/inverse
export SOURCE_WRAPPER=$CSGINVERSE/source_wrapper.pl
export CSGLOG=log.tmp
source $($SOURCE_WRAPPER functions common) || die "$SOURCE_WRAPPER functions common failed" 
export CSGXMLFILE=../cg.xml

average_tables() {
  outfile=$1
  tmp=$(mktemp)
  cp $2 $outfile
  shift 2
  n=1    
  while [ $1 ]; do
    paste $1 $outfile > $tmp
    awk "{printf(\"%s %.16e %s\\n\", \$1,($n*\$2 + \$5)/($n+1),\$3);}" $tmp > $outfile
    n=$((n+1))
    shift    
  done
  rm $tmp    
}

average_imc() { 
 octave $CSGINVERSE/imcdata_from_blocks.octave
 sed -e '/^[#@]/d' $name.gmc.block > $1.gmc
 sed -e '/^[#@]/d' $name.imc.block > $1.imc
}

calc_dpot_ibm() {
  update_POT="$($SOURCE_WRAPPER update ibm_pot)" \
    || die "${0##*/}: $SOURCE_WRAPPER update ibm_pot failed" 

  run_or_exit ${update_POT} ${name}.dist.tgt \
    $1.dist.new ${name}.pot.cur $1.dpot.new
}

all_dist=""
all_dpot=""
for block in $(seq 1 $nblocks); do
  echo "skipping block $block"
  all_dpot="$all_dpot ${name}_no_$block.dpot.new"
  all_dist="$all_dist ${name}_$block.dist.new"
  
  case $method in
    ibm)
      in_dist=""
      for i in $(seq 1 $nblocks | sed "/^${block}\$/d"); do
        in_dist="$in_dist ${name}_$i.dist.new"
      done
      #begin_block $block
      average_tables ${name}_no_$block.dist.new $in_dist
      calc_dpot_ibm ${name}_no_$block
      #end_block $block
    ;;
  imc)
    seq 1 $nblocks | sed "/^${block}\$/d" > $name.blocks
    average_imc  ${name}_no_$block
    $CSGINVERSE/solve_octave.sh ${name}_no_$block $name.pot.cur
    ;;
  esac
done

case $method in
  ibm)
    average_tables ${name}.dist.new $all_dist
    calc_dpot_ibm ${name}
    ;;
  imc)
    seq 1 $nblocks > $name.blocks
    average_imc ${name}
    $CSGINVERSE/solve_octave.sh ${name} $name.pot.cur
    ;;
esac


~/src/csg/scripts/csg_call.sh tables jackknife $name.dpot.err CG-CG.dpot.new $all_dpot
#case  "$method" in
#  ibm)
#    ;;
#  imc)
#    ;;
#esac
#run_or_exit ${update_POT} ${name}.dist.tgt ${name}.dist.new ${name}.pot.cur ${name}.dpot.tmp