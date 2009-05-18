#!/bin/bash

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
    cp $outfile $outfile.$n
    paste $1 $outfile > $tmp
    awk "{printf(\"%s %.16e %s\\n\", \$1,($n*\$2 + \$5)/($n+1),\$3);}" $tmp > $outfile
    n=$((n+1))
    shift    
  done
  rm $tmp    
}

average_matrices() {
  outfile=$1
  tmp=$(mktemp)
  cp $2 $outfile
  shift 2
  n=1
  while [ $1 ]; do
    paste $1 $2 > $tmp
    awk "{for(i=1;i<=NF/2;i++) {printf(\"%e \", ($n*\$(i) + \$(i + NF/2))/($n+1));} printf(\"\n\")}" > $outfile
    n=$((n+1))
    shift    
  done
  rm $tmp    
}

calc_dpot_ibm() {
  update_POT="$($SOURCE_WRAPPER update ibm_pot)" \
    || die "${0##*/}: $SOURCE_WRAPPER update ibm_pot failed" 

  run_or_exit ${update_POT} ${name}.dist.tgt \
    $1.dist.new ${name}.pot.cur $1.dpot.new
}

name="CG-CG"
method="ibm"
nblocks=16

all_dist=""
all_dpot=""
for block in $(seq 1 $nblocks); do
  echo "skipping block $block"
  all_dpot="$all_dpot ${name}_no_$block.dpot.new"
  all_dist="$all_dist ${name}_$block.dist.new"

  in_dist=""
  for i in $(seq 1 $nblocks | sed "/^${block}\$/d"); do
    in_dist="$in_dist ${name}_$i.dist.new"
  done
  #begin_block $block
  average_tables ${name}_no_$block.dist.new $in_dist
  calc_dpot_ibm ${name}_no_$block
  #end_block $block
done

average_tables ${name}.dist.new $all_dist
calc_dpot_ibm ${name}

~/src/csg/scripts/csg_call.sh tables jackknife $name.dpot.err CG-CG.dpot.new $all_dpot
#case  "$method" in
#  ibm)
#    ;;
#  imc)
#    ;;
#esac
#run_or_exit ${update_POT} ${name}.dist.tgt ${name}.dist.new ${name}.pot.cur ${name}.dpot.tmp