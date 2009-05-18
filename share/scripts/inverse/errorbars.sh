#!/bin/bash

echo "calculate error-bars"

add_block() {
    paste $1 $2 | \
      sed -e '/^[#@]/d' | \
      awk '{for (i=1;i<=NF/2;i++){printf "%e ",$i+$(i+NF/2);}printf"\n";}' \
      > $3
}

begin_block_ibm() {
  rm -Rf $name.dist.avg
}

end_block_ibm() {
  update_POT="$($SOURCE_WRAPPER update ibm_pot)" \
    || die "${0##*/}: $SOURCE_WRAPPER update ibm_pot failed" 

  run_or_exit ${update_POT} ${name}.dist.tgt \
    ${name}_${1}.dist.avg ${name}.pot.cur ${name}_$block.dpot.new
}

process_block_ibm() {
  if [ -e ${name}.dpot.avg ]; then
    paste ${name}.dpot.tmp ${name}.dpot.avg awk '{print $1,$5+$,$3}'
  else
    cp ${name}.dpot.tmp ${name}.dpot.avg
  fi
  
}

name="CG-CG"
method="ibm"
nblocks=4

for block in $(seq 1 $nblocks); do
  echo "skipping block $block"
  
  begin_block $block

  for i in $(seq 1 $nblocks); do
    if [ "$i" = "$block" ]; then
      continue
    fi
    process_block $i
  done

  end_block $block
done

#case  "$method" in
#  ibm)
#    ;;
#  imc)
#    ;;
#esac
#run_or_exit ${update_POT} ${name}.dist.tgt ${name}.dist.new ${name}.pot.cur ${name}.dpot.tmp