#!/bin/bash

i=0
for step in step_?? step_???; do
  if [ -f $step/CG-CG.dist.new ]; then
    paste $step/CG-CG.dist.tgt $step/CG-CG.dist.new | sed -e 's/nan/0.0/' | awk "{sum+=(\$5-\$2);}END{print $i,sum;}"
  fi
  i=$((i+1))
done

