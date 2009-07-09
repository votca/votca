#!/bin/bash
i=0
for step in step_?? step_???; do
  if [ -e $step/CG-CG.dist.new ]; then
    echo -n "$i "
    paste $step/CG-CG.dist.tgt $step/CG-CG.dist.new | sed -e 's/nan/0.0/' | awk '{old = new; new=$1; sum+=($5-$2)**2;i++;}END{print sqrt(sum*(new-old));}'
  fi
  i=$((i+1))
done

