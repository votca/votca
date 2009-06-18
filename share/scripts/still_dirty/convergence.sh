#!/bin/bash

for step in step_?? step_???; do
  if [ -e $step/CG-CG.dist.new ]; then
    paste $step/CG-CG.dist.tgt $step/CG-CG.dist.new | sed -e 's/nan/0.0/' | awk '{sum+=($5-$2);}END{print sum;}'
  fi
done

