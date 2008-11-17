#! /bin/bash

base_file_dir="base_files"
step_file="steps"
ignore="no"

long_step_name="./$base_file_dir/$step_file"

if [ -z "$1" ]; then
   $0 --help
   exit 1
fi

if [ "$1" == "--help" ]; then
   echo Usage: ${0##*/} \[show\|clean\|all\|rm\|NUMBER\]
   exit 1
fi

if [ "$1" == "--ignore" ]; then
   ignore="yes"
   shift
fi

if [ ! -r "$long_step_name" ]; then
   echo Step file \"$long_step_name\" not readable
   exit 1
fi

if [ "$1" == "show" ]; then
   echo Step order:
   cat -n $long_step_name
   exit
fi

if [ "$1" == "clean" ]; then
   echo Make clean
   rm -r [0-9]*
   exit
fi

if [ "$1" == "rm" ]; then
   if [ -n "${2//[0-9]}" ]; then
      echo Expected a number !
      exit 1
   fi
   echo Clean step $2
   step="$(awk "(NR==$2){print \$0}" $long_step_name)"
   dirname=$(printf "%02i_%s" $2 $step)
   rm -r $dirname
   exit
fi

make_step() {
      mydirname="$(printf "%02i_%s" $number $step)"
      echo Doing: $mydirname
      if [ -d $mydirname ]; then
         echo $mydirname already exist
         if [ -e $mydirname/done ]; then
            return
         else
            echo Nothing seem to be done - check by hand
            if [ "$ignore" = "no" ]; then
               exit 1
            fi
         fi
      fi
      mkdir $mydirname
      cp $base_file_dir/$step/*[^~] $mydirname
      cd $mydirname
      ls
      ./${step}.sh
      [[ $? -ne 0 ]] && echo Error at step $step && exit 1
      cd ..
}

if [ "$1" == "all" ]; then
   number=1
   for step in $(cat $long_step_name); do
      make_step
      ((number++))
   done
   exit
fi

if [ -n "${1//[0-9]}" ]; then
   echo Expected a number !
   exit 1
fi

step="$(awk "(NR==$1){print \$0}" $long_step_name)"
number=$1

if [ -z "$step" ]; then
   echo step not found
   exit 1
fi

make_step
