#! /bin/bash
#
# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version %version%
postadd plot script, send a certain plot script to gnuplot

Usage: ${0##*/} infile outfile

Used external packages: gnuplot
EOF
   exit 0
fi

start_gnuplot_pipe() {
  eval "exec ${fd}> gnuplot_lock"
  if flock -n -x $fd; then
    rm -rf gnuplot_pipe
    mkfifo gnuplot_pipe
    while true; do
      if read <gnuplot_pipe; then
        echo -e "$REPLY"
        [ "$REPLY" = "exit" ] && break
      fi
    done | $gnuplot $opts
  fi
}

[[ -z $1 || -z $2 ]] && die "${0##*/}: Missing arguments"

[[ -f $2 ]] && die "${0##*/}: $2 is already there"
do_external postadd dummy "$1" "$2"

fd=$(csg_get_interaction_property inverse.post_add_options.plot.fd "8")
is_int "$fd" || die "${0##*/}: inverse.post_add_options.plot.fd should be a number, but I got $fd"

gnuplot=$(csg_get_property cg.inverse.gnuplot_bin "gnuplot")
[ -n "$(type -p $gnuplot)" ] || die "${0##*/}: gnuplot binary '$gnuplot' not found"

opts=$(csg_get_interaction_property --allow-empty inverse.post_add_options.plot.gnuplot_opts)

script=$(csg_get_interaction_property inverse.post_add_options.plot.script)
[ -f "$script" ] || die "${0##*/}: plot script '$script' is not there, did you forget to add it to cg.inverse.filelist?"

what_to_kill="$(csg_get_interaction_property --allow-empty inverse.post_add_options.plot.kill)"

msg "Plotting '$script' using $gnuplot"
if [ -z "${what_to_kill}" ]; then
  cd $(get_main_dir)
  start_gnuplot_pipe &
  #wait for gnuplot_pipe
  sleep 1
  cd - > /dev/null

  #gnuplot is in laststep_dir
  echo "cd '$PWD'" > $(get_main_dir)/gnuplot_pipe || die "piping to gnuplot_pipe failed"

  cat $script > $(get_main_dir)/gnuplot_pipe || die "piping to gnuplot_pipe failed"
else
  killall $what_to_kill
  $gnuplot $opts $script
fi
