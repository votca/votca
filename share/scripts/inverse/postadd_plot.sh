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

Usage: ${0##*/}

Used external packages: gnuplot
EOF
   exit 0
fi

gnuplot=$(csg_get_property cg.inverse.gnuplot.bin)
[ -n "$(type -p $gnuplot)" ] || die "${0##*/}: gnuplot binary '$gnuplot' not found"

opts=$(csg_get_interaction_property --allow-empty inverse.post_add_options.plot.gnuplot_opts)

script=$(csg_get_interaction_property inverse.post_add_options.plot.script)
[[ -f $(get_main_dir)/$script ]] || die "${0##*/}: plot script '$script' is not in maindir"

what_to_kill="$(csg_get_interaction_property --allow-empty inverse.post_add_options.plot.kill)"

msg "Plotting '$script' using $gnuplot"
if [[ -z ${what_to_kill} ]]; then
  [[ -z $(type -p flock) ]] && die "${0##*/}: could not find flock needed by gnuplot pipe"
  [[ -z $(type -p mkfifo) ]] && die "${0##*/}: could not find mkfifo needed by gnuplot pipe"
  (
    flock -n 7 || exit 0 #see flock man page to understand this trick
    cd $(get_main_dir)
    rm -rf gnuplot_pipe gnuplot_pipe.log
    msg "Creating gnuplot_pipe ..."
    mkfifo gnuplot_pipe
    while true; do
      if read <gnuplot_pipe; then
        echo -e "$REPLY"
        echo -e "$REPLY" >> gnuplot_pipe.log
        [[ $REPLY = "exit" ]] && break
      fi
    done | $gnuplot $opts &
    while true; do
      if [[ -z $(ps -o pid= -p "${CSG_MASTER_PID}") ]]; then
	echo "exit" > $(get_main_dir)/gnuplot_pipe
	rm -rf gnuplot_pipe gnuplot_pipe.log gnuplot_pipe.lock
        exit
      fi
      sleep 1 #lowers the load
    done &
    sleep 1 #wait for gnuplot_pipe
    cd - > /dev/null
  ) 7> $(get_main_dir)/gnuplot_pipe.lock

  #gnuplot is in laststep, move to current one
  echo "cd '$PWD'" > $(get_main_dir)/gnuplot_pipe || die "piping to gnuplot_pipe failed"

  #name pipe accept only one command at the time, for i in $(cat ); do echo $i > pipe; done would do the same
  echo "load '$(get_main_dir)/$script'" > $(get_main_dir)/gnuplot_pipe || die "piping to gnuplot_pipe failed"
else
  [[ -z $(type -p killall) ]] && die "${0##*/}: could not find killall needed to kill gnuplot"
  killall $what_to_kill
  $gnuplot $opts "$(get_main_dir)/$script"
fi
