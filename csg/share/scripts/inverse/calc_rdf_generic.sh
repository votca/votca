#! /bin/bash
#
# Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

show_help () {
  cat <<EOF
${0##*/}, version %version%
This script implemtents statistical analysis for the iterative Boltzmann inversion
using generic csg tools (csg_stat)

With --only-intra-nb intramolecular non-bonded distributions are saved as
.dist-intra.new.

Usage: ${0##*/} [--help] [--only-intra-nb]
EOF
}

only_intra_nb=false
while [[ $# -gt 0 ]]
do
    key="$1"

    case "${key}" in
    --only-intra-nb)
        only_intra_nb=true
        shift  # past argument
        ;;
    --help)
        show_help
        exit 0
        ;;
    *)
        die "unknown argument ${key}"
        ;;
    esac
done

name="$(csg_get_interaction_property name)"
sim_prog="$(csg_get_property cg.inverse.program)"

topol=$(csg_get_property --allow-empty "cg.inverse.${sim_prog}.rdf.topol")
[[ -z $topol ]] && topol=$(csg_get_property "cg.inverse.${sim_prog}.topol")
[[ -f $topol ]] || die "${0##*/}: topol file '${topol}' not found, possibly you have to add it to cg.inverse.filelist"

traj=$(csg_get_property --allow-empty "cg.inverse.${sim_prog}.rdf.traj")
[[ -z $traj ]] && traj=$(csg_get_property_substitute "cg.inverse.${sim_prog}.traj")
[[ -f $traj ]] || die "${0##*/}: traj file '${traj}' not found"

maps=
#always try to find mapping files
if : ; then
  mapping=$(csg_get_property --allow-empty "cg.inverse.${sim_prog}.rdf.map")
  [[ -z $mapping ]] && mapping="$(csg_get_property --allow-empty cg.inverse.map)"
  #fail if we have bonded interaction, but no mapping file
  [[ -n $(csg_get_property --allow-empty cg.bonded.name) && -z $mapping ]] && die "Mapping file for bonded interaction needed"
  for map in ${mapping}; do
    [[ -f "$(get_main_dir)/${map}" ]] || die "${0##*/}: Mapping file '${map}' for bonded interaction not found in maindir"
    maps+="$(get_main_dir)/${map};"
  done
fi

equi_time=$(csg_get_property "cg.inverse.${sim_prog}.equi_time")
if [[ ${CSG_RUNTEST} ]] && csg_calc "${equi_time}" ">" "0"; then
  msg --color blue --to-stderr "Automatically setting equi_time to 0, because CSG_RUNTEST was set"
  equi_time=0
fi

first_frame=$(csg_get_property "cg.inverse.${sim_prog}.first_frame")
if [[ ${CSG_RUNTEST} ]] && csg_calc "${first_frame}" ">" "0"; then
  msg --color blue --to-stderr "Automatically setting first_frame to 0, because CSG_RUNTEST was set"
  first_frame=0
fi

if [[ ${only_intra_nb} = "true" ]]; then
  intra_opts="--only-intra-nb"
  dist_type="dist-intra"
  suffix="_intra"
  message=" only intramolecular correlations"
else
  intra_opts=""
  dist_type="dist"
  suffix=""
  message=""
fi

with_errors=$(csg_get_property "cg.inverse.${sim_prog}.rdf.with_errors")
if [[ ${with_errors} == "yes" ]]; then
  suffix="${suffix}_with_errors"
  block_length=$(csg_get_property "cg.inverse.${sim_prog}.rdf.block_length")
  if [[ ${CSG_RUNTEST} ]] && csg_calc "${block_length}" ">" "2"; then
    msg --color blue --to-stderr "Automatically setting block_length to 2, because CSG_RUNTEST was set"
    block_length=2
  fi
  error_opts="--block-length ${block_length}"
  ext_opt="${dist_type}.block"
else
  suffix="$suffix"
  ext_opt="${dist_type}.new"
fi

tasks=$(get_number_tasks)
#rdf calculation is maybe done already in a different interaction
if is_done "rdf_calculation${suffix}"; then
  echo "rdf calculation is already done"
else
  msg "Calculating rdfs with csg_stat using ${tasks} tasks ${message}"
  # do not put quotes around arguments with values ($error_opts)!
  # this will give a codacy warning :/
  critical csg_stat --nt "${tasks}" --options "${CSGXMLFILE}" --top "${topol}" \
    --trj "${traj}" --begin "${equi_time}" --first-frame "${first_frame}" ${error_opts} \
    "${intra_opts}" --ext "${ext_opt}" ${maps:+--cg ${maps}}
  mark_done "rdf_calculation${suffix}"
fi

# improve new RDF at low values
improve_dist_near_core_new="$(csg_get_interaction_property improve_dist_near_core.new)"
if [[ $improve_dist_near_core_new == "true" && $dist_type == "dist" ]]; then  # do not do this for dist_intra
  if ! is_done "${name}_${dist_type}_rdf_improve"; then
    fit_start_g="$(csg_get_interaction_property improve_dist_near_core.fit_start_g)"
    fit_end_g="$(csg_get_interaction_property improve_dist_near_core.fit_end_g)"
    if [[ ${with_errors} = "yes" ]]; then
      # improve all blocks
      files_to_do=${name}_*.${dist_type}.block
    else
      # improve the one file
      files_to_do=${name}.${dist_type}.new
    fi
    for i in ${files_to_do}; do
      [[ -f $i ]] || die "${0##*/}: Could not find ${i} after running csg_stat, that usually means the blocksize (cg.inverse.${sim_prog}.rdf.block_length) is too big."
      do_external dist improve_near_core --in="${i}" --out="${i}.improved" --gmin="$fit_start_g" --gmax="$fit_end_g"
      critical mv "${i}" "${i}.raw"
      critical mv "${i}.improved" "${i}"
    done
    mark_done "${name}_${dist_type}_rdf_improve"
  fi
fi

if [[ ${with_errors} = "yes" ]]; then
  if ! is_done "${name}_${dist_type}_rdf_average"; then
    for i in ${name}_*.${dist_type}.block; do
      [[ -f $i ]] || die "${0##*/}: Could not find ${name}_*.dist.block after running csg_stat, that usually means the blocksize (cg.inverse.${sim_prog}.rdf.block_length) is too big."
    done
    msg "Calculating average rdfs and its errors for interaction ${name}"
    # do not put quotes around strings that should expand to multiple files
    # this will give a codacy warning :/
    do_external table average --output "${name}.${dist_type}.new" "${name}"_*."${dist_type}".block
    mark_done "${name}_${dist_type}_rdf_average"
  fi
fi
