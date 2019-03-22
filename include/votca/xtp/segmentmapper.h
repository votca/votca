/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef VOTCA_XTP_SEGMENTMAPPER_H
#define VOTCA_XTP_SEGMENTMAPPER_H

#include <votca/xtp/topology.h>
#include <votca/xtp/qmmolecule.h>
#include <votca/tools/property.h>
#include <votca/xtp/logger.h>
#include <votca/csg/pdbwriter.h>

namespace votca
{
namespace xtp
{

class SegmentMapper
{
public:

    SegmentMapper(Logger& log): _log(log){};

    void LoadMappingFile(const std::string& mapfile);

    QMMolecule map(const Segment& seg, QMState state)const;

private:

    struct FragInfo{
        std::vector<double> _weights;
        std::vector<QM_atom_id> _qmatom_ids;
        std::vector<MD_atom_id> _mdatom_ids;
        std::vector<int> _qm_local_frame;
    };

    struct Seginfo{
        std::pair<int, int> minmax;
        std::vector< MD_atom_id > mdatoms;
        std::vector< FragInfo> fragments;
        bool map2md;
        std::string segname;
        std::vector<double> weights;
        std::vector< QM_atom_id > qmatoms;
        std::map<std::string,std::string> coordfiles;
    };

    int FindVectorIndexFromAtomId(int atomid,
				const std::vector<QMAtom*>& fragment_qmatoms)const;

    void ParseFragment(Seginfo& seginfo, const tools::Property& frag);

    template<class T>
    Eigen::Vector3d CalcWeightedPos(const FragInfo& frag,
                                    const T& atoms)const;

    void PlaceQMonMD(const std::vector<QMAtom*>& fragment_qmatoms,
                     const std::vector<const Atom*>& fragment_mdatoms)const;

    void MapQMonMD(const FragInfo& frag, const std::vector<QMAtom*>& fragment_qmatoms,
                   const std::vector<const Atom*>& fragment_mdatoms)const;

    Logger &_log;
    std::pair<int, int> CalcResidueRange(const Segment& seg)const;
    std::pair<int, int> CalcResidueRange(const std::vector<MD_atom_id>& seg)const;

    QM_atom_id StringToQMIndex(const std::string& qm_string)const;

    MD_atom_id StringToMDIndex(const std::string& md_string)const;

    std::map<std::string, Seginfo > _segment_info;

};

} // namespace xtp
} // namespace votca

#endif  // VOTCA_XTP_SEGMENTMAPPER_H
