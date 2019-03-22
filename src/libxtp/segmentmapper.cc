/*
 *            Copyright 2016 The MUSCET Development Team
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/segmentmapper.h>


namespace votca {
namespace xtp {


void SegmentMapper::ParseFragment(Seginfo& seginfo, const tools::Property& frag){
    tools::Tokenizer tok_qm_atoms(frag.get("qmatoms").as<std::string>(), " \t\n");
    std::vector<std::string> qm_atoms=tok_qm_atoms.ToVector();
    tools::Tokenizer tok_md_atoms( frag.get("mdatoms").as<std::string>(), " \t\n");
    std::vector<std::string> md_atoms=tok_md_atoms.ToVector();

    if(md_atoms.size()!=qm_atoms.size()){
        throw std::runtime_error("Mapping for segment "+seginfo.segname+
                " fragment "+frag.get("name").as<std::string>()+" does not have same numbers of md and qmatoms. "
                "If you want to leave a qmatom out, place a ':' instead");
    }

    tools::Tokenizer tok_weights(frag.get("weights").as<std::string>(), " \t\n");
    std::vector<double> weights;
    tok_weights.ConvertToVector(weights);

    if(md_atoms.size()!=weights.size()){
        throw std::runtime_error("Mapping for segment "+seginfo.segname+
                " fragment "+frag.get("name").as<std::string>()+" does not have same numbers of md and weights. "
                "If you want to leave a qmatom out, place a '0' instead");
    }

    FragInfo mapfragment;
    std::vector<int> qmatom_ids;

    for(unsigned i=0;i<qm_atoms.size();i++){
        const std::string& qm_string=qm_atoms[i];
        const std::string& md_string=md_atoms[i];
        const double& weight=weights[i];
        MD_atom_id md_result=StringToMDIndex(md_string);
        seginfo.mdatoms.push_back(md_result);

        if(qm_string==":"){continue;}
        QM_atom_id qm_result=StringToQMIndex(qm_string);
        qmatom_ids.push_back(qm_result.first);
        seginfo.qmatoms.push_back(qm_result);
        if(Atom::GetElementFromMDName(md_result.second)!=qm_result.second){
            XTP_LOG_SAVE(logINFO, _log) << "WARNING: mdatom'"<<md_result.second
                    <<"' and qmatom '"<<qm_result.second<<"' do not have same element" << std::flush;
        }
        mapfragment._qmatom_ids.push_back(qm_result);
        mapfragment._weights.push_back(weight);
        mapfragment._mdatom_ids.push_back(md_result);
    }

    tools::Tokenizer tok_frame(frag.get("localframe").as<std::string>(), " \t\n");
    std::vector<int> frame;
    tok_frame.ConvertToVector(frame);
    if(frame.size()>3){
        throw std::runtime_error("Local frame for segment "+seginfo.segname+
                " fragment "+frag.get("name").as<std::string>()+" has more than 3 atoms, please specify only up to three atoms");
    }else if (frame.empty()){
        throw std::runtime_error("No local frame for segment "+seginfo.segname+
                " fragment "+frag.get("name").as<std::string>()+" specified");
    }
    for(int atomid:frame){
        if(std::find (qmatom_ids.begin(), qmatom_ids.end(), atomid)==qmatom_ids.end()){
            throw std::runtime_error("Atom "+std::to_string(atomid)+" in local frame cannot be found in qmatoms.");
        }
    }

    mapfragment._qm_local_frame=frame;
    seginfo.fragments.push_back(mapfragment);
}

void SegmentMapper::LoadMappingFile(const std::string& mapfile){

  tools::Property topology_map;
  tools::load_property_from_xml(topology_map, mapfile);

  std::string molkey = "topology.molecules.molecule";
  std::vector<tools::Property*> molecules = topology_map.Select(molkey);
  std::string segkey = "segments.segment";
  for (tools::Property* mol : molecules) {
    std::vector<tools::Property*> segments = mol->Select(segkey);
    for (tools::Property* seg : segments) {
        Seginfo seginfo;

        std::string coordfile_key="qmcoords_*";
        std::vector<tools::Property*> files = seg->Select(coordfile_key);
        for(tools::Property* file:files){
            seginfo.coordfiles[file->name()]=file->as<std::string>();
        }
        std::string segname = seg->get("name").as<std::string>();
        seginfo.segname=segname;
        std::string fragkey = "fragments.fragment";

        seginfo.map2md=seg->ifExistsReturnElseThrowRuntimeError<bool>("map2md");

        std::vector<tools::Property*> fragments = seg->Select(fragkey);
        for (tools::Property* frag : fragments) {
            ParseFragment(seginfo, *frag);
        }

        int qm_atom_min_id = std::min_element(seginfo.qmatoms.begin(), seginfo.qmatoms.end(),[](const QM_atom_id& a, const QM_atom_id & b) {
            return a.first < b.first;
        })->first;            
        if(qm_atom_min_id!=0){
            throw std::runtime_error("QM atoms for segment "+seginfo.segname+" do not start at zero index. Each segment should have its own qm coordinate file");
        }

        seginfo.minmax=CalcResidueRange(seginfo.mdatoms);
        _segment_info[segname]=seginfo;

        }
      }
    }


  std::pair<int,std::string> SegmentMapper::StringToQMIndex(const std::string& qm_string)const{
      tools::Tokenizer tok(qm_string,":");
      std::vector<std::string> result=tok.ToVector();
      if(result.size()!=2){
        throw std::runtime_error("Entry "+qm_string+" is not properly formatted.");
      }
      return std::pair<int,std::string>(std::stoi(result[0]),result[1]);
  }

  MD_atom_id SegmentMapper::StringToMDIndex(const std::string& md_string)const{
      tools::Tokenizer tok(md_string,":");
      std::vector<std::string> result=tok.ToVector();
      if(result.size()!=3){
        throw std::runtime_error("Entry "+md_string+" is not properly formatted.");
      }
      return MD_atom_id(std::stoi(result[0]),result[2]);
  }



std::pair<int,int> SegmentMapper::CalcResidueRange(const std::vector<MD_atom_id>& seg)const{
 int max_res_id=std::min_element(seg.begin(),seg.end(),[](const MD_atom_id& a,const MD_atom_id& b){
     return a.first<b.first;})->first;
 int min_res_id=std::min_element(seg.begin(),seg.end(),[](const MD_atom_id& a,const MD_atom_id& b){
     return a.first<b.first;})->first;
 return std::pair<int,int>(min_res_id,max_res_id);
}

std::pair<int,int> SegmentMapper::CalcResidueRange(const Segment& seg)const{
    int max_res_id=std::min_element(seg.begin(),seg.end(),[](const Atom& a,const Atom&b){
        return a.getResnr()<b.getResnr();
    })->getResnr();

    int min_res_id=std::min_element(seg.begin(),seg.end(),[](const Atom& a,const Atom&b){
        return a.getResnr()<b.getResnr();
    })->getResnr();
    return std::pair<int,int>(min_res_id,max_res_id);
}

template<class T>
     Eigen::Vector3d SegmentMapper::CalcWeightedPos(const FragInfo& frag,const T& atoms)const{
    Eigen::Vector3d qm_pos=Eigen::Vector3d::Zero();
    double qm_weight=0.0;
    for(unsigned i=0;i<atoms.size();i++){
        qm_pos+=atoms[i]->getPos()*frag._weights[i];
        qm_weight+=frag._weights[i];
    }
    return qm_pos=qm_pos/qm_weight;
}

void SegmentMapper::PlaceQMonMD(const std::vector<QMAtom*>& fragment_qmatoms,
                                const std::vector<const Atom*>& fragment_mdatoms) const{
    for (unsigned i = 0; i < fragment_qmatoms.size(); i++) {
        const Atom* a = fragment_mdatoms[i];
        QMAtom* b = fragment_qmatoms[i];
        b->setPos(a->getPos());
    }
}


int SegmentMapper::FindVectorIndexFromAtomId(int atomid,
                                             const std::vector<QMAtom*>& fragment_qmatoms)const{
    unsigned i=0;
        for(;i<fragment_qmatoms.size();i++){
            if(fragment_qmatoms[i]->getId()==atomid){
                break;
            }
        }
    return i;
}

void SegmentMapper::MapQMonMD(const FragInfo& frag, const std::vector<QMAtom*>& fragment_qmatoms,
                              const std::vector<const Atom*>& fragment_mdatoms)const{
    std::vector<Eigen::Vector3d> local_qm_frame;
    std::vector<Eigen::Vector3d> local_md_frame;
    for(int id:frag._qm_local_frame){
        int i=FindVectorIndexFromAtomId(id,fragment_qmatoms);
        local_qm_frame.push_back(fragment_qmatoms[i]->getPos());
        local_md_frame.push_back(fragment_mdatoms[i]->getPos());
    }

    int symmetry=frag._qm_local_frame.size();
    Eigen::Vector3d qm_com=CalcWeightedPos(frag, fragment_qmatoms);
    Eigen::Vector3d md_com=CalcWeightedPos(frag, fragment_mdatoms);

    Eigen::Vector3d shift_qm2md=md_com-qm_com;

    Eigen::Matrix3d rot_qm=Eigen::Matrix3d::Identity();
    Eigen::Matrix3d rot_md=Eigen::Matrix3d::Identity();

    //building local frame
    if(symmetry>1){
        //middle atom is the origin
        Eigen::Vector3d x_qm=local_qm_frame[0]-local_qm_frame[1];
        Eigen::Vector3d x_md=local_md_frame[0]-local_md_frame[1];
        Eigen::Vector3d y_qm;
        Eigen::Vector3d y_md;
        Eigen::Vector3d z_qm;
        Eigen::Vector3d z_md;
        if(symmetry==3){
            y_qm=local_qm_frame[2]-local_qm_frame[1];
            y_md=local_md_frame[2]-local_md_frame[1];
            z_qm=x_qm.cross(y_qm);
            z_md=x_md.cross(y_md);
            y_qm=z_qm.cross(x_qm);
            y_md=z_md.cross(x_md);

        }else{
            Eigen::Vector3d unit=Eigen::Vector3d::UnitX();
            if(std::abs(unit.dot(x_qm)/x_qm.norm())<1e-6){
                unit=Eigen::Vector3d::UnitY();
            }
            y_qm=(x_qm.cross(unit));
            if(std::abs(unit.dot(x_md)/x_md.norm())<1e-6){
                unit=Eigen::Vector3d::UnitX();
            }
            y_md=(x_md.cross(unit));
            z_qm=x_qm.cross(y_qm);
            z_md=x_md.cross(y_md);
        }

        if(x_qm.squaredNorm()<1e-18 || y_qm.squaredNorm()<1e-18 || z_qm.squaredNorm()<1e-18 ){
            throw std::runtime_error("QM basis vectors are very small, choose different local basis");
        }
        rot_qm.col(0)=x_qm.normalized();
        rot_qm.col(1)=y_qm.normalized();
        rot_qm.col(2)=z_qm.normalized();

        
        if(x_md.squaredNorm()<1e-18 || y_md.squaredNorm()<1e-18 || z_md.squaredNorm()<1e-18 ){
            throw std::runtime_error("MD basis vectors are very small, choose different local basis");
        }
        rot_md.col(0)=x_md.normalized();
        rot_md.col(1)=y_md.normalized();
        rot_md.col(2)=z_md.normalized();
    }
    Eigen::Matrix3d rotateQM2MD=rot_md*rot_qm.transpose();
    for(QMAtom* atom:fragment_qmatoms){
        atom->Translate(shift_qm2md);
        atom->Rotate(rotateQM2MD,md_com);
    }

}

QMMolecule SegmentMapper::map(const Segment& seg,QMState state)const{
   Seginfo seginfo=_segment_info.at(seg.getName());

   if(int(seginfo.mdatoms.size())!=seg.size()){
        throw std::runtime_error("Segment '"+seg.getName()+
                "' does not contain the same number of atoms as mapping file: "+
                std::to_string(seginfo.mdatoms.size())+" vs. "+std::to_string(seg.size()));
    }

   std::pair<int,int> minmax_map=seginfo.minmax;
   std::pair<int,int> minmax=CalcResidueRange(seg);

    if((minmax_map.first-minmax_map.second)!=(minmax.first-minmax.second)){
        throw std::runtime_error("Residue range for segment "+seg.getName()+":"+std::to_string(seg.getId())+
                " and the mapping do not agree: Segment["+std::to_string(minmax.first)+","
                +std::to_string(minmax.second)+"] Map["+std::to_string(minmax_map.first)+","
                +std::to_string(minmax_map.second)+"]");
    }

   int residueoffset=minmax.first-minmax_map.first;
   if (seginfo.coordfiles.count("qmcoords_"+state.ToString()) == 0){
    throw std::runtime_error("Could not find a coordinate file for segment/state: qmcoords_"+state.ToString());
   }
   std::string template_coordinates_filename=
           seginfo.coordfiles.at("qmcoords_"+state.ToString());

   QMMolecule Result(seg.getName(),seg.getId());
   Result.LoadFromXYZ(template_coordinates_filename);

   if(int(seginfo.qmatoms.size())!=Result.size()){
        throw std::runtime_error("QMSegment '"+seg.getName()+
                "' does not contain the same number of atoms as mapping file: "+
                std::to_string(seginfo.qmatoms.size())+" vs. "+std::to_string(Result.size()));
    }

   for(FragInfo& frag:seginfo.fragments){
       for (MD_atom_id& id:frag._mdatom_ids){
           id.first+=residueoffset;
       }
       

       std::vector<QMAtom*> fragment_qmatoms;
       for(const QM_atom_id& id:frag._qmatom_ids){
           if(id.second!=Result[id.first].getElement()){
               throw std::runtime_error("Element of mapping atom "+std::to_string(id.first)
                       +":"+id.second+" does not agree with Element of parsed Atom"+Result[id.first].getElement());
           }
           fragment_qmatoms.push_back(&Result[id.first]);
       }
       std::vector<const Atom*> fragment_mdatoms;
       for(const MD_atom_id& id:frag._mdatom_ids){
           const Atom* atom=seg.getAtom(id);
           if(atom==nullptr){
               throw std::runtime_error("Could not find an atom with name "+
                       std::to_string(id.first)+":"+id.second
                       +" in segment "+seg.getName());
           }
           fragment_mdatoms.push_back(atom);
       }

       if(seginfo.map2md){
           PlaceQMonMD(fragment_qmatoms, fragment_mdatoms);
       }else{
        MapQMonMD(frag, fragment_qmatoms, fragment_mdatoms);
       }
       
       
   }

   Result.getPos();
   return Result;
}
  

}  // namespace xtp
}  // namespace votca
