/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_CSG_XYZREADER_H
#define	__VOTCA_CSG_XYZREADER_H

#include <string>
#include <iostream>
#include <fstream>
#include <votca/csg/topologyreader.h>
#include <votca/csg/trajectoryreader.h>
#include <votca/tools/constants.h>
namespace votca { namespace csg {

/**
    \brief class for reading xyz files

    This class provides the TrajectoryReader + Topology reader interface
    for xyz files

*/
class XYZReader :
    public TrajectoryReader, public TopologyReader
{
    public:
        XYZReader() {}
        ~XYZReader() {}

        /// open a topology file
        bool ReadTopology(std::string file, Topology &top);

        /// open a trejectory file
        bool Open(const std::string &file);
        /// read in the first frame
        bool FirstFrame(Topology &top);
        /// read in the next frame
        bool NextFrame(Topology &top);


        template <class T>
        void ReadFile(T &container){
           if(!ReadFrame<true,T>(container)){
               throw std::runtime_error("Reading xyz file failed");
           }
        }

        void Close();

    private:

        template<class T>
        int getContainerSize(T& container){
            return container.size();
        }

        int getContainerSize(Topology& container){
            return container.BeadCount();
        }
        
        template <bool topology, class T>
        void AddAtom(T &container,std::string name,int id, const tools::vec& pos){
            Eigen::Vector3d pos2=pos.toEigen()*tools::conv::ang2bohr;
            container.push_back(id,name,pos2);
        }

        template <bool topology, class T>
        void AddAtom(Topology &container,std::string name,int id, const tools::vec& pos){
        Bead *b;
        tools::vec posnm=pos*tools::conv::ang2nm;
        if(topology)
            b = container.CreateBead(1, name+boost::lexical_cast<string>(id),
                    (container.GetOrCreateBeadType(name)), 0, 0, 0);
        else
            b = container.getBead(id);
            b->setPos(posnm);
        }


        template <bool topology, class T>
        bool ReadFrame(T &container);

        std::ifstream _fl;

        int _line;
};

template <bool topology, class T>
inline bool XYZReader::ReadFrame(T &container)
{
    string line;
    getline(_fl, line); ++_line;
    if(!_fl.eof()) {
        // read the number of atoms
        Tokenizer tok1(line," \t");
        std::vector<std::string> line1;
        tok1.ToVector(line1);
        if(line1.size()!=1){
          throw std::runtime_error("First line of xyz file should contain number of atoms/beads, nothing else.");
        }
        int natoms = boost::lexical_cast<int>(line1[0]);
        if(!topology && natoms !=getContainerSize(container))
            throw std::runtime_error("number of beads in topology and trajectory differ");

        // the title line
        getline(_fl, line); ++_line;

        // read atoms
        for(int i=0; i<natoms; ++i) {
            getline(_fl, line); ++_line;
            if(_fl.eof())
                throw std::runtime_error("unexpected end of file in xyz file");
            vector<string> fields;
            Tokenizer tok(line, " ");
            tok.ToVector(fields);

            if(fields.size() != 4)
                throw std::runtime_error("invalide line " +
                        boost::lexical_cast<string>(_line) +
                        " in xyz file\n" + line);


            tools::vec pos=vec(
                    boost::lexical_cast<double>(fields[1]),
                    boost::lexical_cast<double>(fields[2]),
                    boost::lexical_cast<double>(fields[3]));

            AddAtom<topology,T>(container,fields[0],i,pos);

        }
    }
    return !_fl.eof();;
}
}}

#endif

