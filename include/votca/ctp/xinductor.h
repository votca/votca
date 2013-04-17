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

#ifndef __XINDUCTOR__H
#define	__XINDUCTOR__H

#include <votca/ctp/topology.h>
#include <votca/ctp/apolarsite.h>
#include <votca/ctp/xinteractor.h>
#include <votca/ctp/xjob.h>
#include <votca/tools/thread.h>
#include <votca/tools/mutex.h>


// TODO Pass sharpness (aDamp) to XInteractor
// TODO Only one wSOR required, outsource to XJob
// TODO No need for _top
// TODO No need for segment containers?

namespace votca { namespace ctp {
    
    
class XInductor
{
public:
    
    // Use via XInductor(...) plus XInductor::Evaluate(...)
    
    XInductor()
            : _induce(0),          _induce_intra_pair(0),
              _wSOR_N(0.5),        _wSOR_C(0.5),
              _epsTol(0.001),      _maxIter(512),
              _maverick(true),     _top(NULL), _id(-1)
            { _actor = XInteractor(NULL); };
    
    XInductor(bool induce,   bool induce_intra_pair, int subthreads,
              float wSOR,    double epsTol,          int maxIter,
              bool maverick, Topology *top,          int thread_id)
            : _induce(induce),         _induce_intra_pair(induce_intra_pair),
              _subthreads(subthreads), _wSOR_N(wSOR),       _wSOR_C(wSOR),   
              _epsTol(epsTol),         _maxIter(maxIter),
              _maverick(maverick),     _top(top), _id(thread_id) 
            { _actor = XInteractor(top); };
            
   ~XInductor() {};
   
   
    // +++++++++++++++++++++++++++ //
    // Induction Op    (Subthread) //
    // +++++++++++++++++++++++++++ //
   
    class InduWorker : public Thread
    {
    public:

        InduWorker(int id, Topology *top, XInductor *forker)
                    : _id(id), _top(top), _forker(forker)
                  { _actor = XInteractor(top); };

       ~InduWorker() {};

        void Run(void) {

            if (_switch_energy_induce == 1) {
                while (_forker->NextChunkTodo(_id,1)) {
                    this->InterInduce();
                    _forker->OpenChunks(_chunk1, _chunk2);
                }
            }
            else if (_switch_energy_induce == 0) {
                while (_forker->NextChunkTodo(_id,0)) {
                    this->InterEnergy();
                    _forker->OpenChunks(_chunk1, _chunk2);
                }
            }
            else {
                assert(false);
            }


        }

        void InitSpheres(vector< Segment* > *vsegs_cut1,
                         vector< Segment* > *vsegs_cut2,
                         vector< vector<APolarSite*> > *vvpoles_cut1,
                         vector< vector<APolarSite*> > *vvpoles_cut2) {
            _vsegs_cut1 = vsegs_cut1;
            _vsegs_cut2 = vsegs_cut2;
            _vvpoles_cut1 = vvpoles_cut1;
            _vvpoles_cut2 = vvpoles_cut2;

            _vvpoles = vvpoles_cut1;
        }

        void SetSwitch(int energy_induce) {
            _switch_energy_induce = energy_induce;
            if (energy_induce == 0) {                
                _actor.ResetEnergy();
                
                _E_Pair_Pair = 0.0;
                _E_Pair_Sph1 = 0.0;
                _E_Sph1_Sph1 = 0.0;
                
                _E_f_C_non_C      = 0.0;
                _E_f_non_C_non_C  = 0.0;
                _E_f_C_C          = 0.0;
                _E_m_C            = 0.0;
                _E_m_non_C        = 0.0;
            }            
        }

        void Setc1c2(int c1, int c2) {
            _chunk1 = c1;
            _chunk2 = c2;
        }

        void Setnx12ny12(int nx1, int nx2, int ny1, int ny2) {
            _nx1 = nx1;
            _nx2 = nx2;
            _ny1 = ny1;
            _ny2 = ny2;
        }

        void PrintInfo() {
            printf("%1d : nx1 %2d nx2 %2d ny1 %2d ny2 %2d\n",
                    _id, _nx1, _nx2, _ny1, _ny2);
        }

        void InterInduce() {

            for (int i = _nx1;                     i < _nx2; ++i) {
            for (int j = (i >= _ny1) ? i+1 : _ny1; j < _ny2; ++j) {

                for (pit1 = (*_vvpoles)[i].begin(); pit1 < (*_vvpoles)[i].end(); ++pit1) {
                for (pit2 = (*_vvpoles)[j].begin(); pit2 < (*_vvpoles)[j].end(); ++pit2) {

                    _actor.FieldIndu(*(*pit1), *(*pit2));
                }}
            }}
        }

        double       GetEPairPair() { return _E_Pair_Pair; }
        double       GetEPairSph1() { return _E_Pair_Sph1; }
        double       GetESph1Sph1() { return _E_Sph1_Sph1; }
        
        double       GetE_f_C_non_C()           { return _E_f_C_non_C; }
        double       GetE_f_non_C_non_C()       { return _E_f_non_C_non_C; }
        double       GetE_f_C_C()               { return _E_f_C_C; }
        double       GetE_m_C()                 { return _E_m_C; }
        double       GetE_m_non_C()             { return _E_m_non_C; }
        
        XInteractor &GetActor()     { return _actor; }

        void InterEnergy() {
            
            double e_f_12_21 = 0.0;
            double e_m_12    = 0.0;
            double e_m_21    = 0.0;
            
            if (_forker->_job->getType() == "site") {
                
                for (int i = _nx1;                     i < _nx2; ++i) {
                for (int j = (i >= _ny1) ? i+1 : _ny1; j < _ny2; ++j) {

                    // Site-non-site interaction
                    if (this->_forker->_job->getSiteId() == (*_vsegs_cut1)[i]->getId()
                     || this->_forker->_job->getSiteId() == (*_vsegs_cut1)[j]->getId()) {

                        for (pit1 = (*_vvpoles_cut1)[i].begin();
                             pit1 < (*_vvpoles_cut1)[i].end();
                             ++pit1) {
                        for (pit2 = (*_vvpoles_cut1)[j].begin();
                             pit2 < (*_vvpoles_cut1)[j].end();
                             ++pit2) {

                            e_f_12_21        = _actor.E_f(*(*pit1),*(*pit2));
                            e_m_12           = _actor.E_m(*(*pit1),*(*pit2));
                            e_m_21           = _actor.E_m(*(*pit2),*(*pit1));

                            _E_Pair_Sph1     += e_f_12_21 + e_m_12 + e_m_21;

                            _E_f_C_non_C += e_f_12_21;
                            if (this->_forker->_job->getSiteId() == (*_vsegs_cut1)[i]->getId()) {
                                _E_m_C += e_m_12;
                                _E_m_non_C += e_m_21;
                            }
                            else {
                                _E_m_C += e_m_21;
                                _E_m_non_C += e_m_12;
                            }
                        }}
                    }

                    // Non-site-non-site interaction
                    else {
                        for (pit1 = (*_vvpoles_cut1)[i].begin();
                             pit1 < (*_vvpoles_cut1)[i].end();
                             ++pit1) {
                        for (pit2 = (*_vvpoles_cut1)[j].begin();
                             pit2 < (*_vvpoles_cut1)[j].end();
                             ++pit2) {

                            e_f_12_21        = _actor.E_f(*(*pit1),*(*pit2));
                            e_m_12           = _actor.E_m(*(*pit1),*(*pit2));
                            e_m_21           = _actor.E_m(*(*pit2),*(*pit1));

                            _E_Sph1_Sph1     += e_f_12_21 + e_m_12 + e_m_21;

                            _E_f_non_C_non_C += e_f_12_21;
                            _E_m_non_C       += e_m_12;
                            _E_m_non_C       += e_m_21;
                        }}
                    }
                }}
            } // end procedure - job type "site"
            
            else if (_forker->_job->getType() == "pair") {
                
                for (int i = _nx1;                     i < _nx2; ++i) {
                for (int j = (i >= _ny1) ? i+1 : _ny1; j < _ny2; ++j) {

                    // Pair-non-pair interaction
                    if ( (this->_forker->_job->getSeg1Id() == (*_vsegs_cut1)[i]->getId()
                       || this->_forker->_job->getSeg2Id() == (*_vsegs_cut1)[i]->getId())
                     ^   (this->_forker->_job->getSeg1Id() == (*_vsegs_cut1)[j]->getId()
                       || this->_forker->_job->getSeg2Id() == (*_vsegs_cut1)[j]->getId())) {

                        for (pit1 = (*_vvpoles_cut1)[i].begin();
                             pit1 < (*_vvpoles_cut1)[i].end();
                             ++pit1) {
                        for (pit2 = (*_vvpoles_cut1)[j].begin();
                             pit2 < (*_vvpoles_cut1)[j].end();
                             ++pit2) {

                            e_f_12_21        = _actor.E_f(*(*pit1),*(*pit2));
                            e_m_12           = _actor.E_m(*(*pit1),*(*pit2));
                            e_m_21           = _actor.E_m(*(*pit2),*(*pit1));

                            _E_Pair_Sph1     += e_f_12_21 + e_m_12 + e_m_21;

                            _E_f_C_non_C += e_f_12_21;
                            if (  this->_forker->_job->getSeg1Id() == (*_vsegs_cut1)[i]->getId()
                               || this->_forker->_job->getSeg2Id() == (*_vsegs_cut1)[i]->getId() ) {
                                _E_m_C += e_m_12;
                                _E_m_non_C += e_m_21;
                            }
                            else {
                                _E_m_C += e_m_21;
                                _E_m_non_C += e_m_12;
                            }
                        }}
                    }
                    
                    // Pair-pair interaction
                    else if ( (this->_forker->_job->getSeg1Id() == (*_vsegs_cut1)[i]->getId()
                            || this->_forker->_job->getSeg2Id() == (*_vsegs_cut1)[i]->getId())
                         &&   (this->_forker->_job->getSeg1Id() == (*_vsegs_cut1)[j]->getId()
                            || this->_forker->_job->getSeg2Id() == (*_vsegs_cut1)[j]->getId())) {
                        for (pit1 = (*_vvpoles_cut1)[i].begin();
                             pit1 < (*_vvpoles_cut1)[i].end();
                             ++pit1) {
                        for (pit2 = (*_vvpoles_cut1)[j].begin();
                             pit2 < (*_vvpoles_cut1)[j].end();
                             ++pit2) {
                            
                            e_f_12_21        = _actor.E_f(*(*pit1),*(*pit2));
                            e_m_12           = _actor.E_m(*(*pit1),*(*pit2));
                            e_m_21           = _actor.E_m(*(*pit2),*(*pit1));
                            
                            _E_Pair_Pair    += e_f_12_21 + e_m_12 + e_m_21;
                            
                            _E_f_C_C        += e_f_12_21;
                            _E_m_C          += e_m_12;
                            _E_m_C          += e_m_21;
                            
                        }}
                    }

                    // Non-pair-non-pair interaction
                    else {
                        for (pit1 = (*_vvpoles_cut1)[i].begin();
                             pit1 < (*_vvpoles_cut1)[i].end();
                             ++pit1) {
                        for (pit2 = (*_vvpoles_cut1)[j].begin();
                             pit2 < (*_vvpoles_cut1)[j].end();
                             ++pit2) {

                            e_f_12_21        = _actor.E_f(*(*pit1),*(*pit2));
                            e_m_12           = _actor.E_m(*(*pit1),*(*pit2));
                            e_m_21           = _actor.E_m(*(*pit2),*(*pit1));

                            _E_Sph1_Sph1     += e_f_12_21 + e_m_12 + e_m_21;

                            _E_f_non_C_non_C += e_f_12_21;
                            _E_m_non_C       += e_m_12;
                            _E_m_non_C       += e_m_21;
                        }}
                    }
                }}         
            } // end procedure - job type "pair"
            
            else { assert(false); }
            
            
        }
        

      private:

          int                                   _id;
          Topology                             *_top;
          XJob                                 *_job;
          XInductor                            *_forker;
          XInteractor                           _actor;
          vector< vector<APolarSite*> >        *_vvpoles;

          vector< Segment* >                   *_vsegs_cut1;
          vector< Segment* >                   *_vsegs_cut2;
          vector< vector<APolarSite*> >        *_vvpoles_cut1;
          vector< vector<APolarSite*> >        *_vvpoles_cut2;

          int _nx1, _nx2;
          int _ny1, _ny2;
          int _chunk1, _chunk2;

          int _switch_energy_induce;

          vector< Segment* >               ::iterator      seg1;
          vector< Segment* >               ::iterator      seg2;
          vector< vector<APolarSite*> >    ::iterator      sit1;
          vector< vector<APolarSite*> >    ::iterator      sit2;
          vector< APolarSite* >            ::iterator      pit1;
          vector< APolarSite* >            ::iterator      pit2;

          double _E_Pair_Pair;
          double _E_Pair_Sph1;
          double _E_Sph1_Sph1;
          
          double _E_f_C_non_C;          // interaction central <> Sph1, Sph2
          double _E_f_non_C_non_C;      // interaction Sph1    <> Sph1
          double _E_f_C_C;              // interaction central <> central
          double _E_m_C;                // induction work central
          double _E_m_non_C;            // induction work not central
    };
    
    
    // +++++++++++++++++++++++++++ //
    // Manage Subthreads           //
    // +++++++++++++++++++++++++++ //
    
    void        ClearTodoTable() {
        for (int i = 0; i < _xy_done.size(); ++i) {
        for (int j = 0; j < _xy_done[i].size(); ++j) {
                _xy_done[i][j] = false;
        }}
    }

    void        InitChunks() {
        
        _nx1.clear();
        _nx2.clear();
        _ny1.clear();
        _ny2.clear();
        
        _xy_done.clear();
        _chunks_avail.clear();

        int T = this->_subthreads;                  // Threads
        int C = T * 2;                              // Chunks
        int N = _polsPolSphere.size();              // Elements
        int nr = N % C;                             // Rest size
        int nt = (N-nr) / C;                        // Chunk size
        
        assert (N == C*nt + nr);

        for (int id = 0; id < C+1; ++id) {
            _nx1.push_back( vector<int>(C+1,0) );
            _nx2.push_back( vector<int>(C+1,0) );
            _ny1.push_back( vector<int>(C+1,0) );
            _ny2.push_back( vector<int>(C+1,0) );

            _xy_done.push_back( vector<bool>(C+1,false) );
            _chunks_avail.push_back(true);
        }

        for (int col = 0; col < C+1; ++col) {
            for (int row = 0; row < C+1; ++row) {

                if (col < row) {
                    _nx1[row][col] = 0; _ny1[row][col] = 0;
                    _nx2[row][col] = 0; _ny2[row][col] = 0;
                }
                else if (col == C && row == C) {
                    _nx1[row][col] = C*nt;
                    _nx2[row][col] = C*nt + nr;
                    _ny1[row][col] = C*nt;
                    _ny2[row][col] = C*nt + nr;
                }
                else if (col == C && row < C) {
                    _nx1[row][col] = row*nt;
                    _nx2[row][col] = (row+1)*nt;
                    _ny1[row][col] = C*nt;
                    _ny2[row][col] = C*nt + nr;
                }
                else {
                    _nx1[row][col] = row*nt;
                    _nx2[row][col] = (row+1)*nt;
                    _ny1[row][col] = col*nt;
                    _ny2[row][col] = (col+1)*nt;
                }
            }
        }

        if (T > 1 && _maverick) {
            printf("\n\nTHREAD %1d MESH LOAD: "
                   "NST%1d C%1d N%1d nt%1d nr%1d\n", _id, T, C, N, nt, nr);
            for (int id = 0; id < C+1; ++id) {
                for (int run = 0; run < C+1; ++run) {
                    printf("--------+");
                }
                printf("\n");
                for (int run = 0; run < C+1; ++run) {
                    printf("%3d %3d |", _nx1[id][run], _ny1[id][run]);
                }
                printf("\n");
                for (int run = 0; run < C+1; ++run) {
                    printf("%3d %3d |", _nx2[id][run], _ny2[id][run]);
                }
                printf("\n");
            }
            for (int run = 0; run < C+1; ++run) {
                printf("--------+");
            }
            printf("\n");
        }
    }

    void        OpenChunks(int c1, int c2) {
        _chunks_avail[c1] = true;
        _chunks_avail[c2] = true;            
    }

    bool        NextChunkTodo(int indu_id, int switch_energy_induce) {

        _alloc_chunk.Lock();

        bool todo = false;

        while (true) {

            for (int i = 0; i < _xy_done.size(); ++i) {
            for (int j = i; j < _xy_done[i].size(); ++j) {
                if (!_xy_done[i][j]) {
                    todo = true;
                    if (_chunks_avail[i] && _chunks_avail[j]) {

                        if (switch_energy_induce == 1) {
                            _chunks_avail[i] = false;
                            _chunks_avail[j] = false;
                        }

                        _xy_done[i][j] = true;
                        _indus[indu_id]->Setc1c2(i,j);
                        _indus[indu_id]->Setnx12ny12(_nx1[i][j],
                                                     _nx2[i][j],
                                                     _ny1[i][j],
                                                     _ny2[i][j]);

                        _alloc_chunk.Unlock();
                        return todo;
                    }
                    else { ; }
                }
            }}

            if (!todo) { break; }

        }

        _alloc_chunk.Unlock();
        return todo;
    }
    
    
    // +++++++++++++++++++++++++++ //
    // Energy Computation          //
    // +++++++++++++++++++++++++++ //
    
    void        Evaluate(XJob *job, 
                         vector< Segment* >            &segsPolSphere,
                         vector< Segment* >            &segsOutSphere,
                         vector< vector<APolarSite*> > &polsPolSphere,
                         vector< vector<APolarSite*> > &polsOutSphere,
                         vector< vector<APolarSite*> > &polarSites,
                         vector< vector<APolarSite*> > &polarSites_job);
    int         Induce(int state, XJob *job);
    double      Energy(int state, XJob *job);
    double      EnergyStatic(int state, XJob *job);
    
private:    
    
    Topology                     *_top;
    XJob                         *_job;
    int                           _id; // Thread ID
    
    // Control options
    bool                          _induce;
    bool                          _induce_intra_pair;
    bool                          _maverick;
    int                           _subthreads;
    
    // Interaction parameters
    double                        _aDamp;
    
    // Convergence parameters
    float                         _wSOR_N;
    float                         _wSOR_C;
    double                        _epsTol;
    int                           _maxIter;
    
    // Polar-site containers
    vector< Segment* >            _segsPolSphere;  // Segments    in c/o 0-1
    vector< Segment* >            _segsOutSphere;  // Segments    in c/0 1-2
    vector< vector<APolarSite*> > _polsPolSphere;  // Polar sites in c/o 0-1
    vector< vector<APolarSite*> > _polsOutSphere;  // Polar sites in c/o 1-2
    vector< vector<APolarSite*> > _polarSites;     // Copy of top polar sites
    vector< vector<APolarSite*> > _polarSites_job; // Adapted to job specs 
    XInteractor                   _actor;

    // Manage induction workers
    vector< InduWorker* >         _indus;
    Mutex                         _alloc_chunk;
    vector< bool >                _chunks_avail;
    vector< vector<bool> >        _xy_done;
    vector< vector<int> >         _nx1;
    vector< vector<int> >         _nx2;
    vector< vector<int> >         _ny1;
    vector< vector<int> >         _ny2;

};  
    



}}


#endif