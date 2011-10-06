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

#include "md2qm_observer.h"
#include <votca/csg/nblist.h>
#include <votca/ctp/qmnblist.h>

MD2QMObserver::MD2QMObserver()
{}


MD2QMObserver::~MD2QMObserver()
{}

void MD2QMObserver::Initialize(QMTopology &qmtop, Property &opts)
{
    _qmtop = &qmtop;
    _save.Open(qmtop, _out);
    if(_save.FramesInDatabase() > 0) {
        throw std::runtime_error("Database already contains frame information. "
                "Appending a database is not supported yet"); }
}

void MD2QMObserver::BeginCG(Topology *top, Topology *top_atom)
{
    _qmtop->Initialize(*top);
}

void MD2QMObserver::EvalConfiguration(Topology *top, Topology *top_atom)
{
    if(_qmtop->getDatabaseId() != 0) {
        throw std::runtime_error("writing several frames to state file not yet supported, please use --nframes=1");
    _qmtop->Update(*top); }
    else {
      // cout << "  reading/writing the output from/to " << _out << endl;
      _save.WriteFrame();
    }
}

void MD2QMObserver::EndCG()
{
    _save.Close();
}

