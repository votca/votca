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
#ifndef __VOTCA_TOOLS_DATABASE_H
#define __VOTCA_TOOLS_DATABASE_H

#include <string>
#include <sqlite3.h>
#include "statement.h"

namespace votca { namespace tools {

using namespace std;

class Database
{
public:
	Database();
	~Database();

	sqlite3 *getSQLiteDatabase() { return _db; }

	void Open(string file, int flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE);
	void Close(void);

	void Exec(string sql);

	Statement *Prepare(string sql);

protected:
	sqlite3 *_db;
};

}}

#endif
