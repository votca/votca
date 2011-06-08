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
#ifndef __VOTCA_TOOLS_STATEMENT_H
#define __VOTCA_TOOLS_STATEMENT_H

#include <sqlite3.h>

namespace votca { namespace tools {

class Database;

class Statement
{
public:
	~Statement();

	template<typename T>
		void Bind(int col, const T &value);

	template<typename T>
	T Column(int col);

	int Step();
	void Reset();

	sqlite3_stmt *getSQLiteStatement() { return _stmt; }
protected:
	Statement(sqlite3_stmt *stmt)
		: _stmt(stmt) {}
    sqlite3_stmt *_stmt;

	friend class Database;
};

}}

#endif
