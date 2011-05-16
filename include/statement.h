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

#include <sqlite3.h>

namespace votca { namespace tools {

class Database;

class Statement
{
public:
	~Statement();

	template<typename T>
	void Bind(int i, const T &value);

	void Step();
	void Reset();

	sqlite3_stmt *getSQLiteStatement() { return _stmt; }
protected:
	Statement(sqlite3_stmt *stmt)
		: _stmt(stmt) {}
    sqlite3_stmt *_stmt;

	friend class Database;
};

inline void Statement::~Statement()
{
	sqlite3_finalize(_stmt);
}

template<>
inline void Statement::Bind(int i, const int &value)
{
	sqlite3_bind_int(_stmt, i, value);
}

template<>
inline void Statement::Bind(int i, const double &value)
{
	sqlite3_bind_double(_stmt, i, value);
}

template<>
inline void Statement::Bind(int i, const string &value)
{
    sqlite3_bind_text(_stmt, i, value.c_str(), -1, NULL);;
}

inline void Statement::Step()
{
	sqlite3_step(_stmt);
}

inline void Statement::Reset()
{
	sqlite3_reset(_stmt);
}

}}
