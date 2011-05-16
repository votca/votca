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

class Database
{
public:
	Database();
	~Database();

	slite3 getSQLiteDatabase() { return _db; }

	void Open(string file, int flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE)
	void Close(void);

	void Exec(string sql);

	Statement *Prepare(string sql);

protected:
	sqlite3 *_db;
};

inline Database::Database()
	: _db(NULL)
{}

inline Database::~Database()
{
	Close();
}

inline void Database::Open(string file, int flags)
{
	int ret = sqlite3_open_v2(file.c_str(),&_db,SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("cannot open database " + file);
}

inline void Database::Close()
{
    if(_db)
        sqlite3_close(_db);
}

inline void Database::Exec(string sql)
{
    char *error;
    int ret = sqlite3_exec(_db, sql.c_str(), NULL, NULL,  &error);
    if(ret != SQLITE_OK)
        throw std::runtime_error(string("cannot create frame table:\n") + error);
}

inline Statement *Database::Prepare(string sql)
{
    sqlite3_stmt *stmt;
    ret = sqlite3_prepare_v2(_db,
            sql.c_str(), -1, &stmt, NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("prepare insert frame statement failed");
    return new Statement(stmt);
}

}}
