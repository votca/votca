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
#include <stdexcept>
#include <votca/tools/database.h>

namespace votca {
namespace tools {

Database::Database() : _db(NULL) {}

Database::~Database() { Close(); }

void Database::Open(string file, int flags) {
  int ret = sqlite3_open_v2(file.c_str(), &_db, flags, NULL);
  if (ret != SQLITE_OK)
    throw std::runtime_error("cannot open database " + file);
}

void Database::OpenHelper(string file) {
  int ret = sqlite3_open_v2(file.c_str(), &_db, SQLITE_OPEN_READWRITE, NULL);
  if (ret != SQLITE_OK) {
    Open(file, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE);
    onCreate();
  }
  // TODO: onUpgrade
}

void Database::onCreate() {
  throw std::runtime_error(
      "database is opened via OpenHelper but onCreate is not imlemented");
}

void Database::onUpgrade(int oldVersion, int newVersion) {
  throw std::runtime_error(
      "database is opened via OpenHelper but onUpgrade is not implemented");
}

void Database::Close() {
  if (_db) sqlite3_close(_db);
}

void Database::Exec(string sql) {
  char *error;
  int ret = sqlite3_exec(_db, sql.c_str(), NULL, NULL, &error);
  if (ret != SQLITE_OK)
    throw std::runtime_error(string("sql execute failed\n") + error +
                             "\nSQL: " + sql);
}

Statement *Database::Prepare(string sql) {
  // char *error;
  sqlite3_stmt *stmt;
  int ret = sqlite3_prepare_v2(_db, sql.c_str(), -1, &stmt, NULL);
  if (ret != SQLITE_OK)
    throw std::runtime_error(string("prepare statement failed") +
                             "\nSQL: " + sql);
  return new Statement(stmt);
}

int Database::LastInsertRowId() { return sqlite3_last_insert_rowid(_db); }

}  // namespace tools
}  // namespace votca
