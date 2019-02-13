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
#ifndef __VOTCA_TOOLS_DATABASE_H
#define __VOTCA_TOOLS_DATABASE_H

#include "statement.h"
#include <sqlite3.h>
#include <string>

namespace votca {
namespace tools {

using namespace std;

/**
 *  \brief  SQLite Database wrapper
 *
 *
 */
class Database {
 public:
  Database();
  ~Database();

  sqlite3 *getSQLiteDatabase() { return _db; }

  /**
   *  \brief Helper function for opening / creating / upgrading a sqlite3
   * database
   *  @param file
   *
   *  This function helps to manage a sqlite database. It is inspired by the
   *  SQLiteOpenHelper from the Android API.
   *
   *  If this function is used to open a database, and the database
   *  doesn't exist already, onCreate is called. If it exists but in an
   *  older Version, onUpgrade is called to bring it to the most current
   * version.
   */
  void OpenHelper(string file);
  void Open(string file,
            int flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE);
  void Close(void);

  /**
   * \brief called by OpenHelper if database doesnt exists
   *
   * All necessary SQL statements to create the database if it does
   * not exist should be executed here. This is only called if the
   * database is opened via OpenHelper.
   */
  virtual void onCreate();

  /**
   * \brief called if database has an older version
   *
   * If the version specified does not match the database version,
   * onUpgrade is called to make necessary changes. This is only called if the
   * database is opened via OpenHelper.
   */
  virtual void onUpgrade(int oldVersion, int newVersion);

  void Exec(string sql);

  Statement *Prepare(string sql);

  int LastInsertRowId();
  void BeginTransaction() { Exec("BEGIN TRANSACTION;"); }
  void EndTransaction() { Exec("END TRANSACTION;"); }
  void CommitTransaction() { Exec("COMMIT TRANSACTION;"); }
  void RollbackTransaction() { Exec("ROLLBACK TRANSACTION;"); }

 protected:
  sqlite3 *_db;
};

}  // namespace tools
}  // namespace votca

#endif
