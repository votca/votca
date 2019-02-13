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
#ifndef __VOTCA_TOOLS_STATEMENT_H
#define __VOTCA_TOOLS_STATEMENT_H

#include "lexical_cast.h"
#include <sqlite3.h>
#include <stdexcept>
#include <string>

namespace votca {
namespace tools {

class Database;

/**
 * \brief Wrapper for sqlite prepared statements
 *
 * The Statement class wraps the interface to sqlite_stmt. It checks
 * for some of the basic errors and throws an exception in case
 * one occurs.
 */
class Statement {
 public:
  ~Statement();

  /**
   * \brief bind a value to prepared statement
   * @param col column number, sqlite starts counting with 1
   * @param value value
   */
  template <typename T>
  void Bind(int col, const T &value);

  /**
   * \brief read a column after a select statement was executed
   * @param col column number, sqlite starts counting with 0 here
   */
  template <typename T>
  T Column(int col);

  /**
   * \brief perform a step
   * @return sqlite return code, see sqlite manual for details
   */
  int Step();

  /**
   * \brief perform an insert step
   *
   * This is basically just a call to Step and does an additional
   * error check if SQLITE_DONE was returned (= insert successful).
   * If not an exception is thrown.
   */
  int InsertStep();

  /**
   * \brief reset the statment to perform another insert or query
   */
  void Reset();

  sqlite3_stmt *getSQLiteStatement() { return _stmt; }

 protected:
  Statement(sqlite3_stmt *stmt) : _stmt(stmt) {}
  sqlite3_stmt *_stmt;

  friend class Database;
};

inline int Statement::InsertStep() {
  int ret = Step();
  if (ret != SQLITE_DONE)
    throw std::runtime_error(
        "Statment::Step did not return SQLITE_DONE. Return code was " +
        boost::lexical_cast<std::string>(ret) +
        "\n"
        "This might be caused be a failed insert statement to the database.");
  return ret;
}

}  // namespace tools
}  // namespace votca

#endif
