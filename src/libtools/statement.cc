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
#include <string>
#include <votca/tools/statement.h>

namespace votca {
namespace tools {

using namespace std;

Statement::~Statement() { sqlite3_finalize(_stmt); }

template <>
void Statement::Bind(int col, const int &value) {
  if (sqlite3_bind_int(_stmt, col, value) != SQLITE_OK)
    throw std::runtime_error("sqlite_bind failed");
}

template <>
void Statement::Bind(int col, const double &value) {
  if (sqlite3_bind_double(_stmt, col, value) != SQLITE_OK)
    throw std::runtime_error("sqlite_bind failed");
}

template <>
int Statement::Column<int>(int col) {
  return sqlite3_column_int(_stmt, col);
}

template <>
double Statement::Column<double>(int col) {
  return sqlite3_column_double(_stmt, col);
}

template <>
string Statement::Column<string>(int col) {
  return string((const char *)sqlite3_column_text(_stmt, col));
}

template <>
void Statement::Bind(int col, const string &value) {
  if (sqlite3_bind_text(_stmt, col, value.c_str(), -1, NULL) != SQLITE_OK)
    throw std::runtime_error("sqlite_bind failed");
}

int Statement::Step() { return sqlite3_step(_stmt); }

void Statement::Reset() { sqlite3_reset(_stmt); }

}  // namespace tools
}  // namespace votca
