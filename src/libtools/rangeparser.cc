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

#include <boost/lexical_cast.hpp>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <votca/tools/rangeparser.h>
#include <votca/tools/tokenizer.h>

namespace votca {
namespace tools {

RangeParser::RangeParser()
//: _has_begin(false) , _has_end(false)
{}

void RangeParser::Parse(string str) {
  // remove all spaces in string
  string::iterator it =
      remove_if(str.begin(), str.end(), bind2nd(equal_to<char>(), ' '));
  str = string(str.begin(), it);

  // tokenize string
  Tokenizer tok(str, ",");
  Tokenizer::iterator bl;

  for (bl = tok.begin(); bl != tok.end(); ++bl) ParseBlock(*bl);

  //    list<block_t>::iterator iter;
  //    for(iter=_blocks.begin();iter!=_blocks.end();++iter) {
  //        cout << (*iter)._begin << ":" << (*iter)._stride << ":" <<
  //        (*iter)._end << endl;
  //    }
}

void RangeParser::ParseBlock(string str) {
  Tokenizer tokenizer(str, ":");
  vector<string> toks;

  block_t block;
  block._stride = 1;

  tokenizer.ToVector(toks);
  if (toks.size() > 3 || toks.size() < 1) {
    throw runtime_error("invalid range");
  }

  block._begin = block._end = ToNumber(toks[0]);

  if (toks.size() == 2) block._end = ToNumber(toks[1]);

  if (toks.size() == 3) {
    block._stride = ToNumber(toks[1]);
    block._end = ToNumber(toks[2]);
  }

  if (block._begin * block._stride > block._end * block._stride) {
    throw runtime_error(
        string("invalid range " + str +
               ": begin, end and stride do not form a closed interval"));
  }

  _blocks.push_back(block);
}

int RangeParser::ToNumber(string str) { return boost::lexical_cast<int>(str); }

RangeParser::iterator& RangeParser::iterator::operator++() {
  _current += (*_block)._stride;
  if (_current > (*_block)._end) {
    ++_block;
    if (_block != _parent->_blocks.end())
      _current = (*_block)._begin;
    else
      _current = -1;
  }
  return *this;
}

}  // namespace tools
}  // namespace votca
