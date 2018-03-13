/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE histogramnew_test
#include <boost/test/unit_test.hpp>
#include <boost/lexical_cast.hpp>
#include <exception>
#include <iostream>
#include <votca/tools/histogramnew.h>
#include <votca/tools/table.h>

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(histogramnew_test)

BOOST_AUTO_TEST_CASE(create_test) {
  HistogramNew hn;
}

BOOST_AUTO_TEST_CASE(init_test) {
  HistogramNew hn;
  double min_v = 1.2;
  double max_v = 102.0;
  hn.Initialize(min_v,max_v);
}

BOOST_AUTO_TEST_CASE(step_test) {
  HistogramNew hn;
  double min_v = 1;
  double max_v = 9;
  hn.Initialize(min_v,max_v,8);
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(hn.getStep()*10),10);
}

BOOST_AUTO_TEST_CASE(nbins_test) {
  HistogramNew hn;
  double min_v = 1;
  double max_v = 9;
  hn.Initialize(min_v,max_v,8);
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(hn.getNBins()*10),80);
}

BOOST_AUTO_TEST_CASE(Process_test){
  HistogramNew hn;
  double min_v = 0.0;
  double max_v = 10.0;
  hn.Initialize(min_v,max_v,10);
  vector<double> data;
  for(double x=0;x<10;++x){
    data.push_back(x);
  }
  hn.ProcessRange(data.begin(),data.end());
  hn.Process(4.5);
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(hn.getStep()*10),10);
  auto dat = hn.data();
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(dat.y(0)),1);
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(dat.y(1)),1);
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(dat.y(2)),1);
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(dat.y(3)),1);
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(dat.y(4)),2);
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(dat.y(5)),1);
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(dat.y(6)),1);
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(dat.y(7)),1);
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(dat.y(8)),1);
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(dat.y(9)),1);
} 

BOOST_AUTO_TEST_CASE(minmax_test){
  HistogramNew hn;
  double min_v = 0.0;
  double max_v = 10.0;
  hn.Initialize(min_v,max_v,10);
  vector<double> data;
  for(double x=0;x<9;++x){
    data.push_back(x);
  }
  hn.ProcessRange(data.begin(),data.end());
  hn.Process(4.5);
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(hn.getMinBinVal()),0); 
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(hn.getMaxBinVal()),2); 
} 

BOOST_AUTO_TEST_SUITE_END()
