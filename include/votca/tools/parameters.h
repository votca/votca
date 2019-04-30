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

#ifndef VOTCA_TOOLS_PARAMETERS_H
#define VOTCA_TOOLS_PARAMETERS_H

#include <boost/any.hpp>
#include <cassert>
#include <unordered_map>

namespace votca {
namespace tools {

/**
 * \breif Supported and Standardized parameter types
 **/
enum Parameter {
  Mass,
  Position,
  MoleculeId,
  ResidueId,
  Charge,
  Element,
  Symmetry,
  ResidueType,
  BeadId,
  BeadType,
  MoleculeType
};

/**
 * \brief Provides a means to standardise the constructors of different classes
 *
 * The purpose of this class is for use primarily in io readers, it provides
 * a means to standardize the use of templated classes in the readers.
 *
 * E.g. Say I have two atom classes
 *
 * class Atom1 {
 *   public:
 *     Atom1(Parameters parameters) :
 *        id_(parameters.get<int>(Parameter::BeadId)),
 *        mass_(parameters.get<double>(Parameter::Mass)){};
 *
 *   private:
 *    int id_;
 *    double mass_;
 * };
 *
 * class Atom2 {
 *   public:
 *     Atom2(Parameters parameters) :
 *        id_(parameters.get<int>(Parameter::BeadId)),
 *        element_(parameters.get<string>(Parameter::Element)){};
 *
 *   private:
 *    int id_;
 *    string element_;
 * };
 *
 * Pseudo code below, our file reader has a method that creates the atoms
 * (shown below), now because we have a standardized constructor the templated
 * file method can be used with either class without specialization.
 *
 * class FileReader {
 *
 *   template<class T>
 *   T CreateAtomOrBead {
 *
 *     string element = "C";
 *     int id = 1;
 *     double mass = 12.01;
 *
 *     Parameters parameters;
 *     parameters.set(Parameter::Element,element);
 *     parameters.set(Parameter::Mass,mass);
 *     parameters.set(Parameter::BeadId,id);
 *
 *     T atom_or_bead(parameters);
 *     return atom_or_bead;
 *   }
 *
 * };
 **/
class Parameters {

 public:
  void set(const Parameter parameter, boost::any value);

  template <class T>
  T get(const Parameter parameter) const;

 private:
  std::unordered_map<Parameter, boost::any> parameters;
};

void Parameters::set(const Parameter parameter, boost::any value) {
  parameters[parameter] = value;
}

template <class T>
T Parameters::get(const Parameter parameter) const {
  assert(parameters.count(parameter) &&
         "Parameter is not stored in Parameters class");
  assert(typeid(T) == parameters.at(parameter).type() &&
         "Cannot return boost any value from parameters class because it is "
         "not being cast to the correct type");
  return boost::any_cast<T>(parameters.at(parameter));
}

}  // namespace tools
}  // namespace votca
#endif  // VOTCA_TOOL_PARAMETERS_H
