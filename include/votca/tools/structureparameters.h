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
#pragma once
#ifndef VOTCA_TOOLS_STRUCTUREPARAMETERS_H
#define VOTCA_TOOLS_STRUCTUREPARAMETERS_H

#include <boost/any.hpp>
#include <cassert>
#include <unordered_map>

namespace votca {
namespace tools {

/**
 * \breif Supported and Standardized parameter types
 **/
enum StructureParameter {
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
 *     Atom1(StructureParameters parameters) :
 *        id_(parameters.get<int>(StructureParameter::BeadId)),
 *        mass_(parameters.get<double>(StructureParameter::Mass)){};
 *
 *   private:
 *    int id_;
 *    double mass_;
 * };
 *
 * class Atom2 {
 *   public:
 *     Atom2(StructureParameters parameters) :
 *        id_(parameters.get<int>(StructureParameter::BeadId)),
 *        element_(parameters.get<string>(StructureParameter::Element)){};
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
 *     StructureParameters parameters;
 *     parameters.set(StructureParameter::Element,element);
 *     parameters.set(StructureParameter::Mass,mass);
 *     parameters.set(StructureParameter::BeadId,id);
 *
 *     T atom_or_bead(parameters);
 *     return atom_or_bead;
 *   }
 *
 * };
 **/
class StructureParameters {

 public:
  void set(const StructureParameter parameter, boost::any value) noexcept {
    parameters[parameter] = value;
  }

  bool ParameterExist(StructureParameter parameter) const noexcept {
    return parameters.count(parameter);
  }

  template <class T>
  T get(const StructureParameter parameter) const {
    assert(parameters.count(parameter) &&
           "StructureParameter is not stored in StructureParameters class");
    assert(typeid(T) == parameters.at(parameter).type() &&
           "Cannot return boost any value from parameters class because it is "
           "not being cast to the correct type");
    return boost::any_cast<T>(parameters.at(parameter));
  }

 private:
  std::unordered_map<StructureParameter, boost::any> parameters;
};

}  // namespace tools
}  // namespace votca
#endif  // VOTCA_TOOL_STRUCTUREPARAMETERS_H
