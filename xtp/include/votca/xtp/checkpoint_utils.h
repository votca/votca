/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
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
#ifndef VOTCA_XTP_CHECKPOINT_UTILS_H
#define VOTCA_XTP_CHECKPOINT_UTILS_H

// Standard includes
#include <cstddef>
#include <cstring>
#include <mutex>
#include <string>

// Third party includes
#include <H5Cpp.h>

// VOTCA includes
#include <votca/tools/types.h>

namespace votca {
namespace xtp {

using CptLoc = H5::Group;

namespace checkpoint_utils {

inline H5::DataSpace StrScalar() { return H5::DataSpace(H5S_SCALAR); }

/*
 * =======================  HDF5 THREAD-SAFETY NOTE  =======================
 *
 * IMPORTANT:
 * ----------
 * The HDF5 C++ API (H5::*) as provided by Homebrew and most system builds
 * is NOT thread-safe, even when different threads operate on different
 * files. In particular:
 *
 *   - HDF5 maintains global internal state
 *   - H5::Attribute, H5::DataSet, H5::Group, and H5::H5File call HDF5
 *     close functions (H5Aclose, H5Dclose, H5Gclose, H5Fclose) in their
 *     destructors
 *   - These destructors may run after leaving local scopes and MUST NOT
 *     execute concurrently with other HDF5 calls
 *
 * Failure to serialize *all* HDF5 activity (including destruction) leads
 * to hard-to-debug runtime errors such as:
 *
 *   - "Attribute::~Attribute - H5Aclose failed"
 *   - sporadic std::runtime_error from HDF5
 *   - segmentation faults under multithreading
 *
 * DESIGN DECISION:
 * ----------------
 * We enforce process-wide serialization of all HDF5 access using a single
 * global recursive mutex (Hdf5Mutex).
 *
 *  - The mutex is recursive because HDF5 calls are layered:
 *      CheckpointWriter -> CptTable -> HDF5
 *  - All HDF5-related code must lock this mutex
 *  - High-level operations (e.g. JobTopology::WriteToHdf5 / ReadFromHdf5)
 *    MUST hold the lock for their entire duration so that all HDF5 object
 *    destructors run while the mutex is held
 *
 * RULES FOR DEVELOPERS:
 * --------------------
 * 1) Any code that creates, reads, writes, or closes HDF5 objects MUST
 *    hold Hdf5Mutex().
 *
 * 2) If a function creates HDF5 objects whose destructors run at function
 *    exit, the lock MUST span the entire function scope.
 *
 * 3) Never call HDF5 (directly or indirectly) without holding this lock,
 *    even if the file names differ.
 *
 * 4) Do NOT remove or replace this mutex unless the entire HDF5 stack
 *    is redesigned (e.g. single-writer thread or MPI HDF5).
 *
 * This policy is REQUIRED for correct behavior under multithreading.
 * ========================================================================
 */
inline std::recursive_mutex& Hdf5Mutex() {
  static std::recursive_mutex m;
  return m;
}

// Declare some HDF5 data type inference stuff:
// Adapted from
// https://github.com/garrison/eigen3-hdf5/blob/2c782414251e75a2de9b0441c349f5f18fe929a2/eigen3-hdf5.hpp#L18

template <typename T>
struct InferDataType;

template <>
struct InferDataType<float> {
  static const H5::DataType* get(void) { return &H5::PredType::NATIVE_FLOAT; }
};

template <>
struct InferDataType<double> {
  static const H5::DataType* get(void) { return &H5::PredType::NATIVE_DOUBLE; }
};

template <>
struct InferDataType<int> {
  static const H5::DataType* get(void) { return &H5::PredType::NATIVE_INT; }
};

template <>
struct InferDataType<long int> {
  static const H5::DataType* get(void) { return &H5::PredType::NATIVE_LONG; }
};

template <>
struct InferDataType<unsigned> {
  static const H5::DataType* get(void) { return &H5::PredType::NATIVE_UINT; }
};

template <>
struct InferDataType<std::string> {
  static const H5::DataType* get(void) {

#if (defined(__GNUC__) && defined(__clang__))
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"
#elif (defined(__GNUC__) && !defined(__INTEL_COMPILER))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#elif (defined(__INTEL_COMPILER))
#pragma warning push
#pragma warning(disable : 1682)  // implicit conversion of a 64-bit integral
                                 // type to a smaller integral type
#endif
    static const H5::StrType strtype(H5T_C_S1, H5T_VARIABLE);
#if (defined(__GNUC__) && defined(__clang__))
#pragma clang diagnostic pop
#elif (defined(__GNUC__) && !defined(__INTEL_COMPILER))
#pragma GCC diagnostic pop
#elif (defined(__INTEL_COMPILER))
#pragma warning pop
#endif

    return &strtype;
  }
};

}  // namespace checkpoint_utils
}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_CHECKPOINT_UTILS_H
