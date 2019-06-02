/*
 *            Copyright 2009-2019 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_MULTIARRAY_H
#define VOTCA_XTP_MULTIARRAY_H

#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
namespace votca {
namespace xtp {
typedef boost::multi_array<double, 3> tensor3d;
typedef boost::multi_array<double, 4> tensor4d;

typedef boost::multi_array_types::extent_range range;  //////////////////
typedef tensor3d::index index3d;                       /////////////////////
typedef tensor4d::index index4d;                       /////////////////////

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_MULTIARRAY_H
