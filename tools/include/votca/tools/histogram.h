/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

#ifndef VOTCA_TOOLS_HISTOGRAM_H
#define VOTCA_TOOLS_HISTOGRAM_H

// Standard includes
#include <cmath>
#include <limits>

// Local VOTCA includes
#include "table.h"

namespace votca {
namespace tools {

/**
 *  \brief class to generate histograms
 *
 *  This class produces a bin centered histogram out of a vector of values.
 *  The min and max values that are provided with the initialization function
 *  will place bins at those positions. For example if I want to create a
 *  histogram for 0 to 5 with bins centered at 0 and 5, I could do:
 *
 *      Histogram hn;
 *      double min_value = 0.0;
 *      double max_value = 5.0;
 *      Index number_bins  = 6;
 *      hn.Initialize(min_value, max_value, number_bins);
 *
 *  This would produce a histogram with bins centered at the positions shown.
 *  below. The upper and lower bounds for each bin are also shown:
 *
 *  Bin:         1   2   3   4   5   6
 *  Pos:         0.0 1.0 2.0 3.0 4.0 5.0
 *  Upper Bound: 0.5 1.5 2.5 3.5 4.5 5.5
 *  Lower Bound:-0.5 0.5 1.5 2.5 3.5 4.5
 *
 *  Step Size: 1.0
 *
 *  You also have the option to change the histogram to a periodic bin centered
 *  histogram by setting the periodic value to true;
 *
 *      hn.setPeriodic(true);
 *
 *  This has the effect of rearranging the step size. The new histogram would
 *  look like this:
 *
 *  Bin:         1     2     3     4     5     6     1
 *  Pos:         0.000 0.833 1.667 2.500 3.333 4.167 5.000
 *  Upper Bound: 0.416 1.249 2.082 2.910 3.750 4.580 0.416
 *  Upper Bound: 4.580 0.416 1.249 2.082 2.910 3.750 4.580
 *
 *  Step Size: 0.833
 *
 *  Note that the bin in position 1 cycles back after bin 6 and that the step
 *  size is smaller this is a result of maintaining the same number of bins to
 *  be clear bin 1 is responsible for any values that fall between
 *
 *  0.00 - 0.416 and 4.580 - 5.00
 *
 */
class Histogram {
 public:
  /**
   * \brief Initialize the Histogram
   * @param min lower bound of interval
   * @param max upper bound of interval
   * @param nbins number of bins
   */
  void Initialize(double min, double max, Index nbins);

  /**
   * \brief process a data point
   * \param v value of this point
   * \param scale weighting of this point, bin of v is increased by scale
   * instead of 1
   */
  void Process(const double &v, double scale = 1.0);

  /**
      \brief process a range of data using iterator interface
   */
  template <typename iterator_type>
  void ProcessRange(const iterator_type &begin, const iterator_type &end);

  /**
   * \brief get the lower bound of the histogram intervaö
   * \return lower limit of interval
   */
  double getMin() const { return min_; }

  /**
   * \brief get the upper bound of the histogram intervaö
   * \return upper limit of interval
   */
  double getMax() const { return max_; }

  /**
   * \brief Get number of grid points
   * \return number of grid poitns
   */
  Index getNBins() const { return nbins_; }

  /**
   * \brief Get the count of the bin with the fewest counts
   * \return minimum counts
   */
  double getMinBinVal() const;

  /**
   * \brief Get the count of the bin with the maximum counts
   * \return maximum counts
   */
  double getMaxBinVal() const;
  /**
   * \brief Given the bin number get the Inverval bounds
   * @param[in] bin - pass in the index of the bin
   * \return pair with upper and lower bounds of the interval of the bin
   */
  std::pair<double, double> getInterval(Index bin) const;

  /**
   * \brief get the grid of histogram
   * \return step per bin
   */
  double getStep() const { return step_; }

  /**
   * \brief normalize the histogram that the integral is 1
   */
  void Normalize();

  /**
   * \brief clear all data
   */
  void Clear();

  /**
   * \brief get access to content of histogram
   * \return table object with bins in x and values in y
   */
  Table &data() { return data_; }
  const Table &data() const { return data_; }

  /**
   * \brief set whether interval is periodic
   * \param periodic is periodic
   */
  void setPeriodic(bool periodic) { periodic_ = periodic; }

 private:
  void Initialize_();
  double min_ = 0;
  double max_ = 0;
  double step_ = 0;
  bool periodic_ = false;
  Index nbins_ = 100;
  Table data_;
};

inline std::ostream &operator<<(std::ostream &out, Histogram &h) {
  out << h.data();
  return out;
}

template <typename iterator_type>
inline void Histogram::ProcessRange(const iterator_type &begin,
                                    const iterator_type &end) {
  for (iterator_type iter = begin; iter != end; ++iter) {
    Process(*iter);
  }
}
}  // namespace tools
}  // namespace votca
#endif  // VOTCA_TOOLS_HISTOGRAM_H
