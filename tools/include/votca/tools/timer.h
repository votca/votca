#pragma once

#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <iomanip>
#include <string>

/*!
 * \brief class used to time things.
 *
 *  To start the timer create an object, the elapsed duration can be queried
 *  anytime either as a chrono duration, a string or a full message.
 */
class Timer {
 public:
  Timer();
  Timer(std::string id);
  std::chrono::microseconds elapsedTime();
  std::string elapsedTimeAsString(int precision = 3);
  std::string elapsedTimeAsMessage();

 private:
  std::chrono::high_resolution_clock::time_point start_;
  std::string timer_id_ = "";

  std::string durationToString(std::chrono::microseconds ms, int precision = 3);
};

#endif