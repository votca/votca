#pragma once

#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>

/*!
 * \brief class used to time things.
 *
 *  To start the timer create an object, the elapsed duration can be queried
 *  anytime either as a chrono duration, a string or a full message.
 */
class Timer {
 public:
  Timer() { start = std::chrono::high_resolution_clock::now(); }

  Timer(std::string id) {
    start = std::chrono::high_resolution_clock::now();
    timer_id = id;
  }

  std::chrono::microseconds elapsedTime() {
    return std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::high_resolution_clock::now() - start);
  }

  std::string elapsedTimeAsString(int precision = 3) {
    std::stringstream msg;
    msg << duration2String(elapsedTime(), precision);
    return msg.str();
  }

  std::string elapsedTimeAsMessage() {
    std::stringstream msg;
    msg << "\nElapsed time";
    if (timer_id != "") {
      msg << " since start of timer" << timer_id;
    }
    msg << ": " << duration2String(elapsedTime(), 3) << std::endl;
    return msg.str();
  }

 private:
  std::chrono::high_resolution_clock::time_point start;
  std::string timer_id = "";

  std::string duration2String(std::chrono::microseconds ms, int precision = 3) {
    std::stringstream msg;
    using namespace std::chrono;
    using days = duration<int, std::ratio<86400>>;
    auto d = duration_cast<days>(ms);
    ms -= d;
    auto h = duration_cast<hours>(ms);
    ms -= h;
    auto m = duration_cast<minutes>(ms);
    ms -= m;
    auto s = duration_cast<seconds>(ms);
    ms -= s;

    long fs_count = 0;
    int printPrecision = 3;
    if (precision > 6) {
      fs_count = ms.count();
      printPrecision = 9;
    } else if (precision > 3) {
      fs_count = duration_cast<microseconds>(ms).count();
      printPrecision = 6;
    } else {
      fs_count = duration_cast<milliseconds>(ms).count();
      printPrecision = 3;
    }

    char fill = msg.fill('0');
    if (d.count()) {
      msg << d.count() << "d ";
    }
    if (d.count() || h.count()) {
      msg << std::setw(2) << h.count() << "h ";
    }
    if (d.count() || h.count() || m.count()) {
      msg << std::setw(d.count() || h.count() ? 2 : 1) << m.count() << "m ";
    }
    msg << std::setw(d.count() || h.count() || m.count() ? 2 : 1) << s.count();
    if (fs_count > 0) {
      msg << "." << std::setw(printPrecision) << fs_count;
    }
    msg << "s";
    msg.fill(fill);
    return msg.str();
  }
};

#endif