

#include "votca/tools/timer.h"

#include <ostream>
#include <sstream>

Timer::Timer() { start_ = std::chrono::high_resolution_clock::now(); }

Timer::Timer(std::string id) {
  start_ = std::chrono::high_resolution_clock::now();
  timer_id_ = id;
}

std::chrono::microseconds Timer::elapsedTime() {
  return std::chrono::duration_cast<std::chrono::microseconds>(
      std::chrono::high_resolution_clock::now() - start_);
}

std::string Timer::elapsedTimeAsString(int precision) {
  std::stringstream msg;
  msg << durationToString(elapsedTime(), precision);
  return msg.str();
}

std::string Timer::elapsedTimeAsMessage() {
    std::stringstream msg;
    msg << "\nElapsed time";
    if (timer_id_ != "") {
      msg << " since start of timer" << timer_id_;
    }
    msg << ": " << durationToString(elapsedTime(), 3) << std::endl;
    return msg.str();
  }

std::string Timer::durationToString(std::chrono::microseconds ms, int precision) {
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