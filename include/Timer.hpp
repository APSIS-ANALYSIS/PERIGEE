#ifndef TIMER_HPP
#define TIMER_HPP
// ============================================================================
// Timer.hpp
// ----------------------------------------------------------------------------
// This file defines a steady-clock timer for measuring elapsed wall-clock time.
//
// Author: Ju Liu
// Date Created: June 1 2026
// ============================================================================

#include <chrono>

namespace SYS_T
{

  class Timer
  {
    public:
      Timer() = default;
      ~Timer() = default;

      void Reset() noexcept
      {
        startedAt = std::chrono::steady_clock::time_point{};
        stoppedAt = std::chrono::steady_clock::time_point{};
        is_running = false;
        has_started = false;
      }

      void Start() noexcept
      {
        startedAt = std::chrono::steady_clock::now();
        stoppedAt = std::chrono::steady_clock::time_point{};
        is_running = true;
        has_started = true;
      }

      void Stop() noexcept
      {
        if (is_running)
        {
          stoppedAt = std::chrono::steady_clock::now();
          is_running = false;
        }
      }

      double get_sec() const noexcept
      {
        if (!has_started)
        {
          return 0.0;
        }

        const auto end_time = is_running ? std::chrono::steady_clock::now() : stoppedAt;
        return std::chrono::duration<double>(end_time - startedAt).count();
      }

      bool IsRunning() const noexcept
      {
        return is_running;
      }

    private:
      std::chrono::steady_clock::time_point startedAt{};
      std::chrono::steady_clock::time_point stoppedAt{};
      bool is_running{false};
      bool has_started{false};
  };

}

#endif
