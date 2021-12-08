#include "Sys_Tools.hpp"

double SYS_T::gen_randomD_closed(const double &min, const double &max)
{
  return ( rand() % 1000001 ) * 1.0e-6 * (max - min) + min;
}

double SYS_T::gen_randomD_open(const double &min, const double &max)
{
  return ( rand() % 999998 + 1 ) * 1.0e-6 * (max - min) + min;
}

int SYS_T::gen_randomI_closed(const int &min, const int &max)
{
  return ( rand() % (max - min + 1)) + min;
}


void SYS_T::get_memory_stats (MemoryStats &stats)
{
  stats.VmPeak = stats.VmSize = stats.VmHWM = stats.VmRSS = 0;

#if defined(__linux__)
  std::ifstream file("/proc/self/status");
  std::string line;
  std::string name;
  while (!file.eof())
  {
    file >> name;
    if (name == "VmPeak:")
      file >> stats.VmPeak;
    else if (name == "VmSize:")
      file >> stats.VmSize;
    else if (name == "VmHWM:")
      file >> stats.VmHWM;
    else if (name == "VmRSS:")
    {
      file >> stats.VmRSS;
      break; //this is always the last entry
    }

    getline(file, line);
  }
#endif
}

SYS_T::Timer::Timer()
{
  startedAt = 0;
  stoppedAt = 0;
}

SYS_T::Timer::~Timer()
{}

void SYS_T::Timer::Start()
{
  startedAt = clock();
}

void SYS_T::Timer::Stop()
{
  stoppedAt = clock();
}

void SYS_T::Timer::Reset()
{
  startedAt = 0;
  stoppedAt = 0;
}

double SYS_T::Timer::get_sec() const
{
  return (double)(stoppedAt - startedAt)/(double)CLOCKS_PER_SEC;
}

// EOF
