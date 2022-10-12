// adapted from Kitware kwsys, with BSD 3-Clause license
// https://gitlab.kitware.com/utils/kwsys/-/blob/master/SystemInformation.cxx

// Tested with:
// Windows (g++, clang++, icx, cl)
// MacOS (g++, clang++, icpc)
// Linux (g++, clang++, icpx)

// Compiler OS-detection macros
// https://sourceforge.net/p/predef/wiki/OperatingSystems/

#include <vector>
#include <cassert>
#include <bitset>
#include <limits>
#include <string>
#include <cstring>
#include <set>
#include <thread>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#elif defined (__APPLE__)
#include <sys/sysctl.h>
#elif defined(__OpenBSD__) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__DragonFly__)
#include <sys/sysctl.h>
#elif defined(__hpux)
#include <sys/param.h>
#include <sys/mpctl.h>
#elif defined(__HAIKU__)
#include <OS.h>
#elif __has_include(<unistd.h>)
#include <unistd.h>
#endif

unsigned int CPUCountWindows();
unsigned int ParseSysCtl();
unsigned int RetrieveInformationFromCpuInfoFile();
unsigned int QueryBSDProcessor();
unsigned int QueryHaikuInfo();
unsigned int QueryHPUXProcessor();
unsigned int QueryProcessorBySysconf();
unsigned int QueryThreads();

std::string ExtractValueFromCpuInfoFile(std::string buffer, const char* word,
  size_t& CurrentPositionInFile, size_t init = 0);

#ifdef __cplusplus
extern "C" {
#endif

unsigned int cpu_count(){

  unsigned int NumberOfPhysicalCPU = 0;

#if defined (_WIN32)
  NumberOfPhysicalCPU = CPUCountWindows();
#elif defined (__APPLE__)
  NumberOfPhysicalCPU = ParseSysCtl();
#elif defined(__OpenBSD__) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__DragonFly__)
  NumberOfPhysicalCPU = QueryBSDProcessor();
#elif defined(__linux) || defined(__CYGWIN__)
  NumberOfPhysicalCPU = RetrieveInformationFromCpuInfoFile();
#elif defined(__QNX__)
  // kwSys uses other kwSys functions for QNX. Is there a QNX library call to do this?
#elif defined(_AIX)
  // https://www.ibm.com/support/pages/determining-how-many-cpus-you-have-under-aix
  // looks like parsing text is required
#elif defined(__hpux)
  NumberOfPhysicalCPU = QueryHPUXProcessor();
#endif

  if (NumberOfPhysicalCPU == 0)
    NumberOfPhysicalCPU = QueryProcessorBySysconf();

  if (NumberOfPhysicalCPU == 0)
    NumberOfPhysicalCPU = QueryThreads();

  return NumberOfPhysicalCPU;

}

#ifdef __cplusplus
}
#endif


unsigned int CPUCountWindows(){

  unsigned int NumberOfPhysicalCPU = 0;

#ifdef _WIN32

  typedef BOOL(WINAPI * GetLogicalProcessorInformationType)(
    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION, PDWORD);
  static GetLogicalProcessorInformationType pGetLogicalProcessorInformation =
    (GetLogicalProcessorInformationType)GetProcAddress(
      GetModuleHandleW(L"kernel32"), "GetLogicalProcessorInformation");

  if (!pGetLogicalProcessorInformation) {
    return 0;
  }

  std::vector<SYSTEM_LOGICAL_PROCESSOR_INFORMATION> ProcInfo;
  {
    DWORD Length = 0;
    DWORD rc = pGetLogicalProcessorInformation(nullptr, &Length);
    assert(rc == 0);
    (void)rc; // Silence unused variable warning
    assert(GetLastError() == ERROR_INSUFFICIENT_BUFFER);
    ProcInfo.resize(Length / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION));
    rc = pGetLogicalProcessorInformation(&ProcInfo[0], &Length);
    assert(rc != 0);
    (void)rc; // Silence unused variable warning
  }

  typedef std::vector<SYSTEM_LOGICAL_PROCESSOR_INFORMATION>::iterator
    pinfoIt_t;
  for (pinfoIt_t it = ProcInfo.begin(); it != ProcInfo.end(); ++it) {
    SYSTEM_LOGICAL_PROCESSOR_INFORMATION PInfo = *it;
    if (PInfo.Relationship != RelationProcessorCore) {
      continue;
    }

    std::bitset<std::numeric_limits<ULONG_PTR>::digits> ProcMask(
      (unsigned long long)PInfo.ProcessorMask);
    unsigned int count = (unsigned int)ProcMask.count();
    if (count == 0) { // I think this should never happen, but just to be safe.
      continue;
    }
    NumberOfPhysicalCPU++;
  }

#endif

  return NumberOfPhysicalCPU;

}


unsigned int RetrieveInformationFromCpuInfoFile(){

  unsigned int NumberOfPhysicalCPU = 0;
  std::string buffer;

  FILE* fd = fopen("/proc/cpuinfo", "r");
  if (!fd) return 0;

  size_t fileSize = 0;
  while (!feof(fd)) {
    buffer += static_cast<char>(fgetc(fd));
    fileSize++;
  }
  fclose(fd);
  buffer.resize(fileSize - 2);
  // Number of logical CPUs (combination of multiple processors, multi-core
  // and SMT)
  size_t pos = buffer.find("processor\t");
  while (pos != std::string::npos) {
    pos = buffer.find("processor\t", pos + 1);
  }

  // Count sockets.
  size_t CurrentPositionInFile;
  std::set<int> PhysicalIDs;
  std::string idc = ExtractValueFromCpuInfoFile(buffer, "physical id", CurrentPositionInFile);
  while (CurrentPositionInFile != std::string::npos) {
    int id = atoi(idc.c_str());
    PhysicalIDs.insert(id);
    idc = ExtractValueFromCpuInfoFile(buffer, "physical id",
      CurrentPositionInFile, CurrentPositionInFile + 1);
  }

  uint64_t NumberOfSockets = PhysicalIDs.size();
  // Physical ids returned by Linux don't distinguish cores.
  // We want to record the total number of cores in NumberOfPhysicalCPU
  // (checking only the first proc)
  std::string Cores = ExtractValueFromCpuInfoFile(buffer, "cpu cores", CurrentPositionInFile);
  if (Cores.empty()) {
    // Linux Sparc is different
    Cores = ExtractValueFromCpuInfoFile(buffer, "ncpus probed", CurrentPositionInFile);
  }
  auto NumberOfCoresPerSocket = (unsigned int)atoi(Cores.c_str());
  NumberOfCoresPerSocket = std::max(NumberOfCoresPerSocket, 1u);
  NumberOfPhysicalCPU = NumberOfCoresPerSocket * (unsigned int)NumberOfSockets;

  return NumberOfPhysicalCPU;

}

unsigned int ParseSysCtl(){

  unsigned int NumberOfPhysicalCPU = 0;

#ifdef __APPLE__

  int N;
  size_t size = sizeof(N);

  if (sysctlbyname("hw.perflevel0.physicalcpu", &N, &size, nullptr, 0) == 0) {
    // Apple Silicon performance core count
    NumberOfPhysicalCPU = N;
  }
  else if (sysctlbyname("hw.physicalcpu", &N, &size, nullptr, 0) == 0) {
    // assumes heterogenous cores e.g. Intel Mac
    NumberOfPhysicalCPU = N;
  }

#endif

  return NumberOfPhysicalCPU;

}

unsigned int QueryHaikuInfo(){

  unsigned int NumberOfPhysicalCPU = 0;

#if defined(__HAIKU__)

  system_info info;
  get_system_info(&info);
  NumberOfPhysicalCPU = info.cpu_count;

#endif

  return NumberOfPhysicalCPU;

}

unsigned int QueryBSDProcessor(){

  unsigned int NumberOfPhysicalCPU = 0;

#if defined(__OpenBSD__) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__DragonFly__)

  int k;
  size_t sz = sizeof(k);
  int ctrl[2] = { CTL_HW, HW_NCPU };

  if (sysctl(ctrl, 2, &k, &sz, nullptr, 0) != 0) {
    return 0;
  }

  NumberOfPhysicalCPU = k;

#endif

  return NumberOfPhysicalCPU;

}


unsigned int QueryHPUXProcessor(){

  unsigned int NumberOfPhysicalCPU = 0;

#if defined(__hpux)
  int c = mpctl(MPC_GETNUMSPUS_SYS, 0, 0);
  if (c <= 0) {
    return 0;
  }

  NumberOfPhysicalCPU = c;
#endif

  return NumberOfPhysicalCPU;

}


unsigned int QueryProcessorBySysconf(){

  unsigned int NumberOfPhysicalCPU = 0;

#if defined(_SC_NPROCESSORS_ONLN)

 long c = sysconf(_SC_NPROCESSORS_ONLN);
 if (c > 0)
    NumberOfPhysicalCPU = static_cast<unsigned int>(c);

#endif

return NumberOfPhysicalCPU;

}

unsigned int QueryThreads(){
  // fallback, doesn't consider hyperthreading

  unsigned int NumberOfLogicalCPU = std::thread::hardware_concurrency();
  unsigned int NumberOfPhysicalCPU = NumberOfLogicalCPU;

  return NumberOfPhysicalCPU;

}


/** Extract a value from the CPUInfo file */
std::string ExtractValueFromCpuInfoFile(std::string buffer, const char* word,
  size_t & CurrentPositionInFile, size_t init)
{

  size_t pos = buffer.find(word, init);
  if (pos != std::string::npos) {
    CurrentPositionInFile = pos;
    pos = buffer.find(':', pos);
    size_t pos2 = buffer.find('\n', pos);
    if (pos != std::string::npos && pos2 != std::string::npos) {
      // It may happen that the beginning matches, but this is still not the
      // requested key.
      // An example is looking for "cpu" when "cpu family" comes first. So we
      // check that
      // we have only spaces from here to pos, otherwise we search again.
      for (size_t i = CurrentPositionInFile + strlen(word); i < pos;
           ++i) {
        if (buffer[i] != ' ' && buffer[i] != '\t') {
          return ExtractValueFromCpuInfoFile(buffer, word, CurrentPositionInFile, pos2);
        }
      }
      buffer.erase(0, pos + 2);
      buffer.resize(pos2 - pos - 2);
      return buffer;
    }
  }
  CurrentPositionInFile = std::string::npos;
  return "";
}
