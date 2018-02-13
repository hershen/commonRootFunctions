#pragma once

#include <memory>
#include <vector>

// Mine
#include "testbeam/EventLoopBase.h"

class TDataContainer;

namespace myFuncs {
namespace testbeam {

class Waveform;

//-----------------------------------------------------------
//
// Must setup the files using setupFiles() member function before using the class.
//
//-----------------------------------------------------------
class EventLoopBothChannels : public EventLoopBase {
public:
  EventLoopBothChannels(const std::string midasFilesPath, const std::string TOFfilesPath, const int runNum)
      : myFuncs::testbeam::EventLoopBase(midasFilesPath, TOFfilesPath, runNum) {}

  bool ProcessMidasEvent(TDataContainer& dataContainer) override final;

private:
  virtual bool processEvent(TDataContainer& dataContainer, Waveform& waveformCh1, Waveform& waveformCh15,
                            Waveform& waveformCombined) = 0;
};

} // namespace testbeam
} // namespace myFuncs
