#pragma once

#include <memory>
#include <vector>

// Mine
#include "testbeam/EventLoopBase.h"
#include "testbeam/Waveform.h"

class TDataContainer;

namespace myFuncs {
namespace testbeam {

//-----------------------------------------------------------
//
// Must setup the files using setupFiles() member function before using the class.
//
//-----------------------------------------------------------
class EventLoopBothChannels : public EventLoopBase {
public:
  EventLoopBothChannels(const std::string midasFilesPath, const std::string TOFfilesPath, const int runNum)
      : myFuncs::testbeam::EventLoopBase(midasFilesPath, TOFfilesPath, runNum) {}

  bool ProcessMidasEvent(TDataContainer &dataContainer) override final;

private:
  virtual bool processEvent(TDataContainer &dataContainer, Waveform &waveformCh1, Waveform &waveformCh15,
                            Waveform &waveformCombined) = 0;

  Waveform m_channel1waveform;
  long m_eventNumChannel1;
};

} // namespace testbeam
} // namespace myFuncs
