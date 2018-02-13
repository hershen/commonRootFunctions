#include "testbeam/EventLoopBothChannels.h"

// Rootana
// Prevent unused variable warning from this file
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#include "TV1730RawData.hxx"
#pragma GCC diagnostic pop

// Mine
#include "testbeam/Waveform.h"
#include "testbeam/testbeamFuncs.h"

using namespace myFuncs::testbeam;

bool EventLoopBothChannels::ProcessMidasEvent(TDataContainer& dataContainer) {
  TV1730RawData* v1730 = dataContainer.GetEventData<TV1730RawData>("V730");
  if (!v1730)
    return false;

  Waveform channel1Waveform;

  for (auto measurement : v1730->GetMeasurements()) {

    // Choose only crystal channels
    const int channel = measurement.GetChannel();
    if (!isCrystalChannel(channel))
      continue;

    Waveform waveform(std::move(measurement.getSamples()), 2.0);

    if (channel == 1) {
      channel1Waveform = waveform;
    } else if (channel == 15) {
      // combine channels
      Waveform channel1aligned = channel1Waveform;
      channel1aligned.timeShift(myFuncs::testbeam::c_channel1_channel15_timeDifference);
      Waveform combinedWaveforms = waveform + channel1aligned;
      return processEvent(dataContainer, channel1Waveform /*ch1*/, waveform /*ch15*/, combinedWaveforms /*combined*/);
    }
  }

  // If you're here, couldn't find channel 1 and channel 15.
  return false;
}
