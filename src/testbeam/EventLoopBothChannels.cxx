#include "testbeam/EventLoopBothChannels.h"

// Rootana
// Prevent unused variable warning from this file
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#include "TV1730RawData.hxx"
#pragma GCC diagnostic pop

// Mine
#include "testbeam/testbeamFuncs.h"

using namespace myFuncs::testbeam;

bool EventLoopBothChannels::ProcessMidasEvent(TDataContainer &dataContainer) {
  TV1730RawData *v1730 = dataContainer.GetEventData<TV1730RawData>("V730");
  if (!v1730)
    return false;

  for (auto measurement : v1730->GetMeasurements()) {

    // Choose only crystal channels
    const int channel = measurement.GetChannel();
    if (!isCrystalChannel(channel))
      continue;

    Waveform waveform(std::move(measurement.getSamples()), 2.0);

    const long eventNum = dataContainer.GetMidasEvent().GetSerialNumber();

    if (channel == 1) {
      m_eventNumChannel1 = eventNum;
      m_channel1waveform = waveform;

    } else if (channel == 15) {

      // combine channels
      if (eventNum == m_eventNumChannel1) {
        Waveform channel1aligned = m_channel1waveform;
        channel1aligned.timeShift(myFuncs::testbeam::c_channel1_channel15_timeDifference);
        Waveform combinedWaveforms = waveform + channel1aligned;
        return processEvent(dataContainer, m_channel1waveform, waveform, combinedWaveforms);
      }

      // This is channel 15 processing
      //
    }
  }

  // If you're here, couldn't find channel 1 and channel 15.
  return false;
}
