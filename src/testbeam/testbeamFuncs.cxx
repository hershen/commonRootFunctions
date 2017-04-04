#include "testbeam/testbeamFuncs.h"

//Root
#include "TGraphErrors.h"
#include "TPad.h"

//Mine
#include "fileFuncs.h"
#include "testbeam/RunDB.h"
#include "histFuncs.h" //for PaveText

namespace myFuncs {
namespace testbeam{

//The rightmost x position of the run + channel + event string
const double runChannelEventXpos = 0.64;

//-----------------------------------------------------------
//getNominalMomentumString
//-----------------------------------------------------------
std::string getNominalMomentumString(const int runNum) {
	int momentum = static_cast<int>(RunDB::instance()[runNum].getNominalBeamMomentum()); //Convert to int to get rid of decimal place
	return std::to_string(momentum) + " MeV/c";
}

//-----------------------------------------------------------
//getCrystalString
//-----------------------------------------------------------
std::string getCrystalString(const int runNum) {
	
	Crystal crystal = RunDB::instance()[runNum].getCrystal();
	if(crystal == Crystal::CsI_Tl_Belle) return "CsI(Tl) Belle";
	else if(crystal == Crystal::CsI_Tl_Babar) return "CsI(Tl) Babar";
	else if(crystal == Crystal::CsI_Ukrainian) return "CsI AMCRYS";
	else if(crystal == Crystal::CsI_Chinese) return "CsI SICCAS";
	else return "Unknown crystal";
}

//-----------------------------------------------------------
//getSourceDistanceString
//-----------------------------------------------------------
std::string getSourceDistanceString(const int runNum) {	
	int distance = static_cast<int>(RunDB::instance()[runNum].getSourceDistance()); //Convert to int to get rid of decimal place
	if(distance < 0 ) return "no ^{60}Co";
	return "^{60}Co " + std::to_string(distance) + " cm away";
}

std::vector<std::string> getFilesRelatedToRun(const std::string pathToFiles, const int runNum, const std::string extension)
{
	const std::vector<std::string> allFilenames = myFuncs::getFilesEndingWith( pathToFiles, extension );
	std::vector<std::string> filenames;
	//Find relevant files
	const std::string runString = std::to_string( runNum );
	for( const auto filename : allFilenames) if (filename.find(runString) != std::string::npos) filenames.push_back(filename);
	
	return filenames;
}

std::unique_ptr<TGraphErrors> getWaveformGraph(const std::vector<double>& voltage, const std::vector<double>& errors) {
	
	//Create time vector
	std::vector<double> times;
	times.reserve(voltage.size());
	for(size_t i = 0; i < voltage.size(); ++i)
		times.push_back(i * 2.0);
	
	std::unique_ptr<TGraphErrors> graph;
	
	if(errors.size() == 0) //no errors
		graph = std::unique_ptr<TGraphErrors>( new TGraphErrors( voltage.size(), times.data(), voltage.data(), nullptr, nullptr ) );
	else //with errors
		graph = std::unique_ptr<TGraphErrors>( new TGraphErrors( voltage.size(), times.data(), voltage.data(), nullptr, errors.data() ) );
	
	graph->SetTitle(";Time (ns); Voltage (ADC counts)");
	graph->SetMarkerSize(0.5);
	
	return graph;
	
}

std::unique_ptr<myFuncs::PaveText> getRunChannelEventPaveText(const int runNum, const int channelNum, const int eventNum) {
	std::unique_ptr<myFuncs::PaveText> pt(new myFuncs::PaveText(runChannelEventXpos) );
	pt->AddText( (getRunParamsString(runNum) + ", channelNum " + std::to_string(channelNum) + ", event "  + std::to_string(eventNum) ).data() );
	
 	return pt;
}


void drawWaveform(const std::vector<double>& voltages, const int runNum, const int channelNum, const size_t eventNum) {
		TCanvas c("c","c",0,0,1200,900);
		auto waveformGraph = myFuncs::testbeam::getWaveformGraph(voltages);
		waveformGraph->Draw("AP");
		
		auto eventInfoText = myFuncs::testbeam::getRunChannelEventPaveText(runNum, channelNum, eventNum);
		eventInfoText->DrawClone();
		
		gPad->WaitPrimitive();
}
			
}//testbeam namespace
}//myFuncs namespace
