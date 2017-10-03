#include "testbeam/testbeamFuncs.h"

//Root
#include "TGraphErrors.h"
#include "TPad.h"
#include "TChain.h"

//Mine
#include "fileFuncs.h"
#include "testbeam/RunDB.h"
#include "histFuncs.h" //for PaveText

namespace myFuncs {
namespace testbeam{

//The rightmost x position of the run + channel + event string
const double runChannelEventXpos = 0.6;

//-----------------------------------------------------------
//getNominalMomentumString
//-----------------------------------------------------------
std::string getNominalMomentumString(const int runNum) {
	int momentum = std::lround(RunDB::instance()[runNum].getNominalBeamMomentum()); //Convert to int to get rid of decimal place
	return std::to_string(momentum) + " MeV/c";
}

//-----------------------------------------------------------
//getMeasuredMomentumString
//-----------------------------------------------------------
std::string getMeasuredMomentumString(const int runNum) {
	int momentum = std::lround(RunDB::instance()[runNum].getMeasuredBeamMomentum()); //Convert to int to get rid of decimal place
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
	int distance = std::lround(RunDB::instance()[runNum].getSourceDistance()); //Convert to int to get rid of decimal place
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
	pt->AddText( (getRunParamsString(runNum) + ", channel " + std::to_string(channelNum) + ", event "  + std::to_string(eventNum) ).data() );
	
 	return pt;
}


void drawWaveform(const std::vector<double>& voltages, const int runNum, const int channelNum, const size_t eventNum, const bool waitPrimitive) {
		TCanvas c("c","c",0,0,1200,900);
		auto waveformGraph = myFuncs::testbeam::getWaveformGraph(voltages);
		waveformGraph->Draw("AP");
		
		auto eventInfoText = myFuncs::testbeam::getRunChannelEventPaveText(runNum, channelNum, eventNum);
		eventInfoText->DrawClone();
		
		if(waitPrimitive)
			gPad->WaitPrimitive();
}

bool inElectronRange(const double beta, const int runNum, const double sigmasAway) {
	const double sigma  = electronBetaSigma.at(std::lround(RunDB::instance()[runNum].getNominalBeamMomentum()));
	return std::abs(1.0 - beta) < sigmasAway * sigma;
}

bool inMuonRange(const double beta, const int runNum, const double sigmasAway) {
	const double mean  = muonBetaMean.at(std::lround(RunDB::instance()[runNum].getNominalBeamMomentum()));
	const double sigma  = muonBetaSigma.at(std::lround(RunDB::instance()[runNum].getNominalBeamMomentum()));
	return std::abs(mean - beta) < sigmasAway * sigma;
}

bool inPionRange(const double beta, const int runNum, const double sigmasAway) {
	const double mean  = pionBetaMean.at(std::lround(RunDB::instance()[runNum].getNominalBeamMomentum()));
	const double sigma  = pionBetaSigma.at(std::lround(RunDB::instance()[runNum].getNominalBeamMomentum()));
	return std::abs(mean - beta) < sigmasAway * sigma;
}

bool isElectron(const double beta, const int runNum) {
	return inElectronRange(beta, runNum, 5.0) and not inMuonRange(beta, runNum, 4.0);	
}
bool isMuon(const double beta, const int runNum) {
		return inMuonRange(beta, runNum, 4.0) and not inPionRange(beta, runNum, 4.0) and not inElectronRange(beta, runNum, 5.0);	
}
bool isPion(const double beta, const int runNum) {
	return inPionRange(beta, runNum, 4.0) and not inMuonRange(beta, runNum, 4.0);	
}

int getV1730waveformLength(const int runNum) {
	if(isCsI(runNum))
		return 5000;
	else 
		return 17500;
}

bool isCsI(const int runNum) {
	const auto crystal = RunDB::instance()[runNum].getCrystal();
	if(crystal == Crystal::CsI_Ukrainian or crystal == Crystal::CsI_Chinese)
		return true;
	else 
		return false;
}

std::map<std::string, double> getGeantFileSimParamers(const std::string& filename) {
	
	const std::vector<std::string> branches{"beamParticle_pdgID", "nominalParticleMomentum_Mev_c", "particleMomentumResolution"};
	int pdg = 0; //Initialized so clang doesn't complain
	double primaryMeanMomentum = 0; //Initialized so clang doesn't complain
	double momentumResolution = 0; //Initialized so clang doesn't complain
	const std::vector<void*> pointers{&pdg, &primaryMeanMomentum, &momentumResolution};
	
	std::unique_ptr<TChain> chain(myFuncs::openChain_setBranch(filename, "simParameters", branches, pointers));
	
	chain->GetEntry(0);
	
	std::map<std::string, double> map;
	map["primaryPdg"] = pdg;
	map["primaryMeanMomentum"] = primaryMeanMomentum;
	map["momentumResolution"] = momentumResolution;
	
	return map;
	
}

}//testbeam namespace
}//myFuncs namespace
