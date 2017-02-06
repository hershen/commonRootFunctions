#include "testbeam/TOFtiming.h"

//Root
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1D.h"

//Mine
#include "fileFuncs.h"
#include "testbeam/constants.h"
#include "testbeam/RunDB.h"
#include "testbeam/testbeamFuncs.h"

using namespace myFuncs::testbeam;

TOFtiming::TOFtiming(const std::string pathToFiles, const int runNum):
  m_eventNumber(-1),
	m_ch4Time(0.0),
	m_ch6Time(0.0),
	m_ch12Time(0.0),
	m_ch13Time(0.0),
	m_ch4Error(0.0),
	m_ch6Error(0.0),
	m_ch12Error(0.0),
	m_ch13Error(0.0),
	m_electronSimpleTOFmean(0.0),
	m_electronSimpleTOFsigma(0.0),
	m_muonSimpleTOFmean(0.0),
	m_muonSimpleTOFsigma(0.0),
	m_pionSimpleTOFmean(0.0),
	m_pionSimpleTOFsigma(0.0),
	m_timeErrorScaling(1.0), // ns
	m_calculatedSimpleTOFmean_sigma(false),
	m_runNum(runNum),
	m_treeName("timingAndWaveforms"),
  m_branchNames({"eventNumber", "TOFch4Time", "TOFch4TimeError", "TOFch6Time", "TOFch6TimeError", "TOFch12Time", "TOFch12TimeError","TOFch13Time", "TOFch13TimeError"})
{
	//Define vector of pointer addresses
	m_pointers = {&m_eventNumber, &m_ch4Time, &m_ch4Error, &m_ch6Time, &m_ch6Error, &m_ch12Time, &m_ch12Error, &m_ch13Time, &m_ch13Error};
	
	//Create chain of files
	m_chain = std::shared_ptr<TChain>(myFuncs::openChain_setBranch( getFilesRelatedToRun(pathToFiles, getRunNum()), m_treeName, m_branchNames, m_pointers));
}


//Do linear fit and return fit result
TFitResultPtr TOFtiming::fitTOF(const bool boundSpeed) const
{
												//x values																								//y values																									 //y values errors
	TGraphErrors graph(3, (std::array<double,3>{getX0(), getX1(), getX2()}).data(), (std::array<double,3>{getT0(), getT1(), getT2()}).data(), 0, (std::array<double,3>{getT0Error(), getT1Error(), getT2Error()}).data() );
	
	//Fit function
	TF1 function("function",m_TOFfunctionString.data(), -4,1); //0 - intercept, 1 - speed
	
	std::string fitOptions = "E M EX0 S Q"; //E - better errors, M - Minos, EX0 - Don't use errors on x values, S - return smaprt ptr, Q - queit
	if(boundSpeed) {
		fitOptions += " B"; //B - use parameter limits in the fit
		function.SetParLimits(1, 0.0, 1);
	}
	
	function.SetParameter(1, 0.9);
	
	return graph.Fit(&function, fitOptions.data()); 
}

//GetEntry entry and perform TOF fit. Return fit result
TFitResultPtr TOFtiming::fitTOF(const Long64_t entry, const bool boundSpeed)
{
	getEntry(entry);
	return fitTOF(boundSpeed);
}

double TOFtiming::getX0() const {
	return -myFuncs::testbeam::c_downstreamTOF2upstreamTOFdistance - RunDB::instance()[getRunNum()].getDownstream2crystalCenterDistance();
}
	
double TOFtiming::getX1() const {
	return -myFuncs::testbeam::c_downstreamTOF2S0distance - RunDB::instance()[getRunNum()].getDownstream2crystalCenterDistance();
}

double TOFtiming::getX2() const {
	return -RunDB::instance()[getRunNum()].getDownstream2crystalCenterDistance();
}
	
bool TOFtiming::isSimpleElectron(const double numSigmas) {
	//Make sure mean and sigma have been calculated
	if(!m_calculatedSimpleTOFmean_sigma) calculateSimpleTOFmean_sigma();
	
	return std::abs(m_electronSimpleTOFmean - getSimpleTOF()) < numSigmas * m_electronSimpleTOFsigma;
}
bool TOFtiming::isSimpleMuon(const double numSigmas) {
	//Make sure mean and sigma have been calculated
	if(!m_calculatedSimpleTOFmean_sigma) calculateSimpleTOFmean_sigma();
		
	return std::abs(m_muonSimpleTOFmean - getSimpleTOF()) < numSigmas * m_muonSimpleTOFsigma;
}

bool TOFtiming::isSimplePion(const double numSigmas) {
	//Make sure mean and sigma have been calculated
	if(!m_calculatedSimpleTOFmean_sigma) calculateSimpleTOFmean_sigma();
		
	return std::abs(m_pionSimpleTOFmean - getSimpleTOF()) < numSigmas * m_pionSimpleTOFsigma;
}

bool TOFtiming::isSimpleElectron(const Long64_t entry, const double numSigmas) {
	getEntry(entry);
	return isSimpleElectron(numSigmas);
}

bool TOFtiming::isSimpleMuon(const Long64_t entry, const double numSigmas) {
	getEntry(entry);
	return isSimpleMuon(numSigmas);
}

bool TOFtiming::isSimplePion(const Long64_t entry, const double numSigmas) {
	getEntry(entry);
	return isSimplePion(numSigmas);
}

TFitResultPtr TOFtiming::fitSimpleTOF(TH1D hist, const double minInitialRange, const double maxInitialRange) {
	
	//Fit using initial range
	TFitResultPtr fitResult = hist.Fit("gaus","NQS","", minInitialRange, maxInitialRange );
	
	//Find +- 2sigma range
	const double minRange = fitResult->Parameter(1) - 2.0 * fitResult->Parameter(2);
	const double maxRange = fitResult->Parameter(1) + 2.0 * fitResult->Parameter(2);
	
	//return fit to +- 2 sigma range
	return hist.Fit("gaus","MENQS","", minRange, maxRange);
}

void TOFtiming::calculateSimpleTOFmean_sigma() {
	
	//If these have already been calculate, return
	if(m_calculatedSimpleTOFmean_sigma) return;
	
	//Save last entry retrieved.
	const Long64_t oldEntry = m_chain ? m_chain->GetReadEntry() : 0;
	
	//Fill simple TOF histogram
	TH1D simpleTOFhist("simpleTOFhist", "simpleTOFhist", std::pow(2,8), -13,0);
	Long64_t entries = getEntries();
	for(Long64_t iEntry = 0; iEntry< entries; ++iEntry) simpleTOFhist.Fill( getSimpleTOF(iEntry) );
	
	TFitResultPtr electronFitResult = fitSimpleTOF(simpleTOFhist, m_minElectronSimpleTOF, m_maxElectronSimpleTOF);
	m_electronSimpleTOFmean = electronFitResult->Parameter(1);
	m_electronSimpleTOFsigma = electronFitResult->Parameter(2);
	
	//Do for muons
	
	//Do for pions
	
	
	//Set flag that values have already been calculated
	m_calculatedSimpleTOFmean_sigma = true;
	
	//Get last entry that was retreived before we did anything.
	getEntry(oldEntry);
}

double TOFtiming::getElectronSimplTOFmean() {
	//Make sure mean and sigma have been calculated
	if(!m_calculatedSimpleTOFmean_sigma) calculateSimpleTOFmean_sigma();
		
	return m_electronSimpleTOFmean;
}
double TOFtiming::getElectronSimplTOFsigma() {
	//Make sure mean and sigma have been calculated
	if(!m_calculatedSimpleTOFmean_sigma) calculateSimpleTOFmean_sigma();
		
	return m_electronSimpleTOFsigma;
}
double TOFtiming::getMuonSimplTOFmean() {
	//Make sure mean and sigma have been calculated
	if(!m_calculatedSimpleTOFmean_sigma) calculateSimpleTOFmean_sigma();
		
	return m_muonSimpleTOFmean;
}
double TOFtiming::getMuonSimplTOFsigma() {
	//Make sure mean and sigma have been calculated
	if(!m_calculatedSimpleTOFmean_sigma) calculateSimpleTOFmean_sigma();
		
	return m_muonSimpleTOFsigma;
}
double TOFtiming::getPionSimplTOFmean() {
	//Make sure mean and sigma have been calculated
	if(!m_calculatedSimpleTOFmean_sigma) calculateSimpleTOFmean_sigma();
		
	return m_pionSimpleTOFmean;
}
double TOFtiming::getPionSimplTOFsigma() {
	//Make sure mean and sigma have been calculated
	if(!m_calculatedSimpleTOFmean_sigma) calculateSimpleTOFmean_sigma();
		
	return m_pionSimpleTOFsigma;
}


