#pragma once

#include <vector>
#include <memory> //For unique ptr

//Root
#include "Rtypes.h"
#include "TChain.h"


//Mine
#include "testbeam/testbeamFuncs.h"

class TFitResultPtr;
class TH1D;

namespace myFuncs {
namespace testbeam{

//-----------------------------------------------------------
//
//-----------------------------------------------------------
class TOFtiming
{
public:
	TOFtiming(const std::vector<std::string>& filenames);
	
	std::shared_ptr<TChain> getChain() {return m_chain;}

	double eventNumber() const {return m_eventNumber;}
	double ch4Time() const {return m_ch4Time;}
	double ch6Time() const {return m_ch6Time;}
	double ch12Time() const {return m_ch12Time;}
	double ch13Time() const {return m_ch13Time;}
	double ch4TimeError() const {return m_ch4Error;}
	double ch6TimeError() const {return m_ch6Error;}
	double ch12TimeError() const {return m_ch12Error;}
	double ch13TimeError() const {return m_ch13Error;}
	
	
	 Long64_t getEntries() const {
		if(!m_chain) return 0;
		return m_chain->GetEntries();
	}
	
	void getEntry(const size_t entry) {
		if(m_chain) m_chain->GetEntry(entry);
	}
	
	//Do linear fit and return fit result
	//Assumes all channel times and errors are filled
	TFitResultPtr getTOF(const bool boundSpeed) const;
	
	//GetEntry entry and perform TOF fit. Return fit result
	TFitResultPtr getTOF(const Long64_t entry, const bool boundSpeed);
	
	//Downstream time - upstream time (ignores S0)
	//Assumes all channel times and errors are filled
	double getSimpleTOF() const { return (m_ch12Time + m_ch13Time)/2. - m_ch4Time; }
	
	//Downstream time - upstream time (ignores S0)
	//Does getEntry before the calculation
	double getSimpleTOF(const Long64_t entry) { 
		getEntry(entry);
		return getSimpleTOF(); 
	}
	
	double getC20() const {return m_C20;}
	double getC21() const {return m_C21;}
	
	double getT0() const {return m_ch4Time + m_C20;}
	double getT1() const {return m_ch6Time + m_C21;}
	double getT2() const {return 0.5 * ( m_ch12Time + m_ch13Time );}
	
	double getT0Error() const {return m_ch4Error;}
	double getT1Error() const {return m_ch6Error;}
	double getT2Error() const {return 0.5 * std::sqrt( m_ch12Error*m_ch12Error + m_ch13Error*m_ch13Error ); }
	
	/*constexpr - might be possible with c++14*/ double getX0() const {return -myFuncs::testbeam::c_downstreamTOF2upstreamTOFdistance;}
	/*constexpr - might be possible with c++14*/ double getX1() const {return -myFuncs::testbeam::c_downstreamTOF2S0distance;}
	/*constexpr - might be possible with c++14*/ double getX2() const {return 0.0;}
	
	//--------------------------------------------------------------------------------------
	//Calculate them by creating a simpleTOF histogram, if they haven't been calculated yet.
	double getElectronSimplTOFmean();
	double getElectronSimplTOFsigma();
	double getMuonSimplTOFmean();
	double getMuonSimplTOFsigma();
	double getPionSimplTOFmean();
	double getPionSimplTOFsigma();
	//--------------------------------------------------------------------------------------
	
	//return | mean(simpleTOF) - simpleTOF | < numSigmas * sigma(simpleTOF)
	bool isSimpleElectron(const double numSigmas);
	bool isSimpleMuon(const double numSigmas);
	bool isSimplePion(const double numSigmas);
	
	//Same, but does getEntry before the calculation
	bool isSimpleElectron(const Long64_t entry, const double numSigmas);
	bool isSimpleMuon(const Long64_t entry, const double numSigmas);
	bool isSimplePion(const Long64_t entry, const double numSigmas);
	
	std::string getTOFfunctionString() const {return m_TOFfunctionString;}
private:
	const double m_C20 = -20.654; //ns. C2 - C0
	const double m_C21 = -18.619; //ns. C2 - C1
	
	const double m_minElectronSimpleTOF = -11.0; // ns
	const double m_maxElectronSimpleTOF = -9.75; // ns
	
	Long64_t m_eventNumber; 
	double m_ch4Time;
	double m_ch6Time;
	double m_ch12Time;
	double m_ch13Time;

	double m_ch4Error;
	double m_ch6Error;
	double m_ch12Error;
	double m_ch13Error;
	
	double m_electronSimpleTOFmean;
	double m_electronSimpleTOFsigma;
	double m_muonSimpleTOFmean;
	double m_muonSimpleTOFsigma;
	double m_pionSimpleTOFmean;		
	double m_pionSimpleTOFsigma;
	
	bool m_calculatedSimpleTOFmean_sigma; //Keeps track if simpleTOF mean and sigma have been calculated for different particles
	
	const std::string m_treeName;
	const std::vector<std::string> m_branchNames;
	const std::string m_TOFfunctionString;
	std::vector<void*> m_pointers;
	const std::vector<std::string> m_filenames;
	std::shared_ptr<TChain> m_chain;
	
	//-----------------------------------------------------------
	//Fit hist with a Gaussian, initially in the range (minInitialRange, maxInitialRange).
	//Fit a new Gaussian to the range +- 2 sigma of the previos range and return fit results.
	//-----------------------------------------------------------
	TFitResultPtr fitSimpleTOF(TH1D hist, const double minInitialRange, const double maxInitialRange);
	
	//-----------------------------------------------------------
	//Calculate the mean and sigma of Gaussians for simpleTOF for electrons, muons and pions.
	//At the end, makes sure the last entry loaded is the same as before the function was called.
	//-----------------------------------------------------------
	void calculateSimpleTOFmean_sigma();
};


}//testbeam namespace
}//myFuncs namespace