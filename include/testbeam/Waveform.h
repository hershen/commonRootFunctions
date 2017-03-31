#pragma once

#include <vector>
#include <cstdint>

namespace myFuncs {
namespace testbeam {

	

	
class Waveform {
public:
	Waveform(const std::vector<uint32_t> samples, const double dt);	
	
	//Get standard deviation between first and last
	double getStd(const unsigned int first, const unsigned int last) const;
	
	//Overloaded - Get standard deviation between m_samples[0] and m_samples[size * 12%]
	inline double getStd() const { return getStd(0, m_samples.size() * 0.12); }
	
	//Get mean between first and last
	double getMean(const unsigned int first, const unsigned int last) const;
	
	//Overloaded - Get mean deviation between m_samples[0] and m_samples[size * 12%]
	double getMean() const { return getMean(0, m_samples.size() * 0.12); }
	
	inline const std::vector<uint32_t>& getSamples() const {return m_samples;}
	std::vector<double>& getSamplesDouble() const;
	
	//Get vector of times.
	std::vector<double>& getTimes() const;
	
	//Return time at maximum and maximum.
	//Calculated with 2nd degree polynomial between first and last
	std::pair<double,double> getMaxPoly2(const unsigned int first, const unsigned int last) const;
	
	//Overloaded - First and last calculated 60-80%
	inline std::pair<double,double> getMaxPoly2() const { return getMaxPoly2( m_samples.size() * 0.6,  m_samples.size() * 0.8 ); }
	
	
	//Get simple amplitude = maximum sample - pedestal
	//Pedestal is taken to be getMean()
	//Because of this, it's not the most efficient because it loops on values again.
	double getSimpleAmplitude() const;
	
	
private:
	std::vector<uint32_t> m_samples;
	mutable std::vector<double> m_samples_double;
	mutable std::vector<double> m_times;
	const double m_dt; 
	
	
};


} //testbeam
} //myFuncs