#pragma once

//STL
#include <string>

//ROOT
#include "TGraphErrors.h"

//Mine
#include "histFuncs.h"

class TPaveText;

//TODO
//Change pointers to regular objects.
//Root fails to properly use the assignment operator for TPaveText (i.e. TPaveText pt; pt = TPaveText(x,y,x2,y2); pt->AddText("AA"); - there's a break segmentation).
//Therefore, either set x1, x2, y1, y2, options of TPaveText manually before using it (might not set TAttText options).
//Or try to see if the missing TAttText(22,0,gStyle->GetTextColor(),gStyle->GetTextFont(),0) from the constructor is the thing that causes the break segmentation
namespace myFuncs {
	

class PaveText;

	
class ParameterComparisonPlot : public TGraphErrors {
public:
	ParameterComparisonPlot();
	~ParameterComparisonPlot();
	
	void addEntry(const std::string& experiment, const double parValue, const double parError);
	
	void Draw(const std::string& options = "");
	
	std::vector<myFuncs::PaveText*> getValueLabels() const {return m_valueLabels;}
	
	void setDrawValueLabels(const bool drawValueLabels) {m_drawValueLabels = drawValueLabels;}
	bool getDrawValueLabels() const {return m_drawValueLabels;}
	
	void setDrawRectangleForLastEntry(const bool drawRectangle) {m_drawRectangleForLastEntry = drawRectangle;}
	bool getDrawRectangleForLastEntry() const {return m_drawRectangleForLastEntry;}
	
private:
	std::vector<std::string> m_experiments;
	std::vector<double> m_parValues;
	std::vector<double> m_parErrors;
	
	//Vector that holds the value lables (value +- error) above the markers
	//Can't be a shared pointer because root will delete them when program ends
	std::vector<myFuncs::PaveText*> m_valueLabels;
	
	//Vector that holds the value lables (value +- error) above the markers
	std::vector<myFuncs::PaveText*> m_yAxisLabels;
	
	TPaveText* m_lastEntryRectangle;
	
	bool m_drawValueLabels;
	
	bool m_drawRectangleForLastEntry;
	
	//Add text labels as y axis labels
	void drawYaxisLabels(const double textLabelSize);
	
	//Add text labels above points
	void drawValueLabels(const double textLabelSize);
	
	double getYlengthPerPoint() const;
	
	//Draw a rectangle that corresponds to last entry value +- it's error that goes all the way
	//from the bottom of the tpad to the top
	void drawRectangleForLastEntry();
};

}//myFuncs namespace