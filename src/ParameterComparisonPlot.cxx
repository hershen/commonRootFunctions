#include "ParameterComparisonPlot.h"

//STL
#include <iostream>

//ROOT
#include "TPad.h"
#include "TPaveText.h"

namespace myFuncs {
	
ParameterComparisonPlot::ParameterComparisonPlot() :
	TGraphErrors(), 
	m_drawValueLabels(true)
	{}

	ParameterComparisonPlot::~ParameterComparisonPlot() {
		
		//Delete value labels
		for(auto& label : m_valueLabels)
			delete label;
	
		//Delete y axis labels
		for(auto& label : m_yAxisLabels)
			delete label;
	}
	
void ParameterComparisonPlot::addEntry(const std::string& experiment, const double parValue, const double parError, const bool drawRectangle) {
	m_experiments.push_back(experiment);
	m_parValues.push_back(parValue);
	m_parErrors.push_back(parError);
	m_drawRectangles.push_back(drawRectangle);
	
}

double ParameterComparisonPlot::getYlengthPerPoint() const {
	const double yAxisRange = this->GetYaxis()->GetXmax() - this->GetYaxis()->GetXmin();
	return yAxisRange / this->GetN();
}

void ParameterComparisonPlot::drawValueLabels(const double textLabelSize) {
	
		const double xAxisRange = this->GetXaxis()->GetXmax() - this->GetXaxis()->GetXmin();

		for(int iPoint = 0; iPoint < this->GetN(); ++iPoint) {
			double x,y;
			this->GetPoint(iPoint, x, y);
			const double error = this->GetErrorX(iPoint);
			
			//This is a memory leak, but I plan to change this into regular objects, so I'm not bothering
			myFuncs::PaveText* pt = new myFuncs::PaveText(x, y - 0.25*getYlengthPerPoint(), x + 0.5*xAxisRange, y + 0.35*getYlengthPerPoint(), "NB");
			
			pt->SetTextSize(textLabelSize * gPad->GetWNDC() / 1.5);
			
			pt->SetTextAlign(kHAlignLeft + kVAlignBottom); //align so bottom of text is touching the central horizontal line

			pt->SetMargin(0.02); //Reduce marge (default 0.05) so that text is closer to left
			
			//Use stringstreams so that there are no trailing zeros in floating numbers
			std::stringstream ssX;
			ssX << x;
			std::stringstream ssError;
			ssError << error;
			
			pt->AddText( ( ssX.str() + "#pm" + ssError.str()).data() );
			pt->Draw();		
			m_valueLabels.push_back(pt);
		}
}

void ParameterComparisonPlot::drawYaxisLabels(const double textLabelSize) {
	//Couldn't get this to work with this->GetYaxis()->ChangeLabel() because root would add/remove labels depending on the number of points

	for(int iPoint = 0; iPoint < this->GetN(); ++iPoint) {
		
		double x,y;
		this->GetPoint(iPoint, x, y);

		//This is a memory leak, but I plan to change this into regular objects, so I'm not bothering
		myFuncs::PaveText* pt = new myFuncs::PaveText(myFuncs::xNdcToUser(0), y - 0.25*getYlengthPerPoint(), myFuncs::xNdcToUser(gPad->GetLeftMargin()), y + 0.35*getYlengthPerPoint(), "NB");
		
		pt->AddText((m_experiments[iPoint]).data());
		
		pt->SetTextSize(textLabelSize * gPad->GetWNDC());
		
		pt->SetTextAlign(kHAlignRight + kVAlignCenter);

		pt->SetMargin(0.02); //Reduce marge (default 0.05) so that text is closer to axis
		
		pt->SetTextFont( this->GetXaxis()->GetLabelFont() );
		
		pt->Draw();
		
		m_yAxisLabels.push_back(pt);
	}	
}

void ParameterComparisonPlot::drawRectangle(const int entry) {
	
	//Get last point's information
	double xParValue, y;
	this->GetPoint(entry, xParValue, y);
	const double xParError = this->GetErrorX(entry);
	
	gPad->Update();

	//Calculate tick marks length in order for rectangle not to be over the tick marks
	//Tick marks of x axis depend on Y axis length!!!
	const double yAxisRange = this->GetYaxis()->GetXmax() - this->GetYaxis()->GetXmin();
	const double xTickLength = this->GetXaxis()->GetTickLength()*yAxisRange;
		
	TPaveText* rectangle = new TPaveText( xParValue - xParError, gPad->GetUymin() + xTickLength, xParValue + xParError, gPad->GetUymax() - xTickLength,"NB" );
	
	rectangle->Draw();
	
	m_rectangles.push_back(rectangle);
}

void ParameterComparisonPlot::Draw(const std::string& options) {
	
	//Populate TGraphErrors
	for(uint iEntry = 0; iEntry < m_experiments.size(); ++iEntry) {
		
		const int entryYvalue = m_experiments.size() - iEntry;
		
		//Parameter value
		this->SetPoint(iEntry, m_parValues[iEntry], entryYvalue);
		
		//Parameter error
		this->SetPointError(iEntry, m_parErrors[iEntry], 0);
	}
	
	
	//////////////////////////////////////
	//Drawing order is important!!!
	//////////////////////////////////////
	
	//Draw axis of TGraph (by the internal histogram
	this->GetHistogram()->Draw();
	
	//Draw rectangles first so that the graph's points are drawn on top of it
	for(uint iEntry = 0; iEntry < m_experiments.size(); ++iEntry) { 
		if(m_drawRectangles[iEntry])
			drawRectangle(iEntry);
	}
	TGraphErrors::Draw(("P" + options).data());
	
	//Keep only X axis main ticks. Remove secondary ticks
	this->GetXaxis()->SetNdivisions((this->GetXaxis()->GetNdivisions() % 100));
	
	//Supress ticks on axis
	this->GetYaxis()->SetTickLength(0);
	
	this->SetMarkerStyle(kFullCircle);
	
	//This might be better set from script file using this class
	this->SetLineWidth(2);

	//Remove root generated y axis labels
	const double textLabelSize = this->GetYaxis()->GetLabelSize();
	this->GetYaxis()->SetLabelSize(0);
	
	//Set y axis labels as experiment names
	drawYaxisLabels(textLabelSize);
	
	//Add labels above datapoints
	if(m_drawValueLabels) 
		drawValueLabels(textLabelSize);
}

}//myFuncs namespace