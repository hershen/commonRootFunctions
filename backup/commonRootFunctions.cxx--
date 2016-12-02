#pragma once

#include <iostream>
#include <string.h>
#include <limits>
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TDatime.h"

using namespace std;
using namespace TMath;

//--------------------------------------------------------------------------------------------
//recieve 2 vectors and add their elements together
template <typename Type> 
std::vector<Type> add2vectors(std::vector<Type> vec1, std::vector<Type> vec2)
{
  std::vector<Type> outputVec;
  if (vec1.size() != vec2.size()) 
  {
      cout << "add2vectors: Vector sizes don't match. vec1 has " << vec1.size() << " elements. vect2 has " << vec2.size() << " elements. Aborting." << endl;
      return outputVec;
  }
  
  for(size_t  i = 0; i < vec1.size(); i++)
  {
    outputVec.push_back(vec1[i] + vec2[i]);    
  }
  return outputVec;
}
//--------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------
//recieve 2 vectors and add their elements together
template <typename Type1, typename Type2> 
double linearInterpolation(Type1 x1, Type2 y1, Type1 x2, Type2 y2, Type1 x)
{
  if( abs(x1 - x2) < 1e-9)
  {
    cout << "linearInterpolation: x1 = x2. Returning DBL_MAX" << endl;
    return std::numeric_limits<double>::max() ;    
  }
    
  return (y2 - y1)/(x2 - x1)*(x - x1) + y1;
}
//--------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------
//Read waveform of 4104 scope.
Int_t read4104waveform(TString fullFileName, std::vector<Double_t> *time,std::vector<Double_t> *ch1Voltage = 0, 
									 std::vector<Double_t> *ch2Voltage = 0,
									 std::vector<Double_t> *ch3Voltage = 0,
									 std::vector<Double_t> *ch4Voltage = 0,Int_t numHeaderLines = 14)
{
  //Sanity checks and clear vectors
  if(!time || (!ch1Voltage && !ch2Voltage && !ch3Voltage && !ch4Voltage) )
  {
    cout << "read4104waveform: non initialized vectors passed. Aborting" << endl;
    return 0;
  }
  else 
  {
    time->clear(); 
    if(ch1Voltage) ch1Voltage->clear();
    if(ch2Voltage) ch2Voltage->clear();
    if(ch3Voltage) ch3Voltage->clear();
    if(ch4Voltage) ch4Voltage->clear();
  }
  
  FILE* file = fopen(fullFileName,"r");
  if(file == NULL)
  {
    cout << "read4104waveform: Could not open file " << fullFileName << ". Aborting" << endl;
    return 0;
  }
  char line[160], dummy[5] = "    ", title1[5] = "    ", title2[5] = "    ", title3[5] = "    ", title4[5] = "    ";
  //Read file
  for(Int_t i=1; i<=numHeaderLines; i++) {0 != fgets(line,160,file);}
  
  //check which channels appear in the file
  0 != fgets(line,160,file);
  sscanf(line,"%4s,%3s,%3s,%3s,%3s",dummy,title1,title2,title3,title4);

  std::vector<Double_t> *firstColumnChannel = 0, *secondColumnChannel = 0, *thirdColumnChannel = 0, *fourthColumnChannel = 0; 
  
  if( !strcmp(title1,"CH4") && ch4Voltage) firstColumnChannel = ch4Voltage; 		// CH4
  else if( !strcmp(title1,"CH3") && ch3Voltage)                                         // Ch3   
  {
    firstColumnChannel = ch3Voltage;
    if( !strcmp(title2,"CH4") && ch4Voltage) secondColumnChannel = ch4Voltage;    	// CH3 CH4
  }
  else if( !strcmp(title1,"CH2") && ch2Voltage)						//CH2
  {
    firstColumnChannel = ch2Voltage;
    if( !strcmp(title2,"CH3") && ch3Voltage) 						//CH2 CH3
    {
      secondColumnChannel = ch3Voltage;
      if( !strcmp(title3,"CH4") && ch4Voltage) thirdColumnChannel = ch4Voltage; 	//CH2 CH3 CH4
    }
    else if(!strcmp(title2,"CH4") && ch4Voltage) secondColumnChannel = ch4Voltage;      //CH2 CH4
  }
  else if( !strcmp(title1,"CH1") && ch1Voltage)						//CH1
  {
    firstColumnChannel = ch1Voltage;
    if( !strcmp(title2,"CH2") && ch2Voltage) 						//CH1 CH2
    {
      secondColumnChannel = ch2Voltage;
      if( !strcmp(title3,"CH3") && ch3Voltage) 					 	//CH1 CH2 CH3
      {
	thirdColumnChannel = ch3Voltage;
	if( !strcmp(title4,"CH4") && ch4Voltage) fourthColumnChannel = ch4Voltage;	//CH1 CH2 CH3 CH4
      }
      else if( !strcmp(title3,"CH4") && ch4Voltage) thirdColumnChannel = ch4Voltage;	//CH1 CH2 CH4 
    }
    else if(!strcmp(title2,"CH3") && ch3Voltage)					//CH1 CH3
    {
      secondColumnChannel = ch3Voltage;
      if( !strcmp(title3,"CH4") && ch4Voltage) thirdColumnChannel = ch4Voltage;		//CH1 CH3 CH4
    }
    else if( !strcmp(title2,"CH4") && ch4Voltage) secondColumnChannel = ch4Voltage;	//CH1 CH4
  }
  else
  {
    cout << "read4104waveform: Requested channels don't corrispond to what's in the file " << fullFileName.Data() <<". Aborting" << endl;
    if(ch1Voltage) ch1Voltage->clear();
    if(ch2Voltage) ch2Voltage->clear();
    if(ch3Voltage) ch3Voltage->clear();
    if(ch4Voltage) ch4Voltage->clear();
    return 0;
  }
  
  Double_t tempTime = 0., column1Voltage = 0., column2Voltage = 0., column3Voltage = 0., column4Voltage = 0.;
  
  while (fgets(line,160,file) != 0 ) 
  {
    //if line begins with whitespace, carry on
    if (isspace(line[0])) continue;
    
    sscanf(line,"%lf,%lf,%lf,%lf,%lf",&tempTime,&column1Voltage,&column2Voltage,&column3Voltage,&column4Voltage);
    
    //sanity check
    if(abs(tempTime*pow(10,9)) > 5e5 || abs(column1Voltage) > 10 || abs(column2Voltage) > 10 || abs(column3Voltage) > 10 || abs(column4Voltage) > 10)
    {
      cout << "Invalid measurement : fullFileName = " << fullFileName.Data() << ", time = " << tempTime*pow(10,9) << " ns, column1Voltage = " << column1Voltage << "V, column2Voltage = " << column2Voltage << "V, column3Voltage = " << column3Voltage << "V, column4Voltage = " << column4Voltage << "V. Aborting" << endl;
      time->clear();  
      if( ch1Voltage) ch1Voltage->clear();
      if( ch2Voltage) ch2Voltage->clear();
      if( ch3Voltage) ch3Voltage->clear();
      if( ch4Voltage) ch4Voltage->clear();
      return 0;
    }
    
    time->push_back(tempTime*pow(10,9));
    if(firstColumnChannel) firstColumnChannel->push_back(column1Voltage);
    if(secondColumnChannel) secondColumnChannel->push_back(column2Voltage);
    if(thirdColumnChannel) thirdColumnChannel->push_back(column3Voltage);
    if(fourthColumnChannel) fourthColumnChannel->push_back(column4Voltage);
  }
  fclose(file);
  return 1;
}
//--------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------
//Check if file outputFullFile exists. If so, create a backup of it with current date and time.
//Return codes:
// 0 - There was a problem
// 1 - outputFullFile did not exist. Did nothing
// 2 - Created copy
Int_t checkIfFileExistsAndCreateBackup(TString outputFullFile)
{
  FileStat_t dummyBuf;
  if( !gSystem->GetPathInfo(outputFullFile.Data(),dummyBuf) )   //0 - file exsits, 1 - file doesn't exist
  {
    TDatime tdatime;
    TString currentDateTString;
    currentDateTString += "_";
    currentDateTString += tdatime.GetYear(); currentDateTString += "_";
    currentDateTString += tdatime.GetMonth(); currentDateTString += "_";
    currentDateTString += tdatime.GetDay(); currentDateTString += "-";
    currentDateTString += tdatime.GetHour(); currentDateTString += "_";
    currentDateTString += tdatime.GetMinute(); currentDateTString += "_";
    currentDateTString += tdatime.GetSecond(); 
    TString backupFileName = outputFullFile;
    backupFileName.Insert(backupFileName.Last('.'),currentDateTString.Data());
    //If file already exists, add '_' charecters to end of it
    while( !gSystem->GetPathInfo(backupFileName.Data(),dummyBuf) ) backupFileName.Insert(backupFileName.Last('.'),"_");  //0 - file exsits, 1 - file doesn't exist
    Int_t returnVal = 0;
//     cout << "outputFullFile.Data() = " << outputFullFile.Data() << ",  backupFileName.Data() = " << backupFileName.Data() << endl;
    returnVal = gSystem->CopyFile(outputFullFile.Data(), backupFileName.Data(), kFALSE);
    if(returnVal != 0)
    {
      cout << "checkIfFileExistsAndCreateBackup: Could not backup file " << outputFullFile.Data() << ". Aborting"  << endl;
      return 0;
    }
    return 2;
  }
  else return 1;
}
//--------------------------------------------------------------------------------------------




//--------------------------------------------------------------------------------------------
// size_t  binarySearch(std::vector<Double_t> sortedVector, Double_t value)
// {
//   if(sortedVector.size() < 1)
//   {
//     cout << "Error in binarySearch. Vector is empty. Returning UINT_MAX and aborting" << endl;
//     return UINT_MAX;    
//   }
//   
//   //Search only for values inside the vector
//   if (value < sortedVector.at(0) ) return -1;
//   if (value > sortedVector.at(sortedVector.size()-1) ) return -1;
//   
//   size_t  upperIdx = sortedVector.size();
//   size_t  lowerIdx = 0;
//   size_t  middleIdx;
//   while(upperIdx - lowerIdx > 1)
//   {
//     middleIdx = (upperIdx + lowerIdx)/2;
//     if(sortedVector.at(middleIdx) - value > 0 ) upperIdx = middleIdx;
//     else lowerIdx = middleIdx;  
//   }
//   if (lowerIdx == 0) //on lower boundary
//   {
//     
//     if(value < (sortedVector.at(0) + sortedVector.at(1))/2) return 0;
//     else return 1;
//   }
//   if (upperIdx >= sortedVector.size()-1) //on upper boundary
//   {
//     if(value < (sortedVector.at(sortedVector.size()-2) + sortedVector.at(sortedVector.size()-1))/2) return sortedVector.size()-2;
//     else return sortedVector.size()-1;
//   }
//   if( (sortedVector.at(lowerIdx) + sortedVector.at(upperIdx))/2. < value) return upperIdx;
//   else return lowerIdx;
// }
//--------------------------------------------------------------------------------------------





/*
    Save generated canvases in a single folder, sorted by date and time
    Note: Need to create base directory by hand.
    Author: Derek Fujimoto
    Date: Dec 2014
*/

#include<TSystem.h>
#include<TDatime.h>
#include<TString.h>
#include<TPad.h>

/// FUNCTION PROTOTYPES ///
#ifndef FigureDump_cxx
void figureDump(TPad* can, TString fileName = "", TString baseDir="/media/sf_PhD/figureDump/");
#endif // FigureDump_cxx

/// FUNCTION: SAVE FIGURES AS GENERATED ///
// Inputs:
//          can: canvas to save
//          baseDir: location to write files and directories
// Outputs: void, writes to disk

//Author - Derek Fujimoto
#ifndef FigureDump_cxx
#define FigureDump_cxx
void figureDump(TPad* can, TString fileName, TString baseDir)
{
    // Variable Declarations
    TDatime* now = new TDatime();       // the time and date of now
    TString dateDir = TString::Format("%04d-%02d-%02d/",now->GetYear(),now->GetMonth(),now->GetDay());    // name of directory to write the day's figures (yyyy_mm_dd)
    
    if( fileName.EqualTo("") ) fileName = TString::Format("%02dh%02dm%02ds",now->GetHour(),now->GetMinute(),now->GetSecond());// name of file (hh_mm_ss)
     
    TString temp;                       // temp in case of multiple file names
    int i = 0;                          // iteration

    // Check for date directory
    dateDir.Prepend(baseDir);
    gSystem->MakeDirectory(dateDir);

    // Update file path
    fileName.Prepend(dateDir);
    checkIfFileExistsAndCreateBackup(fileName+".gif");
    checkIfFileExistsAndCreateBackup(fileName+".C");
    // If file exists - edit file name
//     if(!gSystem->AccessPathName(fileName+".png") || !gSystem->AccessPathName(fileName+".C"))
//     {
//         temp = fileName.Copy().Append("_%d");
//         do
//         {
//             i++;
//             fileName = TString::Format(temp,i);
//         } while(!gSystem->AccessPathName(fileName+".png") || !gSystem->AccessPathName(fileName+".C"));
//     }

    // Save file
    can->SaveAs(fileName+".gif");
    
    
    can->SaveAs(fileName+".C");

    return;
}
#endif // FigureDump_cxx