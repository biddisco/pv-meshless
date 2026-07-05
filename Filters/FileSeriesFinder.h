/*=========================================================================

  Project                 : XdmfGenerator
  Module                  : XdmfFileSeriesFinder.h

  Authors:
     John Biddiscombe     Jerome Soumagne
     biddisco@cscs.ch     soumagne@cscs.ch

  Copyright (C) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing
  1) This copyright notice appears on all copies of source code
  2) An acknowledgment appears with any substantial usage of the code
  3) If this code is contributed to any other open source project, it
  must not be reformatted such that the indentation, bracketing or
  overall style is modified significantly.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  This work has received funding from the European Community's Seventh
  Framework Programme (FP7/2007-2013) under grant agreement 225967 “NextMuSE”

=========================================================================*/
#ifndef FILESERIESFINDER_H
#define FILESERIESFINDER_H

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <vtksys/SystemTools.hxx>
#include <vtksys/Glob.hxx>
#include <vtksys/RegularExpression.hxx>
#include <vtksys/Process.h>
#include <sstream>
#include "vtkSystemIncludes.h"

typedef std::vector<std::string> stringlist;
typedef std::pair< std::string, std::string > maptype;
//
class VTK_EXPORT FileSeriesFinder {
  public:
    FileSeriesFinder(std::string filenamepattern);
    FileSeriesFinder();

    //
    // If loading just one file and no scan is required, use this function
    //
    void SetFileName(const char *filename);

    //
    // These should be called before calling Scan
    //
    void SetPrefixRegEx(const char *prefix_regex);
    void SetTimeRegEx(const char *time_regex);
    void SetBlockRegEx(const char *block_regex);
    void SetBlockSubDirRegEx(const char *block_regex);
    void SetVarRegEx(const char *var_regex);
    void SetText0RegEx(const char *text0_regex);
    void SetExtRegEx(const char *ext_regex);
    //
    // Scan will search for files
    void Scan(const char *inputfile);
    // print out all the files found
    void TestFilenameGeneration();
    // return the filename for a given time.block,var
    std::string GenerateFileName(int T, int B, int V);
    // assume block/var are zero and just use time
    std::string GenerateFileName(int T);
    //
    int GetNumberOfTimeSteps();
    int GetNumberOfBlocks();
    int GetNumberOfVars();

    void GetTimeValues(std::vector<double> &values);
    double GetTimeValue(int index);

  protected:
    std::string GenerateNumericPattern(std::string match);
    std::string GenerateGlobString(unsigned int index, stringlist &patterns, bool regexmode);
    void FindGlobbedSegments(std::string &pattern, stringlist &files, bool numeric);


  public:
    std::string FileNamePattern; 
    std::string PrefixRegEx;
    std::string BlockRegEx;
    std::string BlockSubDirRegEx;
    std::string VarRegEx; 
    std::string ExtRegEx;
    std::string TimeRegEx;
    std::string Text0RegEx;
    //
    std::map<std::string, std::string> RegExMap;
    stringlist patterns, regexmatches, Tstrings, Bstrings, Vstrings;
    std::string regex;
    int TimeIndex, NumberOfTimeSteps;
    int BlockIndex, NumberOfBlocks;
    int VarIndex, NumberOfVars;
    std::string timeform, blockform, varform, filepattern;
};

#endif // FILESERIESFINDER_H
