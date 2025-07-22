#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "TString.h"
#include "TMath.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <tuple>
#include <string>

using namespace std;

namespace constants {

	//----------------------------------------//

	// Kerberos user name
  
	TString UserID = getenv("USER");

	//----------------------------------------//    

    std::tuple<std::string, int, double> myTuple("World", 456, 2.718);

}

#endif