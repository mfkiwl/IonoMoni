#pragma once
#include <string>
#include <algorithm>

struct obs {
    double X, Y, Z;
    double C1[60][2888];
    double C2[60][2888];
    double L1[60][2888];
    double L2[60][2888];
    //double S1[60][2888];   
   // double S2[60][2888];  

    obs();
    void reset();
};
