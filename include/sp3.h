#pragma once
#include <algorithm>

struct sp3 {
    double sYear, sMonth, sDay, sHour, sMinute, sSecond;
    double X[35][3288], Y[35][3288], Z[35][3288];
    double CX[49][3288], CY[49][3288], CZ[49][3288];
    double EX[39][3288], EY[39][3288], EZ[39][3288];
    double RX[27][3288], RY[27][3288], RZ[27][3288];

    sp3(); 
};