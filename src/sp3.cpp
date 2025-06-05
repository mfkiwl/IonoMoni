#include "sp3.h"

sp3::sp3()
    : sYear(0), sMonth(0), sDay(0), sHour(0), sMinute(0), sSecond(0) {
    for (int i = 0; i < 35; ++i) {
        std::fill(X[i], X[i] + 3288, 0.0);
        std::fill(Y[i], Y[i] + 3288, 0.0);
        std::fill(Z[i], Z[i] + 3288, 0.0);
    }

    for (int i = 0; i < 49; ++i) {
        std::fill(CX[i], CX[i] + 3288, 0.0);
        std::fill(CY[i], CY[i] + 3288, 0.0);
        std::fill(CZ[i], CZ[i] + 3288, 0.0);
    }

    for (int i = 0; i < 39; ++i) {
        std::fill(EX[i], EX[i] + 3288, 0.0);
        std::fill(EY[i], EY[i] + 3288, 0.0);
        std::fill(EZ[i], EZ[i] + 3288, 0.0);
    }

    for (int i = 0; i < 27; ++i) {
        std::fill(RX[i], RX[i] + 3288, 0.0);
        std::fill(RY[i], RY[i] + 3288, 0.0);
        std::fill(RZ[i], RZ[i] + 3288, 0.0);
    }
}
