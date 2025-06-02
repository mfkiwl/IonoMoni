#include "XYZtoBLH.h"
#include <cmath>
void XYZtoBLH(double sx, double sy, double sz, double& sb, double& sl) {
    double p = sqrt(sx * sx + sy * sy);
    double theta = atan2(sz * a, p * (1 - f));

    sb = atan2(sz + e2 * a * sin(theta) * sin(theta) * sin(theta), p - e2 * a * cos(theta) * cos(theta) * cos(theta));
    sl = atan2(sy, sx);
}
