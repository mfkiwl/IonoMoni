#include "detect_cycle_slip.h"
#include <cmath>

// Detect cycle slips in GNSS observation data and update arc information accordingly
void detect_cycle_slip(double wlAmb[], double GF[], obs& OBS, int i,
    int ARCS[][3000][2], int arc_min_len, double cs_epoch[][2881]) {

    int j = 1;
    while (j <= ARCS[i][0][0]) {
        int e = ARCS[i][j][0];
        double ave_N[3000] = { 0 };   // Moving average of wide-lane ambiguities
        double sigma2[3000] = { 0 };  // Variance
        double sigma[3000] = { 0 };   // Standard deviation

        // Check the head of the current arc for cycle slips using ambiguity and geometry-free
        while (true) {
            if (e + 1 == ARCS[i][j][1] || e == ARCS[i][j][1]) break;
            double amb0 = wlAmb[e];
            double amb1 = wlAmb[e + 1];
            double amb2 = wlAmb[e + 2];
            double gf0 = GF[e];
            double gf1 = GF[e + 1];
            double gf2 = GF[e + 2];

            // Detect cycle slip at epoch e
            if (fabs(amb0 - amb1) > 1 || fabs(amb2 - amb1) > 1 ||
                fabs(gf0 - gf1) > 1 || fabs(gf2 - gf1) > 1) {
                cs_epoch[i][e] = 1; // Mark as cycle slip
                // Invalidate observation data for this epoch
                OBS.C1[i][e] = OBS.C2[i][e] = OBS.L1[i][e] = OBS.L2[i][e] = 0;
                wlAmb[e] = GF[e] = 0;
                e++;
                ARCS[i][j][0] = e;
            }
            else {
                ARCS[i][j][0] = e;
                break;
            }
        }

        // Remove arc if its length is less than the minimum allowed
        if (ARCS[i][j][1] - ARCS[i][j][0] < arc_min_len) {
            for (int k = ARCS[i][j][0]; k <= ARCS[i][j][1]; ++k) {
                cs_epoch[i][k] = 1;
                OBS.C1[i][k] = OBS.C2[i][k] = OBS.L1[i][k] = OBS.L2[i][k] = 0;
                wlAmb[e] = GF[k] = 0;
            }
            // Shift arcs forward to remove the current short arc
            for (int k = j; k < ARCS[i][0][0]; ++k) {
                ARCS[i][k][0] = ARCS[i][k + 1][0];
                ARCS[i][k][1] = ARCS[i][k + 1][1];
            }
            --ARCS[i][0][0];
            continue;
        }

        // Initialize moving average and variance for the arc
        ave_N[1] = wlAmb[ARCS[i][j][0]];
        sigma2[1] = 0;
        sigma[1] = 0;
        int count = 2;

        // Slide through the arc to detect cycle slips using statistical thresholds
        for (int k = ARCS[i][j][0] + 1; k < ARCS[i][j][1]; ++k) {
            ave_N[count] = ave_N[count - 1] + (wlAmb[k] - ave_N[count - 1]) / count;
            sigma2[count] = sigma2[count - 1] + ((wlAmb[k] - ave_N[count - 1]) * (wlAmb[k] - ave_N[count - 1]) - sigma2[count - 1]) / count;
            sigma[count] = sqrt(sigma2[count]);
            double T = fabs(wlAmb[k + 1] - ave_N[count]);
            double I1 = fabs(GF[k + 1] - GF[k]);

            // If no cycle slip detected, continue
            if (T < 4 * sigma[count] && I1 < 5) {
                count++;
                continue;
            }

            // Otherwise, mark cycle slip at epoch k+1
            cs_epoch[i][k + 1] = 1;

            // Handle arc ending at detected slip or end of arc
            if (k + 1 == ARCS[i][j][1]) {
                if (k + 1 - ARCS[i][j][0] > arc_min_len) {
                    // Mark only the slip epoch as invalid
                    OBS.C1[i][k + 1] = OBS.C2[i][k + 1] = OBS.L1[i][k + 1] = OBS.L2[i][k + 1] = 0;
                    wlAmb[k + 1] = GF[k] = 0;
                    ARCS[i][j][1] = k;
                }
                else {
                    // If resulting arc is too short, invalidate the whole segment and remove arc
                    for (int l = ARCS[i][j][0]; l <= k + 1; ++l) {
                        cs_epoch[i][l] = 1;
                        OBS.C1[i][l] = OBS.C2[i][l] = OBS.L1[i][l] = OBS.L2[i][l] = 0;
                        wlAmb[l] = GF[k] = 0;
                    }
                    // Shift arcs forward to remove short arc
                    for (int x = j; x < ARCS[i][0][0]; ++x) {
                        ARCS[i][x][0] = ARCS[i][x + 1][0];
                        ARCS[i][x][1] = ARCS[i][x + 1][1];
                    }
                    --ARCS[i][0][0];
                    --j;
                }
                break;
            }

            double I2 = fabs(GF[k + 2] - GF[k + 1]);
            // Case 1: No big jump after slip, split the arc
            if (fabs(wlAmb[k + 2] - wlAmb[k + 1]) < 1 && I2 < 1) {
                if (k + 1 - ARCS[i][j][0] > arc_min_len) {
                    cs_epoch[i][k + 1] = 1;
                    // Shift later arcs and split current arc
                    for (int l = ARCS[i][0][0]; l >= j + 1; --l) {
                        ARCS[i][l + 1][0] = ARCS[i][l][0];
                        ARCS[i][l + 1][1] = ARCS[i][l][1];
                    }
                    ++ARCS[i][0][0];
                    ARCS[i][j + 1][0] = k + 1;
                    ARCS[i][j + 1][1] = ARCS[i][j][1];
                    ARCS[i][j][1] = k;
                }
                else {
                    // Invalidate the whole short arc
                    for (int l = ARCS[i][j][0]; l <= k; ++l) {
                        cs_epoch[i][l] = 1;
                        OBS.C1[i][l] = OBS.C2[i][l] = OBS.L1[i][l] = OBS.L2[i][l] = 0;
                        wlAmb[l] = GF[k] = 0;
                    }
                    ARCS[i][j][0] = k + 1;
                    --j;
                }
            }
            // Case 2: Detected jump, split arc after k+1
            else {
                if (k + 1 - ARCS[i][j][0] > arc_min_len) {
                    cs_epoch[i][k + 1] = 1;
                    OBS.C1[i][k + 1] = OBS.C2[i][k + 1] = OBS.L1[i][k + 1] = OBS.L2[i][k + 1] = 0;
                    wlAmb[k + 1] = GF[k] = 0;
                    // Shift and split arcs
                    for (int l = ARCS[i][0][0]; l >= j + 1; --l) {
                        ARCS[i][l + 1][0] = ARCS[i][l][0];
                        ARCS[i][l + 1][1] = ARCS[i][l][1];
                    }
                    ++ARCS[i][0][0];
                    ARCS[i][j + 1][0] = k + 2;
                    ARCS[i][j + 1][1] = ARCS[i][j][1];
                    ARCS[i][j][1] = k;
                }
                else {
                    // Invalidate the whole short arc
                    for (int l = ARCS[i][j][0]; l <= k + 1; ++l) {
                        cs_epoch[i][l] = 1;
                        OBS.C1[i][l] = OBS.C2[i][l] = OBS.L1[i][l] = OBS.L2[i][l] = 0;
                        wlAmb[l] = GF[k] = 0;
                    }
                    ARCS[i][j][0] = k + 2;
                    --j;
                }
            }

            break;
        }

        ++j;
    }
}
