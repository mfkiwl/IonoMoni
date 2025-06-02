#include "get_AATR.h"
#include <iostream>
#include <vector>

#include <cstring>
#include "read_DCB.h"
#include "get_arc.h"
#include "get_Elevation.h"
#include "constants.h"

#include <fstream>   
#include <iomanip>  
#include <filesystem> 
#include <cmath>     
#include "gcfg_ppp.h"

void AATR_G_prepro(obs& OBS, const std::string& stationName, sp3 SP3[], const std::string& txt_output_path, t_gcfg_ppp& gset){
    double f1 = 1575.42e6;
    double f2 = 1227.6e6;
    double lambda_wl = 299792458 / (f1 - f2);
    double MW[33][2881], GF[33][2881], wlAmb[33][2881];
    double GPSP4[35][3000], GPS_smoothedP4[35][3000];
    double GPSL4[35][3000];
    double AATR_LGF[35][3000];
    double AATR_ROT[33][2888];
    double cs_epoch[33][2881];
    int arc_min_len = gset.arc_min_length();
    int interval_size = gset.aatr_interval();  
    memset(AATR_ROT, 0, sizeof(AATR_ROT));
    memset(cs_epoch, 0, sizeof(cs_epoch));
    for (int i = 1; i <= 32; i++) {
        for (int j = 1; j <= 2880; j++) {
            MW[i][j] = lambda_wl * (OBS.L1[i][j] - OBS.L2[i][j]) - (f1 * OBS.C1[i][j] + f2 * OBS.C2[i][j]) / (f1 + f2);
            GF[i][j] = OBS.L1[i][j] - f1 * OBS.L2[i][j] / f2;
            wlAmb[i][j] = MW[i][j] * -1.0;
        }
    }

    for (int i = 1; i <= 32; i++) {
        int len = 2881;
        int num_arcs = 0;
        std::vector<std::vector<int>> arcs;

        Get_arc(MW[i], len, arcs, num_arcs);

        for (int k = 0; k < num_arcs; ++k) {
            if (arcs[k][1] - arcs[k][0] < arc_min_len) {
                for (int e = arcs[k][0]; e <= arcs[k][1]; ++e) {
                    OBS.C1[i][e] = 0;
                    OBS.C2[i][e] = 0;
                    OBS.L1[i][e] = 0;
                    OBS.L2[i][e] = 0;
                }
                arcs.erase(arcs.begin() + k);
                --k;
                --num_arcs;
            }
        }

        int ARCS[33][3000][2];
        ARCS[i][0][0] = num_arcs;

        for (int k = 0; k < num_arcs; ++k) {
            ARCS[i][k + 1][0] = arcs[k][0];
            ARCS[i][k + 1][1] = arcs[k][1];
        }

        int j = 1;
        while (j <= ARCS[i][0][0]) {
            int e = ARCS[i][j][0];
            while (true) {
                if (e + 1 == ARCS[i][j][1] || e == ARCS[i][j][1]) {
                    break;
                }
                double amb0 = wlAmb[i][e];
                double amb1 = wlAmb[i][e + 1];
                double amb2 = wlAmb[i][e + 2];
                double gf0 = GF[i][e];
                double gf1 = GF[i][e + 1];
                double gf2 = GF[i][e + 2];
                double diff01 = std::abs(amb0 - amb1);
                double diff21 = std::abs(amb2 - amb1);
                double gfDiff01 = std::abs(gf0 - gf1);
                double gfDiff21 = std::abs(gf2 - gf1);
                if (diff01 > 1 || diff21 > 1 || gfDiff01 > 1 || gfDiff21 > 1) {
                    MW[i][e] = 0;
                    OBS.C1[i][e] = 0;
                    OBS.C2[i][e] = 0;
                    OBS.L1[i][e] = 0;
                    OBS.L2[i][e] = 0;
                    wlAmb[i][e] = 0;
                    GF[i][e] = 0;
                    e++;
                    ARCS[i][j][0] = e;
                }
                else {
                    ARCS[i][j][0] = e;
                    break;
                }
            }

            if (ARCS[i][j][1] - ARCS[i][j][0] < arc_min_len) {
                for (int k = ARCS[i][j][0]; k <= ARCS[i][j][1]; k++) {
                    OBS.C1[i][k] = 0;
                    OBS.C2[i][k] = 0;
                    OBS.L1[i][k] = 0;
                    OBS.L2[i][k] = 0;
                    MW[i][k] = 0;
                    wlAmb[i][e] = 0;
                    GF[i][k] = 0;
                }
                for (int k = j; k < ARCS[i][0][0]; k++) {
                    ARCS[i][k][0] = ARCS[i][k + 1][0];
                    ARCS[i][k][1] = ARCS[i][k + 1][1];
                }
                ARCS[i][0][0]--;
            }
            double ave_N[3000];
            double sigma2[3000];
            double sigma[3000];
            ave_N[1] = wlAmb[i][ARCS[i][j][0]];
            sigma2[1] = 0;
            sigma[1] = 0;
            int count = 2;
            for (int k = ARCS[i][j][0] + 1; k < ARCS[i][j][1]; k++) {

                ave_N[count] = ave_N[count - 1] + (wlAmb[i][k] - ave_N[count - 1]) / count;
                sigma2[count] = sigma2[count - 1] + ((wlAmb[i][k] - ave_N[count - 1]) * (wlAmb[i][k] - ave_N[count - 1]) - sigma2[count - 1]) / count;
                sigma[count] = sqrt(sigma2[count]);
                double T = std::abs(wlAmb[i][k + 1] - ave_N[count]);
                double I1 = std::abs(GF[i][k + 1] - GF[i][k]);
                if (T < 4 * sigma[count] && I1 < 5) {
                    count++;
                    continue;

                }
                else {
                    if (k + 1 == ARCS[i][j][1]) {
                        if (k + 1 - ARCS[i][j][0] > arc_min_len) {
                            MW[i][k + 1] = 0;
                            OBS.C1[i][k + 1] = 0;
                            OBS.C2[i][k + 1] = 0;
                            OBS.L1[i][k + 1] = 0;
                            OBS.L2[i][k + 1] = 0;
                            wlAmb[i][k + 1] = 0;
                            GF[i][k] = 0;
                            ARCS[i][j][1] = k;

                        }
                        else {
                            for (int l = ARCS[i][j][0]; l <= k + 1; l++) {
                                MW[i][l] = 0;
                                OBS.C1[i][l] = 0;
                                OBS.C2[i][l] = 0;
                                OBS.L1[i][l] = 0;
                                OBS.L2[i][l] = 0;
                                wlAmb[i][l] = 0;
                                GF[i][k] = 0;

                            }
                            for (int x = j; x < ARCS[i][0][0]; x++) {
                                ARCS[i][x][0] = ARCS[i][x + 1][0];
                                ARCS[i][x][1] = ARCS[i][x + 1][1];
                            }
                            j--;
                            ARCS[i][0][0]--;
                        }
                        break;
                    }
                    double I2 = std::abs(GF[i][k + 2] - GF[i][k + 1]);
                    if (std::abs(wlAmb[i][k + 2] - wlAmb[i][k + 1]) < 1 && I2 < 1) {
                        if (k + 1 - ARCS[i][j][0] > arc_min_len) {
                            cs_epoch[i][k + 1] = 1;
                            for (int l = ARCS[i][0][0]; l >= j + 1; l--) {
                                ARCS[i][l + 1][0] = ARCS[i][l][0];
                                ARCS[i][l + 1][1] = ARCS[i][l][1];

                            }
                            ARCS[i][0][0]++;
                            ARCS[i][j + 1][0] = k + 1;
                            ARCS[i][j + 1][1] = ARCS[i][j][1];
                            ARCS[i][j][1] = k;


                        }
                        else {
                            for (int l = ARCS[i][j][0]; l <= k; l++) {
                                MW[i][l] = 0;
                                OBS.C1[i][l] = 0;
                                OBS.C2[i][l] = 0;
                                OBS.L1[i][l] = 0;
                                OBS.L2[i][l] = 0;
                                wlAmb[i][l] = 0;
                                GF[i][k] = 0;
                            }
                            ARCS[i][j][0] = k + 1;
                            j--;

                        }
                    }
                    else {
                        if (k + 1 - ARCS[i][j][0] > arc_min_len) {
                            MW[i][k + 1] = 0;
                            OBS.C1[i][k + 1] = 0;
                            OBS.C2[i][k + 1] = 0;
                            OBS.L1[i][k + 1] = 0;
                            OBS.L2[i][k + 1] = 0;
                            wlAmb[i][k + 1] = 0;
                            GF[i][k] = 0;
                            for (int l = ARCS[i][0][0]; l >= j + 1; l--) {
                                ARCS[i][l + 1][0] = ARCS[i][l][0];
                                ARCS[i][l + 1][1] = ARCS[i][l][1];

                            }
                            ARCS[i][0][0]++;
                            ARCS[i][j + 1][0] = k + 2;
                            ARCS[i][j + 1][1] = ARCS[i][j][1];
                            ARCS[i][j][1] = k;

                        }
                        else {
                            for (int l = ARCS[i][j][0]; l <= k + 1; l++) {
                                MW[i][l] = 0;
                                OBS.C1[i][l] = 0;
                                OBS.C2[i][l] = 0;
                                OBS.L1[i][l] = 0;
                                OBS.L2[i][l] = 0;
                                wlAmb[i][l] = 0;
                                GF[i][k] = 0;
                            }
                            ARCS[i][j][0] = k + 2;
                            j--;
                        }
                    }

                    break;
                }

            }

            j++;
        }
        // Calculate AATR_ROT
        for (int j = 1; j <= 2880; j++) {
            AATR_LGF[i][j] = (c / f1) * OBS.L1[i][j] - (c / f2) * OBS.L2[i][j];

        }

        std::string rot_unit = gset.rot_unit();
        std::transform(rot_unit.begin(), rot_unit.end(), rot_unit.begin(), ::tolower);  
        double rotNumber;

        if (rot_unit == "sec") {
            rotNumber = 30.0 * 1e16 * IONO_COEFF * (1 / (f1 * f1) - 1 / (f2 * f2));
        }
        else {
            rotNumber = 0.5 * 1e16 * IONO_COEFF * (1 / (f1 * f1) - 1 / (f2 * f2));
        }

        for (int j = 2; j <= 2880; j++) {
            if (AATR_LGF[i][j] != 0 && AATR_LGF[i][j - 1] != 0 && cs_epoch[i][j] != 1) {
                AATR_ROT[i][j] = (AATR_LGF[i][j] - AATR_LGF[i][j - 1]) / rotNumber;


            }
        }
        
    }
    
    double AATR_block[25] = { 0 };   // Store the computed result for each epoch block
    int block_index = 0;
    for (int j = 1; j <= 2880; j += interval_size) { 
        double sum_block = 0;
        int count = 0;

        for (int k = j; k < j + interval_size; k++) {
            for (int i = 1; i <= 32; i++) {
                if (AATR_ROT[i][k] != 0) { 
                    double E, A;
                    Get_Elevation(SP3[1].X[i][k] * 1000, SP3[1].Y[i][k] * 1000, SP3[1].Z[i][k] * 1000,
                        OBS.X, OBS.Y, OBS.Z, E, A);

                    double E_rad = E * PI / 180.0;
                    double cos_E = cos(E_rad);
                    double tmp = (Re * cos_E) / (Re + h);
                    double M = 1.0 / sqrt(1.0 - tmp * tmp);
                    double AATR = (1.0 / (M * M)) * AATR_ROT[i][k];

                    sum_block += AATR * AATR;
                    count++;
                }
            }
        }

        if (count > 0)
            AATR_block[block_index] = sqrt(sum_block / count);  
        else
            AATR_block[block_index] = 0;

        block_index++;
    }


    std::string aatr_out_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GPS_AATR.txt";
    std::filesystem::create_directories(std::filesystem::path(aatr_out_path).parent_path());

    std::ofstream out(aatr_out_path);
    if (!out.is_open()) {
        std::cerr << "Failed to open GPS AATR output file" << std::endl;
        return;
    }
    out << std::fixed << std::setprecision(5);
    out << std::setw(6) << "Block" << std::setw(15) << "AATR" << "\n";

    for (int k = 0; k < block_index; ++k) {
        char block_buf[20];
        sprintf(block_buf, "Block %02d:", k + 1);
        double val = 1624.7244 * AATR_block[k];  //  1624.7244 * AATR_block[k]
        out << std::setw(10) << block_buf << std::setw(12) << val << "\n";
    }

    out.close();


}

void AATR_R_prepro(obs& OBS, const std::string& stationName, sp3 SP3[], const std::string& txt_output_path, t_gcfg_ppp& gset){

    double AATR_LGF[27][3000];
    double AATR_ROT[25][2888];
    double cs_epoch[25][2881];
    memset(AATR_ROT, 0, sizeof(AATR_ROT));
    memset(cs_epoch, 0, sizeof(cs_epoch));
    int arc_min_len = gset.arc_min_length();
    int interval_size = gset.aatr_interval();
    // Frequency channel numbers for the 24 GLONASS satellites
    int Fre[25] = { 99, 1, -4, 5, 6, 1, -4, 5, 6, -2, -7, 0, -1, -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2 };
    double f1[25], f2[25], lambda_wl[25];
    double MW[25][2881], GF[25][2881], wlAmb[25][2881];

    for (int i = 1; i <= 24; i++) {
        f1[i] = (1602.0 + Fre[i] * 0.5625) * 1e6;  
        f2[i] = (1246.0 + Fre[i] * 0.4375) * 1e6;  
        lambda_wl[i] = c / (f1[i] - f2[i]);  

        for (int j = 1; j <= 2880; j++) {
            MW[i][j] = lambda_wl[i] * (OBS.L1[i][j] - OBS.L2[i][j]) - (f1[i] * OBS.C1[i][j] + f2[i] * OBS.C2[i][j]) / (f1[i] + f2[i]);
            GF[i][j] = OBS.L1[i][j] - f1[i] * OBS.L2[i][j] / f2[i];
            wlAmb[i][j] = MW[i][j] * -1.0;
        }

        int len = 2881;         
        int num_arcs = 0;     
        std::vector<std::vector<int>> arcs;  
        Get_arc(MW[i], len, arcs, num_arcs);  
        for (int k = 0; k < num_arcs; ++k) {
            if (arcs[k][1] - arcs[k][0] < arc_min_len) {

                for (int e = arcs[k][0]; e <= arcs[k][1]; ++e) {

                    OBS.C1[i][e] = 0;
                    OBS.C2[i][e] = 0;
                    OBS.L1[i][e] = 0;
                    OBS.L2[i][e] = 0; 
                }
                arcs.erase(arcs.begin() + k);
                --k;  
                --num_arcs;  
            }
        }


        int ARCS[25][3000][2];

        ARCS[i][0][0] = num_arcs;

        for (int k = 0; k < num_arcs; ++k) {
            ARCS[i][k + 1][0] = arcs[k][0];
            ARCS[i][k + 1][1] = arcs[k][1];
        }

        int j = 1;
        while (j <= ARCS[i][0][0]) {

            int e = ARCS[i][j][0];
            while (true) {
                if (e + 1 == ARCS[i][j][1] || e == ARCS[i][j][1]) {
                    break;

                }
                double amb0 = wlAmb[i][e];
                double amb1 = wlAmb[i][e + 1];
                double amb2 = wlAmb[i][e + 2];
                double gf0 = GF[i][e];
                double gf1 = GF[i][e + 1];
                double gf2 = GF[i][e + 2];
                double diff01 = std::abs(amb0 - amb1);
                double diff21 = std::abs(amb2 - amb1);
                double gfDiff01 = std::abs(gf0 - gf1);
                double gfDiff21 = std::abs(gf2 - gf1);
                if (diff01 > 1 || diff21 > 1 || gfDiff01 > 1 || gfDiff21 > 1) {
                  
                    MW[i][e] = 0;
                    OBS.C1[i][e] = 0;
                    OBS.C2[i][e] = 0;
                    OBS.L1[i][e] = 0;
                    OBS.L2[i][e] = 0;
                    wlAmb[i][e] = 0;
                    GF[i][e] = 0;
                    e++;
                    ARCS[i][j][0] = e;
                }
                else {
                    ARCS[i][j][0] = e;
                    break;

                }
            }
            if (ARCS[i][j][1] - ARCS[i][j][0] < arc_min_len) {
                for (int k = ARCS[i][j][0]; k <= ARCS[i][j][1]; k++) {
                    OBS.C1[i][k] = 0;
                    OBS.C2[i][k] = 0;
                    OBS.L1[i][k] = 0;
                    OBS.L2[i][k] = 0;
                    MW[i][k] = 0;
                    wlAmb[i][e] = 0;
                    GF[i][k] = 0;

                }
                for (int k = j; k < ARCS[i][0][0]; k++) {
                    ARCS[i][k][0] = ARCS[i][k + 1][0];
                    ARCS[i][k][1] = ARCS[i][k + 1][1];
                }
                ARCS[i][0][0]--;
            }

            double ave_N[3000];
            double sigma2[3000];
            double sigma[3000];
            ave_N[1] = wlAmb[i][ARCS[i][j][0]];
            sigma2[1] = 0;
            sigma[1] = 0;
            int count = 2;
            for (int k = ARCS[i][j][0] + 1; k < ARCS[i][j][1]; k++) {

                ave_N[count] = ave_N[count - 1] + (wlAmb[i][k] - ave_N[count - 1]) / count;
                sigma2[count] = sigma2[count - 1] + ((wlAmb[i][k] - ave_N[count - 1]) * (wlAmb[i][k] - ave_N[count - 1]) - sigma2[count - 1]) / count;
                sigma[count] = sqrt(sigma2[count]);
                double T = std::abs(wlAmb[i][k + 1] - ave_N[count]);
                double I1 = std::abs(GF[i][k + 1] - GF[i][k]);
                if (T < 4 * sigma[count] && I1 < 0.28) {
                    count++;
                    continue;

                }
                else {
                    if (k + 1 == ARCS[i][j][1]) {

                        if (k + 1 - ARCS[i][j][0] > 1) {
                            MW[i][k + 1] = 0;
                            OBS.C1[i][k + 1] = 0;
                            OBS.C2[i][k + 1] = 0;
                            OBS.L1[i][k + 1] = 0;
                            OBS.L2[i][k + 1] = 0;
                            wlAmb[i][k + 1] = 0;
                            GF[i][k] = 0;
                            ARCS[i][j][1] = k;

                        }
                        else {
                            for (int l = ARCS[i][j][0]; l <= k + 1; l++) {
                                MW[i][l] = 0;
                                OBS.C1[i][l] = 0;
                                OBS.C2[i][l] = 0;
                                OBS.L1[i][l] = 0;
                                OBS.L2[i][l] = 0;
                                wlAmb[i][l] = 0;
                                GF[i][k] = 0;

                            }
                            for (int x = j; x < ARCS[i][0][0]; x++) {
                                ARCS[i][x][0] = ARCS[i][x + 1][0];
                                ARCS[i][x][1] = ARCS[i][x + 1][1];
                            }
                            j--;
                            ARCS[i][0][0]--;
                        }
                        break;
                    }
                    double I2 = std::abs(GF[i][k + 2] - GF[i][k + 1]);
                    if (std::abs(wlAmb[i][k + 2] - wlAmb[i][k + 1]) < 1 && I2 < 1) {
                        if (k + 1 - ARCS[i][j][0] > 1) {
                            cs_epoch[i][k + 1] = 1;
                            for (int l = ARCS[i][0][0]; l >= j + 1; l--) {
                                ARCS[i][l + 1][0] = ARCS[i][l][0];
                                ARCS[i][l + 1][1] = ARCS[i][l][1];

                            }
                            ARCS[i][0][0]++;
                            ARCS[i][j + 1][0] = k + 1;
                            ARCS[i][j + 1][1] = ARCS[i][j][1];
                            ARCS[i][j][1] = k;
                            
                        }
                        else {
                            for (int l = ARCS[i][j][0]; l <= k; l++) {
                                MW[i][l] = 0;
                                OBS.C1[i][l] = 0;
                                OBS.C2[i][l] = 0;
                                OBS.L1[i][l] = 0;
                                OBS.L2[i][l] = 0;
                                wlAmb[i][l] = 0;
                                GF[i][k] = 0;
                            }
                            ARCS[i][j][0] = k + 1;
                            j--;

                        }
                    }
                    else {
                        if (k + 1 - ARCS[i][j][0] > 1) {
                            MW[i][k + 1] = 0;
                            OBS.C1[i][k + 1] = 0;
                            OBS.C2[i][k + 1] = 0;
                            OBS.L1[i][k + 1] = 0;
                            OBS.L2[i][k + 1] = 0;
                            wlAmb[i][k + 1] = 0;
                            GF[i][k] = 0;
                            for (int l = ARCS[i][0][0]; l >= j + 1; l--) {
                                ARCS[i][l + 1][0] = ARCS[i][l][0];
                                ARCS[i][l + 1][1] = ARCS[i][l][1];

                            }
                            ARCS[i][0][0]++;
                            ARCS[i][j + 1][0] = k + 2;
                            ARCS[i][j + 1][1] = ARCS[i][j][1];
                            ARCS[i][j][1] = k;

                        }
                        else {
                            for (int l = ARCS[i][j][0]; l <= k + 1; l++) {
                                MW[i][l] = 0;
                                OBS.C1[i][l] = 0;
                                OBS.C2[i][l] = 0;
                                OBS.L1[i][l] = 0;
                                OBS.L2[i][l] = 0;
                                wlAmb[i][l] = 0;
                                GF[i][k] = 0;
                            }
                            ARCS[i][j][0] = k + 2;
                            j--;
                        }
                    }
                    
                    break;
                }

            }

            j++;
            
        }

        for (int j = 1; j <= 2880; j++) {
            AATR_LGF[i][j] = (c / f1[i]) * OBS.L1[i][j] - (c / f2[i]) * OBS.L2[i][j];

        }

        std::string rot_unit = gset.rot_unit();
        std::transform(rot_unit.begin(), rot_unit.end(), rot_unit.begin(), ::tolower);  
        double rotNumber;

        if (rot_unit == "sec") {
            rotNumber = 30.0 * 1e16 * IONO_COEFF * (1 / (f1[i] * f1[i]) - 1 / (f2[i] * f2[i]));
        }
        else {

            rotNumber = 0.5 * 1e16 * IONO_COEFF * (1 / (f1[i] * f1[i]) - 1 / (f2[i] * f2[i]));
        }


        for (int j = 2; j <= 2880; j++) {
            if (AATR_LGF[i][j] != 0 && AATR_LGF[i][j - 1] != 0 && cs_epoch[i][j] != 1) {
                AATR_ROT[i][j] = (AATR_LGF[i][j] - AATR_LGF[i][j - 1]) / rotNumber;


            }
        }

    }

    double AATR_block[25] = { 0 };   
    int block_index = 0;
    for (int j = 1; j <= 2880; j += interval_size) { 
        double sum_block = 0;
        int count = 0;

        for (int k = j; k < j + interval_size; k++) {
            for (int i = 1; i <= 24; i++) {
                if (AATR_ROT[i][k] != 0) { 
                    double E, A;
                    Get_Elevation(SP3[1].X[i][k] * 1000, SP3[1].Y[i][k] * 1000, SP3[1].Z[i][k] * 1000,
                        OBS.X, OBS.Y, OBS.Z, E, A);

                    double E_rad = E * PI / 180.0;
                    double cos_E = cos(E_rad);
                    double tmp = (Re * cos_E) / (Re + h);
                    double M = 1.0 / sqrt(1.0 - tmp * tmp);
                    double AATR = (1.0 / (M * M)) * AATR_ROT[i][k];

                    sum_block += AATR * AATR;
                    count++;
                }
            }
        }

        if (count > 0)
            AATR_block[block_index] = sqrt(sum_block / count);  
        else
            AATR_block[block_index] = 0;

        block_index++;
    }

    std::string aatr_out_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GLO_AATR.txt";
    std::filesystem::create_directories(std::filesystem::path(aatr_out_path).parent_path());

    std::ofstream out(aatr_out_path);
    if (!out.is_open()) {
        std::cerr << "Failed to open GLO AATR output file." << std::endl;
        return;
    }

    out << std::fixed << std::setprecision(5);

    out << std::setw(6) << "Block" << std::setw(15) << "AATR" << "\n";

    for (int k = 0; k < block_index; ++k) {
        char block_buf[20];
        sprintf(block_buf, "Block %02d:", k + 1);
        double val = AATR_block[k];  
        out << std::setw(10) << block_buf << std::setw(12) << val << "\n";
    }

    out.close();


}
void AATR_E_prepro(obs& OBS, const std::string& stationName, sp3 SP3[], const std::string& txt_output_path, t_gcfg_ppp& gset){
    double f1 = 1575.42e6;
    double f2 = 1176.45e6;
    double lambda_wl = 299792458 / (f1 - f2);

    double MW[37][2881], GF[37][2881], wlAmb[37][2881];


    double AATR_LGF[39][3000];
    double AATR_ROT[37][2888];

    double cs_epoch[37][2881];
    int arc_min_len = gset.arc_min_length();
    int interval_size = gset.aatr_interval();
    memset(AATR_ROT, 0, sizeof(AATR_ROT));
    memset(cs_epoch, 0, sizeof(cs_epoch));
    for (int i = 1; i <= 36; i++) {
        for (int j = 1; j <= 2880; j++) {
            MW[i][j] = lambda_wl * (OBS.L1[i][j] - OBS.L2[i][j]) - (f1 * OBS.C1[i][j] + f2 * OBS.C2[i][j]) / (f1 + f2);
            GF[i][j] = OBS.L1[i][j] - f1 * OBS.L2[i][j] / f2;
            wlAmb[i][j] = MW[i][j] * -1.0;
        }
    }

    for (int i = 1; i <= 36; i++) {

        int len = 2881;         
        int num_arcs = 0;      
        std::vector<std::vector<int>> arcs; 
        Get_arc(MW[i], len, arcs, num_arcs);  
        for (int k = 0; k < num_arcs; ++k) {
            if (arcs[k][1] - arcs[k][0] < arc_min_len) {

                for (int e = arcs[k][0]; e <= arcs[k][1]; ++e) {
                    OBS.C1[i][e] = 0;
                    OBS.C2[i][e] = 0;
                    OBS.L1[i][e] = 0;
                    OBS.L2[i][e] = 0;
                }
                arcs.erase(arcs.begin() + k);
                --k;
                --num_arcs;
            }
        }

        int ARCS[37][3000][2];

        ARCS[i][0][0] = num_arcs;

        for (int k = 0; k < num_arcs; ++k) {
            ARCS[i][k + 1][0] = arcs[k][0];
            ARCS[i][k + 1][1] = arcs[k][1];
        }

        int j = 1;
        while (j <= ARCS[i][0][0]) {

            int e = ARCS[i][j][0];
            while (true) {
                if (e + 1 == ARCS[i][j][1] || e == ARCS[i][j][1]) {
                    break;

                }
                double amb0 = wlAmb[i][e];
                double amb1 = wlAmb[i][e + 1];
                double amb2 = wlAmb[i][e + 2];
                double gf0 = GF[i][e];
                double gf1 = GF[i][e + 1];
                double gf2 = GF[i][e + 2];
                double diff01 = std::abs(amb0 - amb1);
                double diff21 = std::abs(amb2 - amb1);
                double gfDiff01 = std::abs(gf0 - gf1);
                double gfDiff21 = std::abs(gf2 - gf1);
                if (diff01 > 1 || diff21 > 1 || gfDiff01 > 1 || gfDiff21 > 1) {
                    MW[i][e] = 0;
                    OBS.C1[i][e] = 0;
                    OBS.C2[i][e] = 0;
                    OBS.L1[i][e] = 0;
                    OBS.L2[i][e] = 0;
                    wlAmb[i][e] = 0;
                    GF[i][e] = 0;
                    e++;
                    ARCS[i][j][0] = e;
                }
                else {
                    ARCS[i][j][0] = e;
                    break;

                }
            }
            if (ARCS[i][j][1] - ARCS[i][j][0] < arc_min_len) {
                for (int k = ARCS[i][j][0]; k <= ARCS[i][j][1]; k++) {
                    OBS.C1[i][k] = 0;
                    OBS.C2[i][k] = 0;
                    OBS.L1[i][k] = 0;
                    OBS.L2[i][k] = 0;
                    MW[i][k] = 0;
                    wlAmb[i][e] = 0;
                    GF[i][k] = 0;

                }
                for (int k = j; k < ARCS[i][0][0]; k++) {
                    ARCS[i][k][0] = ARCS[i][k + 1][0];
                    ARCS[i][k][1] = ARCS[i][k + 1][1];
                }
                ARCS[i][0][0]--;
            }

            double ave_N[3000];
            double sigma2[3000];
            double sigma[3000];
            ave_N[1] = wlAmb[i][ARCS[i][j][0]];
            sigma2[1] = 0;
            sigma[1] = 0;
            int count = 2;
            for (int k = ARCS[i][j][0] + 1; k < ARCS[i][j][1]; k++) {

                ave_N[count] = ave_N[count - 1] + (wlAmb[i][k] - ave_N[count - 1]) / count;
                sigma2[count] = sigma2[count - 1] + ((wlAmb[i][k] - ave_N[count - 1]) * (wlAmb[i][k] - ave_N[count - 1]) - sigma2[count - 1]) / count;
                sigma[count] = sqrt(sigma2[count]);
                double T = std::abs(wlAmb[i][k + 1] - ave_N[count]);
                double I1 = std::abs(GF[i][k + 1] - GF[i][k]);
                if (T < 4 * sigma[count] && I1 < 0.28) {
                    count++;
                    continue;

                }
                else {
                    if (k + 1 == ARCS[i][j][1]) {

                        if (k + 1 - ARCS[i][j][0] > 1) {
                            MW[i][k + 1] = 0;
                            OBS.C1[i][k + 1] = 0;
                            OBS.C2[i][k + 1] = 0;
                            OBS.L1[i][k + 1] = 0;
                            OBS.L2[i][k + 1] = 0;
                            wlAmb[i][k + 1] = 0;
                            GF[i][k] = 0;
                            ARCS[i][j][1] = k;

                        }
                        else {
                            for (int l = ARCS[i][j][0]; l <= k + 1; l++) {
                                MW[i][l] = 0;
                                OBS.C1[i][l] = 0;
                                OBS.C2[i][l] = 0;
                                OBS.L1[i][l] = 0;
                                OBS.L2[i][l] = 0;
                                wlAmb[i][l] = 0;
                                GF[i][k] = 0;

                            }
                            for (int x = j; x < ARCS[i][0][0]; x++) {
                                ARCS[i][x][0] = ARCS[i][x + 1][0];
                                ARCS[i][x][1] = ARCS[i][x + 1][1];
                            }
                            j--;
                            ARCS[i][0][0]--;
                        }
                        break;
                    }
                    double I2 = std::abs(GF[i][k + 2] - GF[i][k + 1]);
                    if (std::abs(wlAmb[i][k + 2] - wlAmb[i][k + 1]) < 1 && I2 < 1) {
                        if (k + 1 - ARCS[i][j][0] > 1) {
                            cs_epoch[i][k + 1] = 1;
                            for (int l = ARCS[i][0][0]; l >= j + 1; l--) {
                                ARCS[i][l + 1][0] = ARCS[i][l][0];
                                ARCS[i][l + 1][1] = ARCS[i][l][1];

                            }
                            ARCS[i][0][0]++;
                            ARCS[i][j + 1][0] = k + 1;
                            ARCS[i][j + 1][1] = ARCS[i][j][1];
                            ARCS[i][j][1] = k;
                           

                        }
                        else {
                            for (int l = ARCS[i][j][0]; l <= k; l++) {
                                MW[i][l] = 0;
                                OBS.C1[i][l] = 0;
                                OBS.C2[i][l] = 0;
                                OBS.L1[i][l] = 0;
                                OBS.L2[i][l] = 0;
                                wlAmb[i][l] = 0;
                                GF[i][k] = 0;
                            }
                            ARCS[i][j][0] = k + 1;
                            j--;

                        }
                    }
                    else {
                        if (k + 1 - ARCS[i][j][0] > 1) {
                            MW[i][k + 1] = 0;
                            OBS.C1[i][k + 1] = 0;
                            OBS.C2[i][k + 1] = 0;
                            OBS.L1[i][k + 1] = 0;
                            OBS.L2[i][k + 1] = 0;
                            wlAmb[i][k + 1] = 0;
                            GF[i][k] = 0;
                            for (int l = ARCS[i][0][0]; l >= j + 1; l--) {
                                ARCS[i][l + 1][0] = ARCS[i][l][0];
                                ARCS[i][l + 1][1] = ARCS[i][l][1];

                            }
                            ARCS[i][0][0]++;
                            ARCS[i][j + 1][0] = k + 2;
                            ARCS[i][j + 1][1] = ARCS[i][j][1];
                            ARCS[i][j][1] = k;

                        }
                        else {
                            for (int l = ARCS[i][j][0]; l <= k + 1; l++) {
                                MW[i][l] = 0;
                                OBS.C1[i][l] = 0;
                                OBS.C2[i][l] = 0;
                                OBS.L1[i][l] = 0;
                                OBS.L2[i][l] = 0;
                                wlAmb[i][l] = 0;
                                GF[i][k] = 0;
                            }
                            ARCS[i][j][0] = k + 2;
                            j--;
                        }
                    }
                    
                    break;
                }

            }

            j++;
            
        }

        for (int j = 1; j <= 2880; j++) {
            AATR_LGF[i][j] = (c / f1) * OBS.L1[i][j] - (c / f2) * OBS.L2[i][j];

        }

        std::string rot_unit = gset.rot_unit();
        std::transform(rot_unit.begin(), rot_unit.end(), rot_unit.begin(), ::tolower);  // 转小写
        double rotNumber;

        if (rot_unit == "sec") {
            rotNumber = 30.0 * 1e16 * IONO_COEFF * (1 / (f1 * f1) - 1 / (f2 * f2));
        }
        else {

            rotNumber = 0.5 * 1e16 * IONO_COEFF * (1 / (f1 * f1) - 1 / (f2 * f2));
        }


        for (int j = 2; j <= 2880; j++) {
            if (AATR_LGF[i][j] != 0 && AATR_LGF[i][j - 1] != 0 && cs_epoch[i][j] != 1) {
                AATR_ROT[i][j] = (AATR_LGF[i][j] - AATR_LGF[i][j - 1]) / rotNumber;


            }
        }
        
    }

    double AATR_block[37] = { 0 };   
    int block_index = 0;
    for (int j = 1; j <= 2880; j += interval_size) { 
        double sum_block = 0;
        int count = 0;

        for (int k = j; k < j + interval_size; k++) {
            for (int i = 1; i <= 36; i++) {
                if (AATR_ROT[i][k] != 0) {
                    double E, A;
                    Get_Elevation(SP3[1].X[i][k] * 1000, SP3[1].Y[i][k] * 1000, SP3[1].Z[i][k] * 1000,
                        OBS.X, OBS.Y, OBS.Z, E, A);

                    double E_rad = E * PI / 180.0;
                    double cos_E = cos(E_rad);
                    double tmp = (Re * cos_E) / (Re + h);
                    double M = 1.0 / sqrt(1.0 - tmp * tmp);
                    double AATR = (1.0 / (M * M)) * AATR_ROT[i][k];

                    sum_block += AATR * AATR;
                    count++;
                }
            }
        }

        if (count > 0)
            AATR_block[block_index] = sqrt(sum_block / count);  
        else
            AATR_block[block_index] = 0;

        block_index++;
    }

    std::string aatr_out_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GAL_AATR.txt";
    std::filesystem::create_directories(std::filesystem::path(aatr_out_path).parent_path());

    std::ofstream out(aatr_out_path);
    if (!out.is_open()) {
        std::cerr << "Failed to open GAL AATR output file." << std::endl;
        return;
    }

    out << std::fixed << std::setprecision(5);
    out << std::setw(6) << "Block" << std::setw(15) << "AATR" << "\n";

    for (int k = 0; k < block_index; ++k) {
        char block_buf[20];
        sprintf(block_buf, "Block %02d:", k + 1);
        double val = AATR_block[k];
        out << std::setw(10) << block_buf << std::setw(12) << val << "\n";
    }

    out.close();

}
void AATR_EX_prepro(obs& OBS, const std::string& stationName, sp3 SP3[], const std::string& txt_output_path, t_gcfg_ppp& gset){
    double f1 = 1575.42e6;
    double f2 = 1176.45e6;
    double lambda_wl = 299792458 / (f1 - f2);
    double MW[37][2881], GF[37][2881], wlAmb[37][2881];


    double AATR_LGF[39][3000];
    double AATR_ROT[37][2888];

    double cs_epoch[37][2881];
    int arc_min_len = gset.arc_min_length();
    int interval_size = gset.aatr_interval();
    memset(AATR_ROT, 0, sizeof(AATR_ROT));
    memset(cs_epoch, 0, sizeof(cs_epoch));
    for (int i = 1; i <= 36; i++) {
        for (int j = 1; j <= 2880; j++) {
            MW[i][j] = lambda_wl * (OBS.L1[i][j] - OBS.L2[i][j]) - (f1 * OBS.C1[i][j] + f2 * OBS.C2[i][j]) / (f1 + f2);
            GF[i][j] = OBS.L1[i][j] - f1 * OBS.L2[i][j] / f2;
            wlAmb[i][j] = MW[i][j] * -1.0;
        }
    }

    for (int i = 1; i <= 36; i++) {

        int len = 2881; 
        int num_arcs = 0;

        std::vector<std::vector<int>> arcs; 
        
        Get_arc(MW[i], len, arcs, num_arcs);  
        for (int k = 0; k < num_arcs; ++k) {
            if (arcs[k][1] - arcs[k][0] < arc_min_len) {

                for (int e = arcs[k][0]; e <= arcs[k][1]; ++e) {
                    
                    OBS.C1[i][e] = 0;
                    OBS.C2[i][e] = 0;
                    OBS.L1[i][e] = 0;
                    OBS.L2[i][e] = 0; 
                }
                arcs.erase(arcs.begin() + k);
                --k;  
                --num_arcs; 
            }
        }


        int ARCS[37][3000][2];

        ARCS[i][0][0] = num_arcs;

        for (int k = 0; k < num_arcs; ++k) {
            ARCS[i][k + 1][0] = arcs[k][0];
            ARCS[i][k + 1][1] = arcs[k][1];
        }

        int j = 1;
        while (j <= ARCS[i][0][0]) {

            int e = ARCS[i][j][0];
            while (true) {
                if (e + 1 == ARCS[i][j][1] || e == ARCS[i][j][1]) {
                    break;

                }
                double amb0 = wlAmb[i][e];
                double amb1 = wlAmb[i][e + 1];
                double amb2 = wlAmb[i][e + 2];
                double gf0 = GF[i][e];
                double gf1 = GF[i][e + 1];
                double gf2 = GF[i][e + 2];
                double diff01 = std::abs(amb0 - amb1);
                double diff21 = std::abs(amb2 - amb1);
                double gfDiff01 = std::abs(gf0 - gf1);
                double gfDiff21 = std::abs(gf2 - gf1);
                if (diff01 > 1 || diff21 > 1 || gfDiff01 > 1 || gfDiff21 > 1) {
                    MW[i][e] = 0;
                    OBS.C1[i][e] = 0;
                    OBS.C2[i][e] = 0;
                    OBS.L1[i][e] = 0;
                    OBS.L2[i][e] = 0;
                    wlAmb[i][e] = 0;
                    GF[i][e] = 0;
                    e++;
                    ARCS[i][j][0] = e;
                }
                else {
                    ARCS[i][j][0] = e;
                    break;

                }
            }
            if (ARCS[i][j][1] - ARCS[i][j][0] < arc_min_len) {
                for (int k = ARCS[i][j][0]; k <= ARCS[i][j][1]; k++) {
                    OBS.C1[i][k] = 0;
                    OBS.C2[i][k] = 0;
                    OBS.L1[i][k] = 0;
                    OBS.L2[i][k] = 0;
                    MW[i][k] = 0;
                    wlAmb[i][e] = 0;
                    GF[i][k] = 0;

                }
                for (int k = j; k < ARCS[i][0][0]; k++) {
                    ARCS[i][k][0] = ARCS[i][k + 1][0];
                    ARCS[i][k][1] = ARCS[i][k + 1][1];
                }
                ARCS[i][0][0]--;
            }

            double ave_N[3000];
            double sigma2[3000];
            double sigma[3000];
            ave_N[1] = wlAmb[i][ARCS[i][j][0]];
            sigma2[1] = 0;
            sigma[1] = 0;
            int count = 2;
            for (int k = ARCS[i][j][0] + 1; k < ARCS[i][j][1]; k++) {

                ave_N[count] = ave_N[count - 1] + (wlAmb[i][k] - ave_N[count - 1]) / count;
                sigma2[count] = sigma2[count - 1] + ((wlAmb[i][k] - ave_N[count - 1]) * (wlAmb[i][k] - ave_N[count - 1]) - sigma2[count - 1]) / count;
                sigma[count] = sqrt(sigma2[count]);
                double T = std::abs(wlAmb[i][k + 1] - ave_N[count]);
                double I1 = std::abs(GF[i][k + 1] - GF[i][k]);
                if (T < 4 * sigma[count] && I1 < 0.28) {
                    count++;
                    continue;

                }
                else {
                    if (k + 1 == ARCS[i][j][1]) {

                        if (k + 1 - ARCS[i][j][0] > 1) {
                            MW[i][k + 1] = 0;
                            OBS.C1[i][k + 1] = 0;
                            OBS.C2[i][k + 1] = 0;
                            OBS.L1[i][k + 1] = 0;
                            OBS.L2[i][k + 1] = 0;
                            wlAmb[i][k + 1] = 0;
                            GF[i][k] = 0;
                            ARCS[i][j][1] = k;

                        }
                        else {
                            for (int l = ARCS[i][j][0]; l <= k + 1; l++) {
                                MW[i][l] = 0;
                                OBS.C1[i][l] = 0;
                                OBS.C2[i][l] = 0;
                                OBS.L1[i][l] = 0;
                                OBS.L2[i][l] = 0;
                                wlAmb[i][l] = 0;
                                GF[i][k] = 0;

                            }
                            for (int x = j; x < ARCS[i][0][0]; x++) {
                                ARCS[i][x][0] = ARCS[i][x + 1][0];
                                ARCS[i][x][1] = ARCS[i][x + 1][1];
                            }
                            j--;
                            ARCS[i][0][0]--;
                        }
                        break;
                    }
                    double I2 = std::abs(GF[i][k + 2] - GF[i][k + 1]);
                    if (std::abs(wlAmb[i][k + 2] - wlAmb[i][k + 1]) < 1 && I2 < 1) {
                        if (k + 1 - ARCS[i][j][0] > 1) {
                            cs_epoch[i][k + 1] = 1;
                            for (int l = ARCS[i][0][0]; l >= j + 1; l--) {
                                ARCS[i][l + 1][0] = ARCS[i][l][0];
                                ARCS[i][l + 1][1] = ARCS[i][l][1];

                            }
                            ARCS[i][0][0]++;
                            ARCS[i][j + 1][0] = k + 1;
                            ARCS[i][j + 1][1] = ARCS[i][j][1];
                            ARCS[i][j][1] = k;
                            

                        }
                        else {
                            for (int l = ARCS[i][j][0]; l <= k; l++) {
                                 MW[i][l] = 0;
                                OBS.C1[i][l] = 0;
                                OBS.C2[i][l] = 0;
                                OBS.L1[i][l] = 0;
                                OBS.L2[i][l] = 0;
                                wlAmb[i][l] = 0;
                                GF[i][k] = 0;
                            }
                            ARCS[i][j][0] = k + 1;
                            j--;

                        }
                    }
                    else {
                        if (k + 1 - ARCS[i][j][0] > 1) {
                            MW[i][k + 1] = 0;
                            OBS.C1[i][k + 1] = 0;
                            OBS.C2[i][k + 1] = 0;
                            OBS.L1[i][k + 1] = 0;
                            OBS.L2[i][k + 1] = 0;
                            wlAmb[i][k + 1] = 0;
                            GF[i][k] = 0;
                            for (int l = ARCS[i][0][0]; l >= j + 1; l--) {
                                ARCS[i][l + 1][0] = ARCS[i][l][0];
                                ARCS[i][l + 1][1] = ARCS[i][l][1];

                            }
                            ARCS[i][0][0]++;
                            ARCS[i][j + 1][0] = k + 2;
                            ARCS[i][j + 1][1] = ARCS[i][j][1];
                            ARCS[i][j][1] = k;

                        }
                        else {
                            for (int l = ARCS[i][j][0]; l <= k + 1; l++) {
                                MW[i][l] = 0;
                                OBS.C1[i][l] = 0;
                                OBS.C2[i][l] = 0;
                                OBS.L1[i][l] = 0;
                                OBS.L2[i][l] = 0;
                                wlAmb[i][l] = 0;
                                GF[i][k] = 0;
                            }
                            ARCS[i][j][0] = k + 2;
                            j--;
                        }
                    }
                    
                    break;
                }

            }

            j++;
            
        }

        for (int j = 1; j <= 2880; j++) {
            AATR_LGF[i][j] = (c / f1) * OBS.L1[i][j] - (c / f2) * OBS.L2[i][j];

        }

        std::string rot_unit = gset.rot_unit();
        std::transform(rot_unit.begin(), rot_unit.end(), rot_unit.begin(), ::tolower);  // 转小写
        double rotNumber;

        if (rot_unit == "sec") {
            rotNumber = 30.0 * 1e16 * IONO_COEFF * (1 / (f1 * f1) - 1 / (f2 * f2));
        }
        else {
         
            rotNumber = 0.5 * 1e16 * IONO_COEFF * (1 / (f1 * f1) - 1 / (f2 * f2));
        }


        for (int j = 2; j <= 2880; j++) {
            if (AATR_LGF[i][j] != 0 && AATR_LGF[i][j - 1] != 0 && cs_epoch[i][j] != 1) {
                AATR_ROT[i][j] = (AATR_LGF[i][j] - AATR_LGF[i][j - 1]) / rotNumber;


            }
        }

    }

    double AATR_block[37] = { 0 }; 
    int block_index = 0;
    for (int j = 1; j <= 2880; j += interval_size) { 
        double sum_block = 0;
        int count = 0;

        for (int k = j; k < j + interval_size; k++) {
            for (int i = 1; i <= 36; i++) {
                if (AATR_ROT[i][k] != 0) { 
                    double E, A;
                    Get_Elevation(SP3[1].X[i][k] * 1000, SP3[1].Y[i][k] * 1000, SP3[1].Z[i][k] * 1000,
                        OBS.X, OBS.Y, OBS.Z, E, A);

                    double E_rad = E * PI / 180.0;
                    double cos_E = cos(E_rad);
                    double tmp = (Re * cos_E) / (Re + h);
                    double M = 1.0 / sqrt(1.0 - tmp * tmp);
                    double AATR = (1.0 / (M * M)) * AATR_ROT[i][k];

                    sum_block += AATR * AATR;
                    count++;
                }
            }
        }

        if (count > 0)
            AATR_block[block_index] = sqrt(sum_block / count);  
        else
            AATR_block[block_index] = 0;

        block_index++;
    }


    std::string aatr_out_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GALX_AATR.txt";
    std::filesystem::create_directories(std::filesystem::path(aatr_out_path).parent_path());

    std::ofstream out(aatr_out_path);
    if (!out.is_open()) {
        std::cerr << "Failed to open GALX AATR output file." << std::endl;
        return;
    }

    out << std::fixed << std::setprecision(5);
    out << std::setw(6) << "Block" << std::setw(15) << "AATR" << "\n";

    for (int k = 0; k < block_index; ++k) {
        char block_buf[20];
        sprintf(block_buf, "Block %02d:", k + 1);
        double val = AATR_block[k];
        out << std::setw(10) << block_buf << std::setw(12) << val << "\n";
    }

    out.close();

}
void AATR_C_prepro(obs& OBS, const std::string& stationName, bool isC7IAllZero, sp3 SP3[], const std::string& txt_output_path, t_gcfg_ppp& gset) {
    double f1 = 1561.098e6;
    double f2;
    int arc_min_len = gset.arc_min_length();
    int interval_size = gset.aatr_interval();
    // Select f2 based on whether C7I is all zeros
    if (isC7IAllZero) {
        f2 = 1268.52e6; 
    }
    else {
        f2 = 1207.140e6;
    }
    double lambda_wl = 299792458 / (f1 - f2);
    double MW[47][2881], GF[47][2881], wlAmb[47][2881];


    double AATR_LGF[49][3000];
    double AATR_ROT[47][2888];

    double cs_epoch[47][2881];
    memset(AATR_ROT, 0, sizeof(AATR_ROT));
    memset(cs_epoch, 0, sizeof(cs_epoch));
    for (int i = 1; i <= 46; i++) {
        for (int j = 1; j <= 2880; j++) {
            MW[i][j] = lambda_wl * (OBS.L1[i][j] - OBS.L2[i][j]) - (f1 * OBS.C1[i][j] + f2 * OBS.C2[i][j]) / (f1 + f2);
            GF[i][j] = OBS.L1[i][j] - f1 * OBS.L2[i][j] / f2;
            wlAmb[i][j] = MW[i][j] * -1.0;
        }
    }
    
    for (int i = 1; i <= 46; i++) {

        int len = 2881;       
        int num_arcs = 0;      

        std::vector<std::vector<int>> arcs;
       
        Get_arc(MW[i], len, arcs, num_arcs);  
        for (int k = 0; k < num_arcs; ++k) {
            if (arcs[k][1] - arcs[k][0] < arc_min_len) {

                for (int e = arcs[k][0]; e <= arcs[k][1]; ++e) {
                  
                    OBS.C1[i][e] = 0;
                    OBS.C2[i][e] = 0;
                    OBS.L1[i][e] = 0;
                    OBS.L2[i][e] = 0; 
                }
                arcs.erase(arcs.begin() + k);
                --k; 
                --num_arcs; 
            }
        }


        int ARCS[47][3000][2];

        ARCS[i][0][0] = num_arcs;

        for (int k = 0; k < num_arcs; ++k) {
            ARCS[i][k + 1][0] = arcs[k][0];
            ARCS[i][k + 1][1] = arcs[k][1];
        }
        int j = 1;
        while (j <= ARCS[i][0][0]) {

            int e = ARCS[i][j][0];
            while (true) {
                if (e + 1 == ARCS[i][j][1] || e == ARCS[i][j][1]) {
                    break;

                }
                double amb0 = wlAmb[i][e];
                double amb1 = wlAmb[i][e + 1];
                double amb2 = wlAmb[i][e + 2];
                double gf0 = GF[i][e];
                double gf1 = GF[i][e + 1];
                double gf2 = GF[i][e + 2];
                double diff01 = std::abs(amb0 - amb1);
                double diff21 = std::abs(amb2 - amb1);
                double gfDiff01 = std::abs(gf0 - gf1);
                double gfDiff21 = std::abs(gf2 - gf1);
                if (diff01 > 1 || diff21 > 1 || gfDiff01 > 1 || gfDiff21 > 1) {
                    MW[i][e] = 0;
                    OBS.C1[i][e] = 0;
                    OBS.C2[i][e] = 0;
                    OBS.L1[i][e] = 0;
                    OBS.L2[i][e] = 0;
                    wlAmb[i][e] = 0;
                    GF[i][e] = 0;
                    e++;
                    ARCS[i][j][0] = e;
                }
                else {
                    ARCS[i][j][0] = e;
                    break;

                }
            }
            if (ARCS[i][j][1] - ARCS[i][j][0] < arc_min_len) {
               
                for (int k = ARCS[i][j][0]; k <= ARCS[i][j][1]; k++) {
                    OBS.C1[i][k] = 0;
                    OBS.C2[i][k] = 0;
                    OBS.L1[i][k] = 0;
                    OBS.L2[i][k] = 0;
                    MW[i][k] = 0;
                    wlAmb[i][e] = 0;
                    GF[i][k] = 0;

                }
                for (int k = j; k < ARCS[i][0][0]; k++) {
                    ARCS[i][k][0] = ARCS[i][k + 1][0];
                    ARCS[i][k][1] = ARCS[i][k + 1][1];
                }
                ARCS[i][0][0]--;
            }


            double ave_N[3000];
            double sigma2[3000];
            double sigma[3000];
            ave_N[1] = wlAmb[i][ARCS[i][j][0]];
            sigma2[1] = 0;
            sigma[1] = 0;
            int count = 2;
            for (int k = ARCS[i][j][0] + 1; k < ARCS[i][j][1]; k++) {

                ave_N[count] = ave_N[count - 1] + (wlAmb[i][k] - ave_N[count - 1]) / count;
                sigma2[count] = sigma2[count - 1] + ((wlAmb[i][k] - ave_N[count - 1]) * (wlAmb[i][k] - ave_N[count - 1]) - sigma2[count - 1]) / count;
                sigma[count] = sqrt(sigma2[count]);
                double T = std::abs(wlAmb[i][k + 1] - ave_N[count]);
                double I1 = std::abs(GF[i][k + 1] - GF[i][k]);
                if (T < 4 * sigma[count] && I1 < 0.28) {
                    count++;
                    continue;

                }
                else {
                    if (k + 1 == ARCS[i][j][1]) {

                        if (k + 1 - ARCS[i][j][0] > 1) {
                            MW[i][k + 1] = 0;
                            OBS.C1[i][k + 1] = 0;
                            OBS.C2[i][k + 1] = 0;
                            OBS.L1[i][k + 1] = 0;
                            OBS.L2[i][k + 1] = 0;
                            wlAmb[i][k + 1] = 0;
                            GF[i][k] = 0;
                            ARCS[i][j][1] = k;

                        }
                        else {
                            for (int l = ARCS[i][j][0]; l <= k + 1; l++) {
                                MW[i][l] = 0;
                                OBS.C1[i][l] = 0;
                                OBS.C2[i][l] = 0;
                                OBS.L1[i][l] = 0;
                                OBS.L2[i][l] = 0;
                                wlAmb[i][l] = 0;
                                GF[i][k] = 0;

                            }
                            for (int x = j; x < ARCS[i][0][0]; x++) {
                                ARCS[i][x][0] = ARCS[i][x + 1][0];
                                ARCS[i][x][1] = ARCS[i][x + 1][1];
                            }
                            j--;
                            ARCS[i][0][0]--;
                        }
                        break;
                    }
                    double I2 = std::abs(GF[i][k + 2] - GF[i][k + 1]);
                    if (std::abs(wlAmb[i][k + 2] - wlAmb[i][k + 1]) < 1 && I2 < 1) {
                        if (k + 1 - ARCS[i][j][0] > 1) {
                            cs_epoch[i][k + 1] = 1;
                            for (int l = ARCS[i][0][0]; l >= j + 1; l--) {
                                ARCS[i][l + 1][0] = ARCS[i][l][0];
                                ARCS[i][l + 1][1] = ARCS[i][l][1];

                            }
                            ARCS[i][0][0]++;
                            ARCS[i][j + 1][0] = k + 1;
                            ARCS[i][j + 1][1] = ARCS[i][j][1];
                            ARCS[i][j][1] = k;
                            
                        }
                        else {
                            for (int l = ARCS[i][j][0]; l <= k; l++) {
                                MW[i][l] = 0;
                                OBS.C1[i][l] = 0;
                                OBS.C2[i][l] = 0;
                                OBS.L1[i][l] = 0;
                                OBS.L2[i][l] = 0;
                                wlAmb[i][l] = 0;
                                GF[i][k] = 0;
                            }
                            ARCS[i][j][0] = k + 1;
                            j--;

                        }
                    }
                    else {
                        if (k + 1 - ARCS[i][j][0] > 1) {
                            MW[i][k + 1] = 0;
                            OBS.C1[i][k + 1] = 0;
                            OBS.C2[i][k + 1] = 0;
                            OBS.L1[i][k + 1] = 0;
                            OBS.L2[i][k + 1] = 0;
                            wlAmb[i][k + 1] = 0;
                            GF[i][k] = 0;
                            for (int l = ARCS[i][0][0]; l >= j + 1; l--) {
                                ARCS[i][l + 1][0] = ARCS[i][l][0];
                                ARCS[i][l + 1][1] = ARCS[i][l][1];

                            }
                            ARCS[i][0][0]++;
                            ARCS[i][j + 1][0] = k + 2;
                            ARCS[i][j + 1][1] = ARCS[i][j][1];
                            ARCS[i][j][1] = k;

                        }
                        else {
                            for (int l = ARCS[i][j][0]; l <= k + 1; l++) {
                                MW[i][l] = 0;
                                OBS.C1[i][l] = 0;
                                OBS.C2[i][l] = 0;
                                OBS.L1[i][l] = 0;
                                OBS.L2[i][l] = 0;
                                wlAmb[i][l] = 0;
                                GF[i][k] = 0;
                            }
                            ARCS[i][j][0] = k + 2;
                            j--;
                        }
                    }
                    
                    break;
                }

            }

            j++;

        }

        for (int j = 1; j <= 2880; j++) {
            AATR_LGF[i][j] = (c / f1) * OBS.L1[i][j] - (c / f2) * OBS.L2[i][j];

        }

        std::string rot_unit = gset.rot_unit();
        std::transform(rot_unit.begin(), rot_unit.end(), rot_unit.begin(), ::tolower);  // 转小写
        double rotNumber;

        if (rot_unit == "sec") {
            rotNumber = 30.0 * 1e16 * IONO_COEFF * (1 / (f1 * f1) - 1 / (f2 * f2));
        }
        else {
   
            rotNumber = 0.5 * 1e16 * IONO_COEFF * (1 / (f1 * f1) - 1 / (f2 * f2));
        }


        for (int j = 2; j <= 2880; j++) {
            if (AATR_LGF[i][j] != 0 && AATR_LGF[i][j - 1] != 0 && cs_epoch[i][j] != 1) {
                AATR_ROT[i][j] = (AATR_LGF[i][j] - AATR_LGF[i][j - 1]) / rotNumber;


            }
        }
        
    }
    
    double AATR_block[47] = { 0 }; 

    int block_index = 0;
    for (int j = 1; j <= 2880; j += interval_size) { 
        double sum_block = 0;
        int count = 0;

        for (int k = j; k < j + interval_size; k++) {
            for (int i = 1; i <= 46; i++) {
                if (AATR_ROT[i][k] != 0) { 
                    double E, A;
                    Get_Elevation(SP3[1].X[i][k] * 1000, SP3[1].Y[i][k] * 1000, SP3[1].Z[i][k] * 1000,
                        OBS.X, OBS.Y, OBS.Z, E, A);

                    double E_rad = E * PI / 180.0;
                    double cos_E = cos(E_rad);
                    double tmp = (Re * cos_E) / (Re + h);
                    double M = 1.0 / sqrt(1.0 - tmp * tmp);
                    double AATR = (1.0 / (M * M)) * AATR_ROT[i][k];

                    sum_block += AATR * AATR;
                    count++;
                }
            }
        }

        if (count > 0)
            AATR_block[block_index] = sqrt(sum_block / count); 
        else
            AATR_block[block_index] = 0;

        block_index++;
    }

    std::string aatr_out_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_BDS_AATR.txt";
    std::filesystem::create_directories(std::filesystem::path(aatr_out_path).parent_path());

    std::ofstream out(aatr_out_path);
    if (!out.is_open()) {
        std::cerr << "Failed to open BDS AATR output file." << std::endl;
        return;
    }

    out << std::fixed << std::setprecision(5);
    out << std::setw(6) << "Block" << std::setw(15) << "AATR" << "\n";

    for (int k = 0; k < block_index; ++k) {
        char block_buf[20];
        sprintf(block_buf, "Block %02d:", k + 1);
        double val = AATR_block[k];
        out << std::setw(10) << block_buf << std::setw(12) << val << "\n";
    }

    out.close();
}

