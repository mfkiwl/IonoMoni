#include "calc_ROTI.h"
#include "calc_elevation.h"
#include "constants.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <filesystem>
#include <iostream>   
#include "gset/gsetproc.h"
#include "gcfg_ppp.h"

#include <vector>
#include "obs.h"  

void calc_roti_GPS(const obs& OBS, const std::string& stationName, const sp3& SP3, bool isC1WAllZero, const std::string& txt_output_path, t_gcfg_ppp& gset) {
    int window_size = gset.roti_window();
    double f1 = 1575.42e6, f2 = 1227.6e6;
    double LGF[35][3000], ROT[33][2888], ROTI[33][2888];
    memset(ROT, 0, sizeof(ROT));
    memset(ROTI, 0, sizeof(ROTI));

    for (int i = 1; i <= 32; i++) {
        for (int j = 1; j <= 2880; j++) {
            LGF[i][j] = (c / f1) * OBS.L1[i][j] - (c / f2) * OBS.L2[i][j];
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

        // Calculate ROT 
        for (int j = 2; j <= 2880; j++)
            if (LGF[i][j] && LGF[i][j - 1])
                ROT[i][j] = (LGF[i][j] - LGF[i][j - 1]) / rotNumber;

        // Calculate ROTI (standard deviation of ROT in sliding window)
        for (int j = window_size; j <= 2880; j++) {
            double sum = 0.0, var = 0.0;
            bool valid = true;

            for (int k = j - window_size + 1; k <= j; k++) {
                if (ROT[i][k] == 0.0) {
                    valid = false;
                    break;
                }
                sum += ROT[i][k];
            }
            if (!valid) continue;

            double mean = sum / window_size;
            for (int k = j - window_size + 1; k <= j; k++)
                var += (ROT[i][k] - mean) * (ROT[i][k] - mean);

            ROTI[i][j] = std::sqrt(var / window_size);
        }
    }

    // Output ROTI and ROT results to file
    std::string gps_roti_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GPS_ROTI.txt";
    std::string gps_rot_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GPS_ROT.txt";

    std::ofstream out(gps_roti_path), rot(gps_rot_path);
    if (!out.is_open() || !rot.is_open()) {
        std::cerr << "Failed to open GPS ROTI/ROT output file." << std::endl;
        return;
    }

    out << std::fixed << std::setprecision(4);
    rot << std::fixed << std::setprecision(4);

    // Output headers
    out << std::setw(12) << "Epoch \\ PRN";
    rot << std::setw(12) << "Epoch \\ PRN";
    for (int i = 1; i <= 32; ++i) {
        char prn_buf[10];
        sprintf(prn_buf, "G%02d", i);
        out << std::setw(11) << prn_buf;
        rot << std::setw(11) << prn_buf;
    }
    out << "\n";
    rot << "\n";

    // Output ROTI/ROT values for each epoch and satellite
    for (int j = 1; j <= 2880; ++j) {
        char epoch_buf[20];
        sprintf(epoch_buf, "Epoch %04d:", j);
        out << std::setw(12) << epoch_buf;
        rot << std::setw(12) << epoch_buf;

        for (int i = 1; i <= 32; ++i) {
            double rti = ROTI[i][j];
            double rt = ROT[i][j];
            out << std::setw(11) << rti;
            rot << std::setw(11) << rt;
        }
        out << "\n";
        rot << "\n";
    }
    out.close();
    rot.close();
}



//// Calculate S4C scintillation index for S1 and S2 channels
//// OBS: observation data, numSats: number of satellites, numEpochs: number of epochs
//// n: window size, outS1/outS2: output file paths for S1 and S2
//void calc_S4C(const obs& OBS, int numSats, int numEpochs, int n, const std::string& outS1, const std::string& outS2) {
//    std::vector<std::vector<double>> S4C_S1(numSats + 1, std::vector<double>(numEpochs + 1, 0.0));
//    std::vector<std::vector<double>> S4C_S2(numSats + 1, std::vector<double>(numEpochs + 1, 0.0));
//
//    // Loop through each satellite
//    for (int prn = 1; prn <= numSats; ++prn) {
//        for (int k = 2 * n; k <= numEpochs; ++k) {
//            std::vector<double> s1_detrended, s2_detrended;
//            bool valid = true;
//            // Sliding window for local mean normalization
//            for (int m = k - n + 1; m <= k; ++m) {
//                if (m - n < 1) { valid = false; break; }
//                double mean_S1 = 0.0, mean_S2 = 0.0;
//                bool v1 = true, v2 = true;
//                for (int j = m - n; j < m; ++j) {
//                    // Transform S/N0 to linear scale (dBHz to linear)
//                    if (OBS.S1[prn][j] <= 0.0) v1 = false;
//                    if (OBS.S2[prn][j] <= 0.0) v2 = false;
//                    double S1_linear = pow(10.0, 0.1 * OBS.S1[prn][j]);
//                    double S2_linear = pow(10.0, 0.1 * OBS.S2[prn][j]);
//                    mean_S1 += S1_linear;
//                    mean_S2 += S2_linear;
//                }
//                // Current epoch S/N0 transformation
//                double S1k = pow(10.0, 0.1 * OBS.S1[prn][m]);
//                double S2k = pow(10.0, 0.1 * OBS.S2[prn][m]);
//                if (!v1 || OBS.S1[prn][m] <= 0.0) { valid = false; break; }
//                if (!v2 || OBS.S2[prn][m] <= 0.0) { valid = false; break; }
//                mean_S1 /= n;
//                mean_S2 /= n;
//                s1_detrended.push_back(S1k / mean_S1);
//                s2_detrended.push_back(S2k / mean_S2);
//            }
//            if (!valid || s1_detrended.size() != n || s2_detrended.size() != n) continue;
//            double mean1 = 0.0, mean_sq1 = 0.0;
//            double mean2 = 0.0, mean_sq2 = 0.0;
//            // Compute variance for S4C
//            for (int i = 0; i < n; ++i) {
//                mean1 += s1_detrended[i];
//                mean_sq1 += s1_detrended[i] * s1_detrended[i];
//                mean2 += s2_detrended[i];
//                mean_sq2 += s2_detrended[i] * s2_detrended[i];
//            }
//            mean1 /= n; mean_sq1 /= n;
//            mean2 /= n; mean_sq2 /= n;
//            S4C_S1[prn][k] = (mean_sq1 - mean1 * mean1 > 0.0) ? std::sqrt(mean_sq1 - mean1 * mean1) : 0.0;
//            S4C_S2[prn][k] = (mean_sq2 - mean2 * mean2 > 0.0) ? std::sqrt(mean_sq2 - mean2 * mean2) : 0.0;
//        }
//    }
//
//    // Output S4C results to TXT files
//    std::ofstream fout1(outS1), fout2(outS2);
//    fout1 << std::fixed << std::setprecision(5);
//    fout2 << std::fixed << std::setprecision(5);
//    for (int j = 1; j <= numEpochs; ++j) {
//        for (int prn = 1; prn <= numSats; ++prn) {
//            fout1 << S4C_S1[prn][j] << (prn < numSats ? "\t" : "");
//            fout2 << S4C_S2[prn][j] << (prn < numSats ? "\t" : "");
//        }
//        fout1 << "\n";
//        fout2 << "\n";
//    }
//    fout1.close(); fout2.close();
//}


void calc_roti_BDS(const obs& OBS, const std::string& stationName, const sp3& SP3, bool isC7IAllZero, const std::string& txt_output_path, t_gcfg_ppp& gset){
    int window_size = gset.roti_window();

    double f1 = 1561.098e6;
    double f2 = isC7IAllZero ? 1268.52e6 : 1207.14e6;
    double LGF[49][3000], ROT[47][2888], ROTI[47][2888];
    memset(ROT, 0, sizeof(ROT));
    memset(ROTI, 0, sizeof(ROTI));

    for (int i = 1; i <= 46; i++) {
        for (int j = 1; j <= 2880; j++)
            LGF[i][j] = (c / f1) * OBS.L1[i][j] - (c / f2) * OBS.L2[i][j];

        std::string rot_unit = gset.rot_unit();
        std::transform(rot_unit.begin(), rot_unit.end(), rot_unit.begin(), ::tolower); 
        double rotNumber;

        if (rot_unit == "sec") {
            rotNumber = 30.0 * 1e16 * IONO_COEFF * (1 / (f1 * f1) - 1 / (f2 * f2));
        }
        else {
            rotNumber = 0.5 * 1e16 * IONO_COEFF * (1 / (f1 * f1) - 1 / (f2 * f2));
        }
        for (int j = 2; j <= 2880; j++)
            if (LGF[i][j] && LGF[i][j - 1])
                ROT[i][j] = (LGF[i][j] - LGF[i][j - 1]) / rotNumber;

        for (int j = window_size; j <= 2880; j++) {
            double sum = 0.0, var = 0.0;
            bool valid = true;

            for (int k = j - window_size + 1; k <= j; k++) {
                if (ROT[i][k] == 0.0) {
                    valid = false;
                    break;
                }
                sum += ROT[i][k];
            }

            if (!valid) continue;

            double mean = sum / window_size;
            for (int k = j - window_size + 1; k <= j; k++)
                var += (ROT[i][k] - mean) * (ROT[i][k] - mean);

            ROTI[i][j] = std::sqrt(var / window_size);
        }
    }

    std::string bds_roti_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_BDS_ROTI.txt";
    std::string bds_rot_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_BDS_ROT.txt";

    std::ofstream out(bds_roti_path), rot(bds_rot_path);
    if (!out.is_open() || !rot.is_open()) {
        std::cerr << "Failed to open BDS ROTI/ROT output file." << std::endl;
        return;
    }

    out << std::fixed << std::setprecision(4);
    rot << std::fixed << std::setprecision(4);

    out << std::setw(12) << "Epoch \\ PRN";
    rot << std::setw(12) << "Epoch \\ PRN";
    for (int i = 1; i <= 46; ++i) {
        char prn_buf[10];
        sprintf(prn_buf, "C%02d", i);  
        out << std::setw(11) << prn_buf;
        rot << std::setw(11) << prn_buf;
    }
    out << "\n";
    rot << "\n";

    for (int j = 1; j <= 2880; ++j) {
        char epoch_buf[20];
        sprintf(epoch_buf, "Epoch %04d:", j);
        out << std::setw(12) << epoch_buf;
        rot << std::setw(12) << epoch_buf;

        for (int i = 1; i <= 46; ++i) {
            double rti = ROTI[i][j];
            double rt = ROT[i][j];
            out << std::setw(11) << rti;
            rot << std::setw(11) << rt;
        }
        out << "\n";
        rot << "\n";
    }

    out.close();
    rot.close();

}

void calc_roti_GLO(const obs& OBS, const std::string& stationName, const sp3& SP3, const std::string& txt_output_path, t_gcfg_ppp& gset){
    int window_size = gset.roti_window();

    int fre[25] = { 99, 1, -4, 5, 6, 1, -4, 5, 6, -2, -7, 0, -1, -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2 };
    double f1[25], f2[25];
    double LGF[27][3000], ROT[25][2888], ROTI[25][2888];
    memset(ROT, 0, sizeof(ROT));
    memset(ROTI, 0, sizeof(ROTI));

    for (int i = 1; i <= 24; i++) {
        f1[i] = (1602.0 + fre[i] * 0.5625) * 1e6;
        f2[i] = (1246.0 + fre[i] * 0.4375) * 1e6;

        for (int j = 1; j <= 2880; j++)
            LGF[i][j] = (c / f1[i]) * OBS.L1[i][j] - (c / f2[i]) * OBS.L2[i][j];

        std::string rot_unit = gset.rot_unit();
        std::transform(rot_unit.begin(), rot_unit.end(), rot_unit.begin(), ::tolower);  
        double rotNumber;

        if (rot_unit == "sec") {
            rotNumber = 30.0 * 1e16 * IONO_COEFF * (1 / (f1[i] * f1[i]) - 1 / (f2[i] * f2[i]));
        }
        else {
            rotNumber = 0.5 * 1e16 * IONO_COEFF * (1 / (f1[i] * f1[i]) - 1 / (f2[i] * f2[i]));
        }
        for (int j = 2; j <= 2880; j++)
            if (LGF[i][j] && LGF[i][j - 1])
                ROT[i][j] = (LGF[i][j] - LGF[i][j - 1]) / rotNumber;

        for (int j = window_size; j <= 2880; j++) {
            double sum = 0.0, var = 0.0;
            bool valid = true;

            for (int k = j - window_size + 1; k <= j; k++) {
                if (ROT[i][k] == 0.0) {
                    valid = false;
                    break;
                }
                sum += ROT[i][k];
            }

            if (!valid) continue;

            double mean = sum / window_size;
            for (int k = j - window_size + 1; k <= j; k++)
                var += (ROT[i][k] - mean) * (ROT[i][k] - mean);

            ROTI[i][j] = std::sqrt(var / window_size);
        }
    }

    std::string glo_roti_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GLO_ROTI.txt";
    std::string glo_rot_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GLO_ROT.txt";

    std::ofstream out(glo_roti_path), rot(glo_rot_path);
    if (!out.is_open() || !rot.is_open()) {
        std::cerr << "Failed to open GLO ROTI/ROT output file." << std::endl;
        return;
    }

    out << std::fixed << std::setprecision(4);
    rot << std::fixed << std::setprecision(4);

    out << std::setw(12) << "Epoch \\ PRN";
    rot << std::setw(12) << "Epoch \\ PRN";
    for (int i = 1; i <= 24; ++i) {
        char prn_buf[10];
        sprintf(prn_buf, "R%02d", i);  
        out << std::setw(11) << prn_buf;
        rot << std::setw(11) << prn_buf;
    }
    out << "\n";
    rot << "\n";

    for (int j = 1; j <= 2880; ++j) {
        char epoch_buf[20];
        sprintf(epoch_buf, "Epoch %04d:", j);
        out << std::setw(12) << epoch_buf;
        rot << std::setw(12) << epoch_buf;

        for (int i = 1; i <= 24; ++i) {
            double rti = ROTI[i][j];
            double rt = ROT[i][j];
            out << std::setw(11) << rti;
            rot << std::setw(11) << rt;
        }

        out << "\n";
        rot << "\n";
    }

    out.close();
    rot.close();

}

void calc_roti_GAL(const obs& OBS, const std::string& stationName, const sp3& SP3, const std::string& txt_output_path, t_gcfg_ppp& gset){
    int window_size = gset.roti_window();
    double f1 = 1575.42e6, f2 = 1176.45e6;
    double LGF[39][3000], ROT[37][2888], ROTI[37][2888];
    memset(ROT, 0, sizeof(ROT));
    memset(ROTI, 0, sizeof(ROTI));

    for (int i = 1; i <= 36; i++) {
        for (int j = 1; j <= 2880; j++)
            LGF[i][j] = (c / f1) * OBS.L1[i][j] - (c / f2) * OBS.L2[i][j];

        std::string rot_unit = gset.rot_unit();
        std::transform(rot_unit.begin(), rot_unit.end(), rot_unit.begin(), ::tolower);  
        double rotNumber;

        if (rot_unit == "sec") {
            rotNumber = 30.0 * 1e16 * IONO_COEFF * (1 / (f1 * f1) - 1 / (f2 * f2));
        }
        else {
            rotNumber = 0.5 * 1e16 * IONO_COEFF * (1 / (f1 * f1) - 1 / (f2 * f2));
        }
        for (int j = 2; j <= 2880; j++)
            if (LGF[i][j] && LGF[i][j - 1])
                ROT[i][j] = (LGF[i][j] - LGF[i][j - 1]) / rotNumber;

        for (int j = window_size; j <= 2880; j++) {
            double sum = 0.0, var = 0.0;
            bool valid = true;

            for (int k = j - window_size + 1; k <= j; k++) {
                if (ROT[i][k] == 0.0) {
                    valid = false;
                    break;
                }
                sum += ROT[i][k];
            }

            if (!valid) continue;

            double mean = sum / window_size;
            for (int k = j - window_size + 1; k <= j; k++)
                var += (ROT[i][k] - mean) * (ROT[i][k] - mean);

            ROTI[i][j] = std::sqrt(var / window_size);
        }
    }

    std::string gal_roti_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GAL_ROTI.txt";
    std::string gal_rot_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GAL_ROT.txt";

    std::ofstream out(gal_roti_path), rot(gal_rot_path);
    if (!out.is_open() || !rot.is_open()) {
        std::cerr << "Failed to open GAL ROTI/ROT output file." << std::endl;
        return;
    }

    out << std::fixed << std::setprecision(4);
    rot << std::fixed << std::setprecision(4);

    out << std::setw(12) << "Epoch \\ PRN";
    rot << std::setw(12) << "Epoch \\ PRN";
    for (int i = 1; i <= 36; ++i) {
        char prn_buf[10];
        sprintf(prn_buf, "E%02d", i);  
        out << std::setw(11) << prn_buf;
        rot << std::setw(11) << prn_buf;
    }
    out << "\n";
    rot << "\n";

    for (int j = 1; j <= 2880; ++j) {
        char epoch_buf[20];
        sprintf(epoch_buf, "Epoch %04d:", j);
        out << std::setw(12) << epoch_buf;
        rot << std::setw(12) << epoch_buf;

        for (int i = 1; i <= 36; ++i) {
            double rti = ROTI[i][j];
            double rt = ROT[i][j];
            out << std::setw(11) << rti;
            rot << std::setw(11) << rt;
        }

        out << "\n";
        rot << "\n";
    }

    out.close();
    rot.close();

}

