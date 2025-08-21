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
#include <string>
#include <limits>

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

void calc_S4C(const obs& OBS,
    int numSats, int numEpochs,
    int n_trend, int L_stat,
    const std::string& txt_output_path,
    const std::string& systemTag,
    int hourWanted)
{
    const double NaN = std::numeric_limits<double>::quiet_NaN();

    std::vector<std::vector<double>> S4C_S1(numSats + 1, std::vector<double>(numEpochs + 1, 0.0));
    std::vector<std::vector<double>> S4C_S2(numSats + 1, std::vector<double>(numEpochs + 1, 0.0));

    for (int prn = 1; prn <= numSats; ++prn)
    {
        // (1) dB-Hz -> linear
        std::vector<double> s1_lin(numEpochs + 1, NaN), s2_lin(numEpochs + 1, NaN);
        for (int k = 1; k <= numEpochs; ++k) {
            double c1 = OBS.S1[prn][k], c2 = OBS.S2[prn][k];
            if (std::isfinite(c1) && c1 > 0.0) s1_lin[k] = std::pow(10.0, 0.1 * c1);
            if (std::isfinite(c2) && c2 > 0.0) s2_lin[k] = std::pow(10.0, 0.1 * c2);
        }

        // (2) Detrending: x[k] = s_lin[k] / mean( s_lin[k-n_trend ... k-1] )
        std::vector<double> x1(numEpochs + 1, NaN), x2(numEpochs + 1, NaN);
        double sumPrev1 = 0.0, sumPrev2 = 0.0;
        int    cntPrev1 = 0, cntPrev2 = 0;

        for (int k = 1; k <= numEpochs; ++k)
        {
            if (k - 1 >= 1) {
                if (std::isfinite(s1_lin[k - 1])) { sumPrev1 += s1_lin[k - 1]; cntPrev1++; }
                if (std::isfinite(s2_lin[k - 1])) { sumPrev2 += s2_lin[k - 1]; cntPrev2++; }
            }
            if (k - n_trend - 1 >= 1) {
                if (std::isfinite(s1_lin[k - n_trend - 1])) { sumPrev1 -= s1_lin[k - n_trend - 1]; cntPrev1--; }
                if (std::isfinite(s2_lin[k - n_trend - 1])) { sumPrev2 -= s2_lin[k - n_trend - 1]; cntPrev2--; }
            }

            if (cntPrev1 == n_trend && std::isfinite(s1_lin[k]) && sumPrev1 > 0.0)
                x1[k] = s1_lin[k] / (sumPrev1 / n_trend);
            if (cntPrev2 == n_trend && std::isfinite(s2_lin[k]) && sumPrev2 > 0.0)
                x2[k] = s2_lin[k] / (sumPrev2 / n_trend);
        }

        // (3) S4C: sliding window of length L_stat on x
        double sumX1 = 0.0, sumX1sq = 0.0; int cntX1 = 0;
        double sumX2 = 0.0, sumX2sq = 0.0; int cntX2 = 0;

        for (int k = 1; k <= numEpochs; ++k)
        {
            if (std::isfinite(x1[k])) { sumX1 += x1[k]; sumX1sq += x1[k] * x1[k]; cntX1++; }
            if (std::isfinite(x2[k])) { sumX2 += x2[k]; sumX2sq += x2[k] * x2[k]; cntX2++; }

            if (k - L_stat >= 1) {
                if (std::isfinite(x1[k - L_stat])) { sumX1 -= x1[k - L_stat]; sumX1sq -= x1[k - L_stat] * x1[k - L_stat]; cntX1--; }
                if (std::isfinite(x2[k - L_stat])) { sumX2 -= x2[k - L_stat]; sumX2sq -= x2[k - L_stat] * x2[k - L_stat]; cntX2--; }
            }

            // First computable epoch: k >= n_trend + L_stat, and the window has L_stat valid samples
            if (k >= n_trend + L_stat)
            {
                if (cntX1 == L_stat) {
                    double m1 = sumX1 / L_stat;
                    double v1 = std::fmax(sumX1sq / L_stat - m1 * m1, 0.0);
                    S4C_S1[prn][k] = (m1 > 0.0) ? std::sqrt(v1) / m1 : 0.0;
                }
                if (cntX2 == L_stat) {
                    double m2 = sumX2 / L_stat;
                    double v2 = std::fmax(sumX2sq / L_stat - m2 * m2, 0.0);
                    S4C_S2[prn][k] = (m2 > 0.0) ? std::sqrt(v2) / m2 : 0.0;
                }
            }
        }
    }

    namespace fs = std::filesystem;
    fs::path inpath(txt_output_path);
    fs::path dir = inpath.parent_path();
    std::error_code ec;
    fs::create_directories(dir, ec);

    std::string stem = inpath.stem().string();
    size_t us = stem.find('_');
    std::string station = (us == std::string::npos) ? stem : stem.substr(0, us);

    // hourWanted: "H01".."H24"
    std::ostringstream hs;
    hs << 'H' << std::setw(2) << std::setfill('0') << std::clamp(hourWanted, 1, 24);
    std::string hourTag = hs.str();

    fs::path pathS1 = dir / (station + "_" + systemTag + "_S4C_S1_" + hourTag + ".txt");
    fs::path pathS2 = dir / (station + "_" + systemTag + "_S4C_S2_" + hourTag + ".txt");

    std::ofstream fout1(pathS1.string()), fout2(pathS2.string());
    if (!fout1.is_open() || !fout2.is_open()) {
        std::cerr << "Failed to open S4C output files: "
            << pathS1 << " / " << pathS2 << "\n";
        return;
    }

    fout1 << std::fixed << std::setprecision(5);
    fout2 << std::fixed << std::setprecision(5);

    auto write_header = [&](std::ofstream& f) {
        f << std::setw(12) << "Epoch \\ PRN";
        for (int i = 1; i <= numSats; ++i) {
            char prn_buf[10];
            std::sprintf(prn_buf, "PRN%02d", i); // If you need to distinguish systems, change to G/C/E/R
            f << std::setw(11) << prn_buf;
        }
        f << "\n";
        };
    write_header(fout1);
    write_header(fout2);

    for (int j = 1; j <= numEpochs; ++j)
    {
        char epoch_buf[20];
        std::sprintf(epoch_buf, "Epoch %04d:", j);
        fout1 << std::setw(12) << epoch_buf;
        fout2 << std::setw(12) << epoch_buf;

        for (int prn = 1; prn <= numSats; ++prn)
        {
            fout1 << std::setw(11) << S4C_S1[prn][j];
            fout2 << std::setw(11) << S4C_S2[prn][j];
        }
        fout1 << "\n";
        fout2 << "\n";
    }

    fout1.close();
    fout2.close();
}


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

