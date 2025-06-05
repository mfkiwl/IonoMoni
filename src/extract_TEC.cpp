#include "get_STEC.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include "read_DCB.h"
#include "extract_arcs.h"
#include "constants.h"
#include <string>
#include <filesystem>
#include <windows.h>
#include <fstream>     
#include <iomanip>     
#include "gcfg_ppp.h"
#include "detect_cycle_slip.h"
#include "calc_elevation.h"
#include "mapping_function.h"
#include <set>


void compute_MW_GF_wlAmb(double L1[], double L2[], double C1[], double C2[],
    double f1, double f2, double lambda_wl,
    double MW[], double GF[], double wlAmb[]) {
    // Compute Melbourne-Wübbena, geometry-free, and wide-lane ambiguity
    for (int j = 1; j <= 2880; j++) {
        MW[j] = lambda_wl * (L1[j] - L2[j]) - (f1 * C1[j] + f2 * C2[j]) / (f1 + f2);
        GF[j] = L1[j] - f1 * L2[j] / f2;
        wlAmb[j] = -MW[j];
    }
}

void clean_short_arcs(std::vector<std::vector<int>>& arcs, obs& OBS, int i, int arc_min_len) {
    for (int k = 0; k < (int)arcs.size(); ++k) {
        if (arcs[k][1] - arcs[k][0] < arc_min_len) {
            for (int e = arcs[k][0]; e <= arcs[k][1]; ++e) {
                OBS.C1[i][e] = OBS.C2[i][e] = OBS.L1[i][e] = OBS.L2[i][e] = 0;
            }
            arcs.erase(arcs.begin() + k);
            --k;
        }
    }
}

void carrier_smooth_pseudorange(double smoothedP4[], double L4[], int ARCS[][2], int arc_count) {
    // Carrier smoothing for pseudorange using phase difference
    for (int j = 1; j <= arc_count; j++) {
        int t = 2;
        for (int k = ARCS[j][0] + 1; k <= ARCS[j][1]; k++) {
            smoothedP4[k] = smoothedP4[k] / t + (smoothedP4[k - 1] + L4[k - 1] - L4[k]) * (t - 1) / t;
            t++;
        }
        for (int k = ARCS[j][0]; k <= ARCS[j][0] + 4; k++) {
            smoothedP4[k] = 0;
        }
    }
}

void split_and_filter_arcs(double* MW, int sat_index, int arc_min_len,
    int ARCS[][3000][2], obs& OBS) {
    // Find arcs and remove short segments
    int len = 2881;
    int num_arcs = 0;
    std::vector<std::vector<int>> arcs;

    extract_arcs(MW, len, arcs, num_arcs);

    for (int k = 0; k < num_arcs; ++k) {
        if (arcs[k][1] - arcs[k][0] < arc_min_len) {
            for (int e = arcs[k][0]; e <= arcs[k][1]; ++e) {
                OBS.C1[sat_index][e] = OBS.C2[sat_index][e] = 0;
                OBS.L1[sat_index][e] = OBS.L2[sat_index][e] = 0;
            }
            arcs.erase(arcs.begin() + k);
            --k;
            --num_arcs;
        }
    }

    ARCS[sat_index][0][0] = num_arcs;
    for (int k = 0; k < num_arcs; ++k) {
        ARCS[sat_index][k + 1][0] = arcs[k][0];
        ARCS[sat_index][k + 1][1] = arcs[k][1];
    }
}

void apply_smoothing_filter(double smoothed[][3000], double phase[][3000],
    int ARCS[][3000][2], int sat_index) {
    // Apply carrier smoothing for each valid arc
    for (int j = 1; j <= ARCS[sat_index][0][0]; j++) {
        int t = 2;
        for (int k = ARCS[sat_index][j][0] + 1; k <= ARCS[sat_index][j][1]; k++) {
            smoothed[sat_index][k] = smoothed[sat_index][k] / t +
                (smoothed[sat_index][k - 1] + phase[sat_index][k - 1] - phase[sat_index][k]) * (t - 1) / t;
            t++;
        }
        for (int k = ARCS[sat_index][j][0]; k <= ARCS[sat_index][j][0] + 4; k++) {
            smoothed[sat_index][k] = 0;
        }
    }
}

void apply_dcb_correction_to_matrix(double matrix[][3000], int sat_count,
    const std::vector<double>& sat_dcb, double rec_dcb) {
    // DCB correction for all satellites and epochs
    for (int i = 1; i <= sat_count; i++) {
        for (int j = 1; j <= 2880; j++) {
            if (matrix[i][j] != 0) {
                matrix[i][j] -= c * sat_dcb[i] * 1e-9;
                matrix[i][j] -= c * rec_dcb * 1e-9;
            }
        }
    }
}

void output_smoothed_stec_txt(const std::string& filepath, const std::string& prefix,
    int max_prn, double smoothde_P4[][3000], double f1, double f2) {
    std::ofstream fout(filepath);
    if (!fout.is_open()) {
        std::cerr << "Failed to open STEC output file: " << filepath << std::endl;
        return;
    }

    fout << std::setw(12) << "Epoch \\ PRN";
    for (int i = 1; i <= max_prn; ++i) {
        char prn_buf[10];
        sprintf(prn_buf, "%s%02d", prefix.c_str(), i);
        fout << std::setw(11) << prn_buf;
    }
    fout << "\n";

    for (int j = 1; j <= 2880; ++j) {
        char epoch_buf[20];
        sprintf(epoch_buf, "Epoch %04d:", j);
        fout << std::setw(12) << epoch_buf;

        for (int i = 1; i <= max_prn; ++i) {
            double result = smoothde_P4[i][j] * (f1 * f1 * f2 * f2) /
                (IONO_COEFF * (f2 * f2 - f1 * f1)) / 1e16;

            if (fabs(result) < 1e-8) result = 0.0;

            fout << std::setw(11) << std::fixed << std::setprecision(5) << result;
        }
        fout << "\n";
    }

    fout.close();
}

void output_smoothed_vtec_txt(
    const std::string& filepath,
    const std::string& prefix,
    int max_prn,
    double smoothde_P4[][3000],
    double f1,
    double f2,
    const obs& OBS,
    const sp3* SP3,
    int mf_type,
    std::function<void(const sp3*, int prn, int epoch, double& x, double& y, double& z)> get_sp3_coord_func
) {
    std::ofstream fout(filepath);
    if (!fout.is_open()) {
        std::cerr << "Failed to open output file: " << filepath << std::endl;
        return;
    }

    fout << std::setw(12) << "Epoch \\ PRN";
    for (int i = 1; i <= max_prn; ++i) {
        char prn_buf[10];
        sprintf(prn_buf, "%s%02d", prefix.c_str(), i);
        fout << std::setw(11) << prn_buf;
    }
    fout << "\n";

    for (int j = 1; j <= 2880; ++j) {
        char epoch_buf[20];
        sprintf(epoch_buf, "Epoch %04d:", j);
        fout << std::setw(12) << epoch_buf;

        for (int i = 1; i <= max_prn; ++i) {
            double stec = smoothde_P4[i][j] * (f1 * f1 * f2 * f2) /
                (IONO_COEFF * (f2 * f2 - f1 * f1)) / 1e16;
            if (fabs(stec) < 1e-8) stec = 0.0;

            double x, y, z;
            get_sp3_coord_func(SP3, i, j, x, y, z);
            x *= 1000; y *= 1000; z *= 1000;

            double sx = OBS.X, sy = OBS.Y, sz = OBS.Z;
            double E_deg, A_deg;
            calc_elevation(x, y, z, sx, sy, sz, E_deg, A_deg);
            double elev_rad = E_deg * PI / 180.0;
            double mf = get_mapping_function(elev_rad, mf_type);

            double vtec = (mf > 0.01) ? (stec / mf) : 0.0;
            if (fabs(vtec) < 1e-8) vtec = 0.0;

            fout << std::setw(11) << std::fixed << std::setprecision(5) << vtec;
        }
        fout << "\n";
    }

    fout.close();
}

void output_smoothed_stec_txt_GLO(const std::string& filepath,
    int max_prn, double stec[][3000],
    const double f1[], const double f2[]) {
    std::ofstream fout(filepath);
    if (!fout.is_open()) {
        std::cerr << "Failed to open GLO STEC output file: " << filepath << std::endl;
        return;
    }

    fout << std::setw(12) << "Epoch \\ PRN";
    for (int i = 1; i <= max_prn; ++i) {
        char prn_buf[10];
        sprintf(prn_buf, "R%02d", i);
        fout << std::setw(11) << prn_buf;
    }
    fout << "\n";

    for (int j = 1; j <= 2880; ++j) {
        char epoch_buf[20];
        sprintf(epoch_buf, "Epoch %04d:", j);
        fout << std::setw(12) << epoch_buf;

        for (int i = 1; i <= max_prn; ++i) {
            double result = stec[i][j] * (f1[i] * f1[i] * f2[i] * f2[i])
                / (IONO_COEFF * (f2[i] * f2[i] - f1[i] * f1[i])) / 1e16;

            if (fabs(result) < 1e-8) result = 0.0;

            fout << std::setw(11) << std::fixed << std::setprecision(5) << result;
        }
        fout << "\n";
    }

    fout.close();
}

void output_smoothed_vtec_txt_GLO(
    const std::string& filepath,
    int max_prn,
    double smoothde_P4[][3000],
    const double f1[],
    const double f2[],
    const obs& OBS,
    const sp3* SP3,
    int mf_type
) {
    std::ofstream fout(filepath);
    if (!fout.is_open()) {
        std::cerr << "Failed to open output file: " << filepath << std::endl;
        return;
    }

    fout << std::setw(12) << "Epoch \\ PRN";
    for (int i = 1; i <= max_prn; ++i) {
        char prn_buf[10];
        sprintf(prn_buf, "R%02d", i);
        fout << std::setw(11) << prn_buf;
    }
    fout << "\n";

    for (int j = 1; j <= 2880; ++j) {
        char epoch_buf[20];
        sprintf(epoch_buf, "Epoch %04d:", j);
        fout << std::setw(12) << epoch_buf;

        for (int i = 1; i <= max_prn; ++i) {
            double f1i = f1[i];
            double f2i = f2[i];

            double stec = smoothde_P4[i][j] * (f1i * f1i * f2i * f2i) /
                (IONO_COEFF * (f2i * f2i - f1i * f1i)) / 1e16;
            if (fabs(stec) < 1e-8) stec = 0.0;

            double x = SP3[1].RX[i][j] * 1000;
            double y = SP3[1].RY[i][j] * 1000;
            double z = SP3[1].RZ[i][j] * 1000;

            double sx = OBS.X, sy = OBS.Y, sz = OBS.Z;
            double E_deg, A_deg;
            calc_elevation(x, y, z, sx, sy, sz, E_deg, A_deg);
            double elev_rad = E_deg * PI / 180.0;
            double mf = get_mapping_function(elev_rad, mf_type);

            double vtec = (mf > 0.01) ? (stec / mf) : 0.0;
            if (fabs(vtec) < 1e-8) vtec = 0.0;

            fout << std::setw(11) << std::fixed << std::setprecision(5) << vtec;
        }

        fout << "\n";
    }

    fout.close();
}


void G_prepro(
    obs& OBS,
    const std::string& stationName,
    bool isC1WAllZero,
    const std::string& dcbFilePath,
    const std::string& txt_output_path,
    t_gcfg_ppp& gset,
    const sp3* SP3,
    int mf_type
) {
    // -------- Main pre-processing workflow for GPS STEC extraction --------
    double f1 = 1575.42e6, f2 = 1227.6e6;
    double lambda_wl = c / (f1 - f2);
    double MW[33][2881], GF[33][2881], wlAmb[33][2881];
    double P4[33][3000], L4[33][3000], smoothedP4[33][3000];
    double cs_epoch[33][2881] = { 0 };
    int ARCS[33][3000][2];

    for (int i = 1; i <= 32; i++) {
        // Compute MW, GF, and wide-lane ambiguity
        compute_MW_GF_wlAmb(OBS.L1[i], OBS.L2[i], OBS.C1[i], OBS.C2[i], f1, f2, lambda_wl,
            MW[i], GF[i], wlAmb[i]);
        // Detect arcs and filter out short ones
        split_and_filter_arcs(MW[i], i, gset.arc_min_length(), ARCS, OBS);
        // Detect and mark cycle slips
        detect_cycle_slip(wlAmb[i], GF[i], OBS, i, ARCS, gset.arc_min_length(), cs_epoch);

        // Calculate geometry-free pseudorange and phase, initialize smoothing array
        for (int j = 1; j <= 2880; ++j) {
            P4[i][j] = OBS.C1[i][j] - OBS.C2[i][j];
            L4[i][j] = (c / f1) * OBS.L1[i][j] - (c / f2) * OBS.L2[i][j];
            smoothedP4[i][j] = P4[i][j];
        }
        // Carrier phase smoothing
        apply_smoothing_filter(smoothedP4, L4, ARCS, i);
    }

    // Apply DCB correction to smoothed pseudorange
    std::vector<double> satelDCB;
    double recDCB;
    std::set<int> missingPRNs;
    if (!readGPSDCB(dcbFilePath, stationName, satelDCB, recDCB, isC1WAllZero, missingPRNs, std::cerr)) {
        throw std::runtime_error("No DCB record found for station: " + stationName);
    }
    apply_dcb_correction_to_matrix(smoothedP4, 32, satelDCB, recDCB);

    for (int prn : missingPRNs) {
        for (int j = 1; j <= 2880; ++j) {
            smoothedP4[prn][j] = 0.0;
        }
    }

    std::string out_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GPS.txt";
    std::string vtec_out_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_VTEC_GPS.txt";
    output_smoothed_stec_txt(out_path, "G", 32, smoothedP4, f1, f2);
    output_smoothed_vtec_txt(
        vtec_out_path, "G", 32, smoothedP4, f1, f2, OBS, SP3, mf_type,
        [](const sp3* sp3_ptr, int prn, int epoch, double& x, double& y, double& z) {
            x = sp3_ptr[1].X[prn][epoch];
            y = sp3_ptr[1].Y[prn][epoch];
            z = sp3_ptr[1].Z[prn][epoch];
        });
}


void R_prepro(
    obs& OBS,
    const std::string& stationName,
    const std::string& dcbFilePath,
    const std::string& txt_output_path,
    t_gcfg_ppp& gset,
    const sp3* SP3,
    int mf_type
) {
    int Fre[25] = { 99, 1, -4, 5, 6, 1, -4, 5, 6, -2, -7, 0, -1, -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2 };
    double f1[25], f2[25], lambda_wl[25];
    double MW[25][2881], GF[25][2881], wlAmb[25][2881];
    double P4[25][3000], L4[25][3000], smoothedP4[25][3000];
    double cs_epoch[25][2881] = { 0 };
    int ARCS[25][3000][2];

    for (int i = 1; i <= 24; i++) {
        f1[i] = (1602.0 + Fre[i] * 0.5625) * 1e6;
        f2[i] = (1246.0 + Fre[i] * 0.4375) * 1e6;
        lambda_wl[i] = c / (f1[i] - f2[i]);

        compute_MW_GF_wlAmb(OBS.L1[i], OBS.L2[i], OBS.C1[i], OBS.C2[i],
            f1[i], f2[i], lambda_wl[i],
            MW[i], GF[i], wlAmb[i]);

        split_and_filter_arcs(MW[i], i, gset.arc_min_length(), ARCS, OBS);
        detect_cycle_slip(wlAmb[i], GF[i], OBS, i, ARCS, gset.arc_min_length(), cs_epoch);

        for (int j = 1; j <= 2880; j++) {
            P4[i][j] = OBS.C1[i][j] - OBS.C2[i][j];
            L4[i][j] = (c / f1[i]) * OBS.L1[i][j] - (c / f2[i]) * OBS.L2[i][j];
            smoothedP4[i][j] = P4[i][j];
        }
        apply_smoothing_filter(smoothedP4, L4, ARCS, i);
    }

    std::vector<double> satelDCB;
    double recDCB;
    std::set<int> missingPRNs;
    if (!readGLODCB(dcbFilePath, stationName, satelDCB, recDCB, missingPRNs)) {
        throw std::runtime_error("No DCB record found for station: " + stationName);
    }
    apply_dcb_correction_to_matrix(smoothedP4, 24, satelDCB, recDCB);

    for (int prn : missingPRNs) {
        for (int j = 1; j <= 2880; ++j) {
            smoothedP4[prn][j] = 0.0;
        }
    }

    std::string out_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GLO.txt";
    std::string vtec_out_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_VTEC_GLO.txt";
    output_smoothed_stec_txt_GLO(out_path, 24, smoothedP4, f1, f2);
    output_smoothed_vtec_txt_GLO(vtec_out_path, 24, smoothedP4, f1, f2, OBS, SP3, mf_type);
}


void E_prepro(
    obs& OBS,
    const std::string& stationName,
    const std::string& dcbFilePath,
    const std::string& txt_output_path,
    t_gcfg_ppp& gset,
    const sp3* SP3,
    int mf_type)
{
    double f1 = 1575.42e6, f2 = 1176.45e6;
    double lambda_wl = c / (f1 - f2);
    double MW[37][2881], GF[37][2881], wlAmb[37][2881];
    double P4[37][3000], L4[37][3000], smoothedP4[37][3000];
    double cs_epoch[37][2881] = { 0 };
    int ARCS[37][3000][2];

    for (int i = 1; i <= 36; i++) {
        compute_MW_GF_wlAmb(OBS.L1[i], OBS.L2[i], OBS.C1[i], OBS.C2[i], f1, f2, lambda_wl,
            MW[i], GF[i], wlAmb[i]);
        split_and_filter_arcs(MW[i], i, gset.arc_min_length(), ARCS, OBS);
        detect_cycle_slip(wlAmb[i], GF[i], OBS, i, ARCS, gset.arc_min_length(), cs_epoch);

        for (int j = 1; j <= 2880; ++j) {
            P4[i][j] = OBS.C1[i][j] - OBS.C2[i][j];
            L4[i][j] = (c / f1) * OBS.L1[i][j] - (c / f2) * OBS.L2[i][j];
            smoothedP4[i][j] = P4[i][j];
        }
        apply_smoothing_filter(smoothedP4, L4, ARCS, i);
    }

    std::vector<double> satelDCB;
    double recDCB;
    std::set<int> missingPRNs;
    bool isC1QAllZero;
    if (!readGALDCB(dcbFilePath, stationName, satelDCB, recDCB, missingPRNs, std::cerr, isC1QAllZero)) {
        throw std::runtime_error("No DCB record found for station: " + stationName);
    }
    apply_dcb_correction_to_matrix(smoothedP4, 36, satelDCB, recDCB);

    for (int prn : missingPRNs) {
        for (int j = 1; j <= 2880; ++j) {
            smoothedP4[prn][j] = 0.0;
        }
    }

    std::string out_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GAL.txt";
    std::string vtec_out_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_VTEC_GAL.txt";
    output_smoothed_stec_txt(out_path, "E", 36, smoothedP4, f1, f2);
    output_smoothed_vtec_txt(
        vtec_out_path, "E", 36, smoothedP4, f1, f2, OBS, SP3, mf_type,
        [](const sp3* sp3_ptr, int prn, int epoch, double& x, double& y, double& z) {
            x = sp3_ptr[1].EX[prn][epoch];
            y = sp3_ptr[1].EY[prn][epoch];
            z = sp3_ptr[1].EZ[prn][epoch];
        });
}



void C_prepro(obs& OBS, const std::string& stationName, bool isC7IAllZero, const std::string& dcbFilePath, const std::string& txt_output_path, t_gcfg_ppp& gset, const sp3* SP3, int mf_type) {
    double f1 = 1561.098e6;
    double f2 = isC7IAllZero ? 1268.52e6 : 1207.14e6;
    double lambda_wl = c / (f1 - f2);
    double MW[47][2881], GF[47][2881], wlAmb[47][2881];
    double P4[47][3000], L4[47][3000], smoothedP4[47][3000];
    double cs_epoch[47][2881] = { 0 };
    int ARCS[47][3000][2];

    for (int i = 1; i <= 46; i++) {
        compute_MW_GF_wlAmb(OBS.L1[i], OBS.L2[i], OBS.C1[i], OBS.C2[i], f1, f2, lambda_wl,
            MW[i], GF[i], wlAmb[i]);
        split_and_filter_arcs(MW[i], i, gset.arc_min_length(), ARCS, OBS);
        detect_cycle_slip(wlAmb[i], GF[i], OBS, i, ARCS, gset.arc_min_length(), cs_epoch);

        for (int j = 1; j <= 2880; ++j) {
            P4[i][j] = OBS.C1[i][j] - OBS.C2[i][j];
            L4[i][j] = (c / f1) * OBS.L1[i][j] - (c / f2) * OBS.L2[i][j];
            smoothedP4[i][j] = P4[i][j];
        }
        apply_smoothing_filter(smoothedP4, L4, ARCS, i);
    }

    std::vector<double> satelDCB;
    double recDCB;
    std::set<int> missingPRNs;
    if (!readBDSDCB(dcbFilePath, stationName, satelDCB, recDCB, isC7IAllZero, missingPRNs)) {
        throw std::runtime_error("No DCB record found for station: " + stationName);
    }
    apply_dcb_correction_to_matrix(smoothedP4, 46, satelDCB, recDCB);

    for (int prn : missingPRNs) {
        for (int j = 1; j <= 2880; ++j) {
            smoothedP4[prn][j] = 0.0;
        }
    }

    std::string out_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_BDS.txt";
    std::string vtec_out_path = txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_VTEC_BDS.txt";
    output_smoothed_stec_txt(out_path, "C", 46, smoothedP4, f1, f2);
    output_smoothed_vtec_txt(
        vtec_out_path, "C", 46, smoothedP4, f1, f2, OBS, SP3, mf_type,
        [](const sp3* sp3_ptr, int prn, int epoch, double& x, double& y, double& z) {
            x = sp3_ptr[1].CX[prn][epoch];
            y = sp3_ptr[1].CY[prn][epoch];
            z = sp3_ptr[1].CZ[prn][epoch];
        });
}

