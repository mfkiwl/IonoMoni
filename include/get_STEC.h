#pragma once
#include "sp3.h"
#include "obs.h"
#include <vector>
#include <string>
#include "gcfg_ppp.h"

// 原始主函数接口
void G_prepro(obs& OBS, const std::string& stationName, bool isC1WAllZero, const std::string& dcbFilePath, const std::string& txt_output_path, t_gcfg_ppp& gset, const sp3* SP3, int mf_type);
void R_prepro(obs& OBS, const std::string& stationName, const std::string& dcbFilePath, const std::string& txt_output_path, t_gcfg_ppp& gset, const sp3* SP3, int mf_type);
void E_prepro(obs& OBS, const std::string& stationName, const std::string& dcbFilePath, const std::string& txt_output_path, t_gcfg_ppp& gset, const sp3* SP3, int mf_type);
void C_prepro(obs& OBS, const std::string& stationName, bool isC7IAllZero, const std::string& dcbFilePath, const std::string& txt_output_path, t_gcfg_ppp& gset, const sp3* SP3, int mf_type);


void compute_MW_GF_wlAmb(double L1[], double L2[], double C1[], double C2[],
    double f1, double f2, double lambda_wl,
    double MW[], double GF[], double wlAmb[]);

void split_and_filter_arcs(double* MW, int sat_index, int arc_min_len,
    int ARCS[][3000][2], obs& OBS);

