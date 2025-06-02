#pragma once
#ifndef ROTI_CALC_H
#define ROTI_CALC_H

#include <string>
#include "obs.h"
#include "sp3.h"
#include "gcfg_ppp.h"

void calc_roti_GPS(const obs& OBS, const std::string& stationName, const sp3& SP3, bool isC1WAllZero, const std::string& txt_output_path, t_gcfg_ppp& gset);
void calc_roti_BDS(const obs& OBS, const std::string& stationName, const sp3& SP3, bool isC7IAllZero, const std::string& txt_output_path, t_gcfg_ppp& gset);
void calc_roti_GLO(const obs& OBS, const std::string& stationName, const sp3& SP3, const std::string& txt_output_path, t_gcfg_ppp& gset);
void calc_roti_GAL(const obs& OBS, const std::string& stationName, const sp3& SP3, const std::string& txt_output_path, t_gcfg_ppp& gset);
void calc_roti_GALX(const obs& OBS, const std::string& stationName, const sp3& SP3, const std::string& txt_output_path, t_gcfg_ppp& gset);


//void calc_S4C(const obs& OBS, int numSats, int numEpochs, int n, const std::string& outS1, const std::string& outS2);

#endif
