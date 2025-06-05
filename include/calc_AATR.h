#pragma once
#include "sp3.h"
#include "obs.h"
#include <string>
#include "gcfg_ppp.h"

void calc_aatr_GPS(obs& OBS, const std::string& stationName, sp3 SP3[], const std::string& txt_output_path, t_gcfg_ppp& gset);
void calc_aatr_GLO(obs& OBS, const std::string& stationName, sp3 SP3[], const std::string& txt_output_path, t_gcfg_ppp& gset);
void calc_aatr_GAL(obs& OBS, const std::string& stationName, sp3 SP3[], const std::string& txt_output_path, t_gcfg_ppp& gset);
void calc_aatr_BDS(obs& OBS, const std::string& stationName, bool isC7IAllZero, sp3 SP3[], const std::string& txt_output_path, t_gcfg_ppp& gset);

