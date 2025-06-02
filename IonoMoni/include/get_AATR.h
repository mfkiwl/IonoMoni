#pragma once
#include "sp3.h"
#include "obs.h"
#include <string>
#include "gcfg_ppp.h"

void AATR_G_prepro(obs& OBS, const std::string& stationName, sp3 SP3[], const std::string& txt_output_path, t_gcfg_ppp& gset);
void AATR_R_prepro(obs& OBS, const std::string& stationName, sp3 SP3[], const std::string& txt_output_path, t_gcfg_ppp& gset);
void AATR_E_prepro(obs& OBS, const std::string& stationName, sp3 SP3[], const std::string& txt_output_path, t_gcfg_ppp& gset);
void AATR_EX_prepro(obs& OBS, const std::string& stationName, sp3 SP3[], const std::string& txt_output_path, t_gcfg_ppp& gset);
void AATR_C_prepro(obs& OBS, const std::string& stationName, bool isC7IAllZero, sp3 SP3[], const std::string& txt_output_path, t_gcfg_ppp& gset);

