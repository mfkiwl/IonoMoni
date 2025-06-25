#pragma once

#include <vector>
#include <string>
#include <set>
#include <memory>
#include "spdlog/spdlog.h"

bool readGPSDCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelliteDCB, double& receiverDCB,
    bool isC1WAllZero, std::set<int>& missingPRNs,
    std::shared_ptr<spdlog::logger> my_logger);

bool readBDSDCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelliteDCB, double& receiverDCB,
    bool isC7IAllZero, std::set<int>& missingPRNs,
    std::shared_ptr<spdlog::logger> my_logger);

bool readGLODCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelliteDCB, double& receiverDCB,
    std::set<int>& missingPRNs,
    std::shared_ptr<spdlog::logger> my_logger);

bool readGALDCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelliteDCB, double& receiverDCB,
    std::set<int>& missingPRNs,
    std::shared_ptr<spdlog::logger> my_logger,
    bool isC1QAllZero); 

