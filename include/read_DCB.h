#pragma once
#include <set>
#include <vector>
#include <string>
#include <iostream>


bool readGPSDCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelliteDCB, double& receiverDCB,
    bool isC1WAllZero, std::set<int>& missingPRNs,
    std::ostream& logger = std::cerr);

bool readBDSDCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelliteDCB, double& receiverDCB,
    bool isC7IAllZero, std::set<int>& missingPRNs,
    std::ostream& logger = std::cerr);

bool readGLODCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelliteDCB, double& receiverDCB,
    std::set<int>& missingPRNs,
    std::ostream& logger = std::cerr);

bool readGALDCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelliteDCB, double& receiverDCB,
    std::set<int>& missingPRNs,
    std::ostream& logger,
    bool isC1QAllZero);


