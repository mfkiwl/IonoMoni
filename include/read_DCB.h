#pragma once
#include <vector>
#include <string>

void readGPSDCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelGFteDCB, double& receiverDCB,
    bool isC1WAllZero);

void readBDSDCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelGFteDCB, double& receiverDCB,
    bool isC7IAllZero);

void readGLODCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelGFteDCB, double& receiverDCB);

void readGALDCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelGFteDCB, double& receiverDCB);

void readGALXDCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelGFteDCB, double& receiverDCB);
