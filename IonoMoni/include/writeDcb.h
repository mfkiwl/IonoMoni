#pragma once
#include <vector>
#include <string>

void writeDCBToFile(const std::string& outputFilePath, const std::vector<double>& satelGFteDCB, double receiverDCB, const std::string& stationName);
