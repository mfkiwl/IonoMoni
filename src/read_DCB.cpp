#include "read_DCB.h"
#include <fstream>
#include <iostream>

void readGPSDCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelliteDCB, double& receiverDCB,
    bool isC1WAllZero) {
    std::ifstream file(dcbFilePath);
    if (!file.is_open()) {
        std::cerr << "Failed to open DCB file: " << dcbFilePath << std::endl;
        return;
    }

    satelliteDCB.assign(33, 0.0);
    receiverDCB = 0.0;

    std::string line;

    while (std::getline(file, line)) {
        if (line.length() < 91) continue;

        // satellite DCB
        if ((isC1WAllZero == false &&
            (line.substr(0, 7) == " DSB  G" || line.substr(0, 7) == " DCB  G") &&
            (line.substr(25, 13).find("C1W  C2W") != std::string::npos) &&
            line.substr(15, 4) == "    ") ||
            (isC1WAllZero &&
                (line.substr(0, 7) == " DSB  G" || line.substr(0, 7) == " DCB  G") &&
                (line.substr(25, 13).find("C1C  C2W") != std::string::npos) &&
                line.substr(15, 4) == "    ")) {

            int prn = std::stoi(line.substr(12, 2));
            double dcbValue = std::stod(line.substr(70, 22));
            if (prn >= 1 && prn <= 32) {
                satelliteDCB[prn] = dcbValue;
            }
        }

        // receiver DCB
        if ((isC1WAllZero == false &&
            (line.substr(0, 12) == " DSB  G    G" || line.substr(0, 12) == " DCB  G    G") &&
            (line.substr(25, 13).find("C1W  C2W") != std::string::npos)) ||
            (isC1WAllZero &&
                (line.substr(0, 12) == " DSB  G    G" || line.substr(0, 12) == " DCB  G    G") &&
                (line.substr(25, 13).find("C1C  C2W") != std::string::npos))) {
            std::string site = line.substr(15, 4);
            if (site == stationName) {
                receiverDCB = std::stod(line.substr(70, 22));
            }
        }
    }

    file.close();
}

void readBDSDCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelliteDCB, double& receiverDCB,
    bool isC7IAllZero) {
    std::ifstream file(dcbFilePath);
    if (!file.is_open()) {
        std::cerr << "Failed to open DCB file: " << dcbFilePath << std::endl;
        return;
    }

    satelliteDCB.assign(47, 0.0);
    receiverDCB = 0.0;
    std::string line;

    while (std::getline(file, line)) {
        if (line.length() < 91) continue;


        if ((isC7IAllZero == false &&
            (line.substr(0, 7) == " DSB  C" || line.substr(0, 7) == " DCB  C") &&
            (line.substr(25, 13).find("C2I  C7I") != std::string::npos) &&
            line.substr(15, 4) == "    ") ||
            (isC7IAllZero &&
                (line.substr(0, 7) == " DSB  C" || line.substr(0, 7) == " DCB  C") &&
                (line.substr(25, 13).find("C2I  C6I") != std::string::npos) &&
                line.substr(15, 4) == "    ")) {

            int prn = std::stoi(line.substr(12, 2));
            double dcbValue = std::stod(line.substr(70, 22));
            if (prn >= 1 && prn <= 46) {
                satelliteDCB[prn] = dcbValue;
            }
        }

        if ((isC7IAllZero == false &&
            (line.substr(0, 12) == " DSB  C    C" || line.substr(0, 12) == " DCB  C    C") &&
            (line.substr(25, 13).find("C2I  C7I") != std::string::npos)) ||
            (isC7IAllZero &&
                (line.substr(0, 12) == " DSB  C    C" || line.substr(0, 12) == " DCB  C    C") &&
                (line.substr(25, 13).find("C2I  C6I") != std::string::npos))) {
            std::string site = line.substr(15, 4);
            if (site == stationName) {
                receiverDCB = std::stod(line.substr(70, 22));
            }
        }
    }

    file.close();
}

void readGLODCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelliteDCB, double& receiverDCB) {
    std::ifstream file(dcbFilePath);
    if (!file.is_open()) {
        std::cerr << "Failed to open DCB file: " << dcbFilePath << std::endl;
        return;
    }

    satelliteDCB.assign(25, 0.0);
    receiverDCB = 0.0;
    std::string line;

    while (std::getline(file, line)) {
        if (line.length() < 91) continue;

 
        if ((line.substr(0, 7) == " DSB  R" || line.substr(0, 7) == " DCB  R") &&
            (line.substr(25, 13).find("C1P  C2P") != std::string::npos) &&
            line.substr(15, 4) == "    ") {

            int prn = std::stoi(line.substr(12, 2));
            double dcbValue = std::stod(line.substr(70, 22));
            if (prn >= 1 && prn <= 24) {
                satelliteDCB[prn] = dcbValue;
            }
        }

  
        if ((line.substr(0, 12) == " DSB  R    R" || line.substr(0, 12) == " DCB  R    R") &&
            (line.substr(25, 13).find("C1P  C2P") != std::string::npos)) {
            std::string site = line.substr(15, 4);
            if (site == stationName) {
                receiverDCB = std::stod(line.substr(70, 22));
            }
        }
    }

    file.close();
}

void readGALDCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelliteDCB, double& receiverDCB) {
    std::ifstream file(dcbFilePath);
    if (!file.is_open()) {
        std::cerr << "Failed to open DCB file: " << dcbFilePath << std::endl;
        return;
    }

    satelliteDCB.assign(37, 0.0);
    receiverDCB = 0.0;
    std::string line;

    while (std::getline(file, line)) {
        if (line.length() < 91) continue;

 
        if ((line.substr(0, 7) == " DSB  E" || line.substr(0, 7) == " DCB  E") &&
            (line.substr(25, 13).find("C1C  C5Q") != std::string::npos) &&
            line.substr(15, 4) == "    ") {

            int prn = std::stoi(line.substr(12, 2));
            double dcbValue = std::stod(line.substr(70, 22));
            if (prn >= 1 && prn <= 36) {
                satelliteDCB[prn] = dcbValue;
            }
        }


        if ((line.substr(0, 12) == " DSB  E    E" || line.substr(0, 12) == " DCB  E    E") &&
            (line.substr(25, 13).find("C1C  C5Q") != std::string::npos)) {
            std::string site = line.substr(15, 4);
            if (site == stationName) {
                receiverDCB = std::stod(line.substr(70, 22));
            }
        }
    }

    file.close();
}

void readGALXDCB(const std::string& dcbFilePath, const std::string& stationName,
    std::vector<double>& satelliteDCB, double& receiverDCB) {
    std::ifstream file(dcbFilePath);
    if (!file.is_open()) {
        std::cerr << "Failed to open DCB file: " << dcbFilePath << std::endl;
        return;
    }

    satelliteDCB.assign(37, 0.0);
    receiverDCB = 0.0;
    std::string line;

    while (std::getline(file, line)) {
        if (line.length() < 91) continue;

        
        if ((line.substr(0, 7) == " DSB  E" || line.substr(0, 7) == " DCB  E") &&
            (line.substr(25, 13).find("C1X  C5X") != std::string::npos) &&
            line.substr(15, 4) == "    ") {

            int prn = std::stoi(line.substr(12, 2));
            double dcbValue = std::stod(line.substr(70, 22));
            if (prn >= 1 && prn <= 36) {
                satelliteDCB[prn] = dcbValue;
            }
        }

        
        if ((line.substr(0, 12) == " DSB  E    E" || line.substr(0, 12) == " DCB  E    E") &&
            (line.substr(25, 13).find("C1X  C5X") != std::string::npos)) {
            std::string site = line.substr(15, 4);
            if (site == stationName) {
                receiverDCB = std::stod(line.substr(70, 22));
            }
        }
    }

    file.close();
}
