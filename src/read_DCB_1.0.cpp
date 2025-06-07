#include "read_DCB_1.0.h"
#include <fstream>
#include <regex>
#include <iostream>

void read_dcb_to_map(
    const std::string& dcbFile,
    std::map<std::string, double>& satelliteDCBMap,
    std::map<std::string, double>& receiverDCBMap)
{
    std::ifstream fin(dcbFile);
    if (!fin.is_open()) {
        std::cerr << "Failed to open DCB file: " << dcbFile << std::endl;
        return;
    }

    std::string line;
    // Regex pattern for satellite DCB entries
    std::regex sat_pattern(R"(DS[B|C]\s+\w+\s+([GREJC]\d{2})\s+([A-Z0-9]{3})\s+([A-Z0-9]{3})\s+[\d: ]+\s+[\d: ]+\s+\w+\s+([-\d.E]+))");
    // Regex pattern for receiver DCB entries
    std::regex rec_pattern(R"((DSB|DCB)\s+([A-Z])\s+[A-Z]\s+([A-Z0-9]{4})\s+([A-Z0-9]{3})\s+([A-Z0-9]{3})\s+[\d: ]+\s+[\d: ]+\s+\w+\s+([-\d.E]+))");

    std::smatch match;
    while (std::getline(fin, line)) {
        // Parse satellite DCB entry
        if (std::regex_search(line, match, sat_pattern)) {
            std::string prn = match[1].str();                  // e.g., G01
            std::string code1 = match[2].str();                // e.g., C1C
            std::string code2 = match[3].str();                // e.g., C1W
            std::string key = prn + "_" + code1 + code2;       // e.g., G01_C1CC1W
            double value = std::stod(match[4].str());           // DCB value (ns)
            satelliteDCBMap[key] = value;

            //std::cout << "[satelliteDCBMap] key = " << key << " value = " << value << " ns" << std::endl;



            continue;
        }

        // Parse receiver DCB entry
        if (std::regex_search(line, match, rec_pattern)) {
            std::string key = match[2].str() + "_" + match[3].str() + "_" + match[4].str() + match[5].str(); // e.g., G_BJFS_C1CC1W
            double value = std::stod(match[6].str());   
            receiverDCBMap[key] = value;              
            //std::cout << "[receiverDCBMap] key = " << key << " value = " << value << " ns" << std::endl;
            continue;
        }

    }

    // Uncomment for debugging:
    // for (const auto& [k, v] : receiverDCBMap)
    //     std::cout << "[receiverDCBMap] " << k << " = " << v << " ns\n";

    fin.close();
}


//double getDCBCorrection(const std::string& satSys, int prn,
//    const std::string& code1, const std::string& code2,
//    const std::map<std::string, double>& satelliteDCBMap,
//    const std::map<std::string, double>& receiverDCBMap,
//    const std::string& receiverName)
//{
//    char satKey[32];
//    sprintf(satKey, "%s%02d_%s%s", satSys.c_str(), prn, code1.c_str(), code2.c_str());
//    std::string recKey = satSys + "_" + receiverName + "_" + code1 + code2;
//
//    double satDCB = 0.0;
//    auto itSat = satelliteDCBMap.find(satKey);
//    if (itSat != satelliteDCBMap.end()) satDCB = itSat->second;
//
//    double recDCB = 0.0;
//    auto itRec = receiverDCBMap.find(recKey);
//    if (itRec != receiverDCBMap.end()) recDCB = itRec->second;
//
//    return satDCB - recDCB;
//}
