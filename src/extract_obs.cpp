#include "extract_obs.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <set>

using namespace std;
using namespace gnut;

namespace gnut {
    void extract_GPS_obs(t_gallobs* gobs, const std::string& station,
        std::vector<std::vector<double>>& GPS_C1,
        std::vector<std::vector<double>>& GPS_C2,
        std::vector<std::vector<double>>& GPS_L1,
        std::vector<std::vector<double>>& GPS_L2,
        std::vector<t_gtime>& epochs,
        std::vector<std::string>& GPS_sats,
        bool& isC1WAllZero,
        std::shared_ptr<spdlog::logger> logger, obs& OBS)
    {
        // Define observation code priorities for GPS
        std::vector<std::string> C1W_C1C_priority = { "C1W", "C1C" };
        std::vector<std::string> C1_old_priority = { "P1", "C1" };
        std::vector<std::string> GPS_C2_priority = { "C2W", "P2", "C2" };
        std::set<std::string> GPS_L1_types = { "L1C", "L1" };
        std::set<std::string> GPS_L2_types = { "L2W", "L2" };

        const size_t num_epochs = 2881;
        const size_t num_sats = 33;

        GPS_C1.assign(num_epochs, std::vector<double>(num_sats, 0.0));
        GPS_C2.assign(num_epochs, std::vector<double>(num_sats, 0.0));
        GPS_L1.assign(num_epochs, std::vector<double>(num_sats, 0.0));
        GPS_L2.assign(num_epochs, std::vector<double>(num_sats, 0.0));

        std::vector<t_gtime> all_epochs = gobs->epochs(station);

        GPS_sats.clear();
        for (int s = 1; s <= 32; ++s)
        {
            std::ostringstream oss;
            oss << "G" << std::setw(2) << std::setfill('0') << s;
            GPS_sats.push_back(oss.str());
        }

        // Variables to store the chosen observation codes for this dataset
        std::string chosen_C1 = "";
        std::string chosen_C2 = "";
        std::string chosen_L1 = "";
        std::string chosen_L2 = "";

        // Select observation types (C1, C2, L1, L2) based on priority and data availability
        for (size_t i = 0; i < all_epochs.size(); ++i)
        {
            std::vector<t_gsatdata> obsvec = gobs->obs(station, all_epochs[i]);
            for (const auto& obs : obsvec)
            {
                if (obs.sat()[0] != 'G') continue; // Only process GPS satellites
                const auto& GFst = obs.obs();

                // Select C1 type
                if (chosen_C1.empty())
                {
                    for (const auto& type : C1W_C1C_priority)
                    {
                        if (std::find(GFst.begin(), GFst.end(), str2gobs(type)) != GFst.end())
                        {
                            chosen_C1 = type;
                            break;
                        }
                    }
                    // Try fallback to old P1/C1 if no C1W/C1C found
                    if (chosen_C1.empty())
                    {
                        bool found = false;
                        for (const auto& type : C1_old_priority)
                        {
                            for (size_t ep = 0; ep < all_epochs.size(); ++ep)
                            {
                                std::vector<t_gsatdata> obsvec = gobs->obs(station, all_epochs[ep]);
                                for (const auto& ob : obsvec)
                                {
                                    if (ob.sat()[0] != 'G') continue;
                                    const auto& oGFst = ob.obs();

                                    if (std::find(oGFst.begin(), oGFst.end(), str2gobs(type)) != oGFst.end())
                                    {
                                        chosen_C1 = type;
                                        found = true;
                                        break;
                                    }
                                }
                                if (found) break;
                            }
                            if (found) break;
                        }
                    }
                }

                // Select C2 type
                if (chosen_C2.empty())
                {
                    bool found = false;
                    for (const auto& type : GPS_C2_priority)
                    {
                        for (size_t ep = 0; ep < all_epochs.size(); ++ep)
                        {
                            std::vector<t_gsatdata> obsvec = gobs->obs(station, all_epochs[ep]);
                            for (const auto& ob : obsvec)
                            {
                                if (ob.sat()[0] != 'G') continue;
                                const auto& oGFst = ob.obs();

                                if (std::find(oGFst.begin(), oGFst.end(), str2gobs(type)) != oGFst.end())
                                {
                                    chosen_C2 = type;
                                    found = true;
                                    break;
                                }
                            }
                            if (found) break;
                        }
                        if (found) break;
                    }
                }

                // Select L1 type
                if (chosen_L1.empty())
                {
                    for (const auto& type : GPS_L1_types)
                    {
                        if (std::find(GFst.begin(), GFst.end(), str2gobs(type)) != GFst.end())
                        {
                            chosen_L1 = type;
                            break;
                        }
                    }
                }

                // Select L2 type
                if (chosen_L2.empty())
                {
                    for (const auto& type : GPS_L2_types)
                    {
                        if (std::find(GFst.begin(), GFst.end(), str2gobs(type)) != GFst.end())
                        {
                            chosen_L2 = type;
                            break;
                        }
                    }
                }

                // If all four types have been selected, no need to continue
                if (!chosen_C1.empty() && !chosen_C2.empty() && !chosen_L1.empty() && !chosen_L2.empty())
                    break;
            }
            if (!chosen_C1.empty() && !chosen_C2.empty() && !chosen_L1.empty() && !chosen_L2.empty())
                break;
        }

        // Set isC1WAllZero according to selected C1 type
        if (chosen_C1 == "C1W" || chosen_C1 == "P1")
            isC1WAllZero = false;
        else
            isC1WAllZero = true;

        logger->info("[{}] Selected GPS types: C1={}, C2={}, L1={}, L2={}", station, chosen_C1, chosen_C2, chosen_L1, chosen_L2);

        if (chosen_C1.empty()) logger->warn("[{}] No C1 observation type found!", station);
        if (chosen_C2.empty()) logger->warn("[{}] No C2 observation type found!", station);
        if (chosen_L1.empty()) logger->warn("[{}] No L1 observation type found!", station);
        if (chosen_L2.empty()) logger->warn("[{}] No L2 observation type found!", station);

        // Fill C1/C2/L1/L2 arrays for all satellites and epochs
        for (size_t i = 0; i < std::min<size_t>(all_epochs.size(), num_epochs - 1); ++i)
        {
            const t_gtime& epo = all_epochs[i];
            std::vector<t_gsatdata> obsvec = gobs->obs(station, epo);

            for (const auto& obs : obsvec)
            {
                std::string prn = obs.sat();
                if (prn[0] != 'G') continue;

                // Find satellite index in GPS_sats vector
                auto it = std::find(GPS_sats.begin(), GPS_sats.end(), prn);
                if (it == GPS_sats.end()) continue;
                size_t j = std::distance(GPS_sats.begin(), it) + 1;
                if (j >= num_sats) continue;

                const auto& GFst = obs.obs();

                if (!chosen_C1.empty() && std::find(GFst.begin(), GFst.end(), str2gobs(chosen_C1)) != GFst.end())
                    GPS_C1[i + 1][j] = obs.getobs(str2gobs(chosen_C1));

                if (!chosen_C2.empty() && std::find(GFst.begin(), GFst.end(), str2gobs(chosen_C2)) != GFst.end())
                    GPS_C2[i + 1][j] = obs.getobs(str2gobs(chosen_C2));
 
                if (!chosen_L1.empty() && std::find(GFst.begin(), GFst.end(), str2gobs(chosen_L1)) != GFst.end())
                    GPS_L1[i + 1][j] = obs.getobs(str2gobs(chosen_L1));

                if (!chosen_L2.empty() && std::find(GFst.begin(), GFst.end(), str2gobs(chosen_L2)) != GFst.end())
                    GPS_L2[i + 1][j] = obs.getobs(str2gobs(chosen_L2));
            }

            epochs.push_back(epo);
        }

        for (int i = 1; i <= 32; i++) {
            for (int j = 1; j <= 2880; j++) {
                OBS.C1[i][j] = GPS_C1[j][i];
                OBS.C2[i][j] = GPS_C2[j][i];
                OBS.L1[i][j] = GPS_L1[j][i];
                OBS.L2[i][j] = GPS_L2[j][i];
            }
        }
    }
    void extract_BDS_obs(
        t_gallobs* gobs, const std::string& station,
        std::vector<std::vector<double>>& BDS_C2,
        std::vector<std::vector<double>>& BDS_C6,
        std::vector<std::vector<double>>& BDS_C7,
        std::vector<std::vector<double>>& BDS_L2,
        std::vector<std::vector<double>>& BDS_L6,
        std::vector<std::vector<double>>& BDS_L7,
        std::vector<t_gtime>& epochs,
        std::vector<std::string>& BDS_sats,
        obs& OBS,
        std::shared_ptr<spdlog::logger> logger,
        bool& isC7IAllZero)
    {
        const size_t num_epochs = 2881;
        const size_t num_sats = 47;

        BDS_C2.assign(num_epochs, std::vector<double>(num_sats, 0.0));
        BDS_C6.assign(num_epochs, std::vector<double>(num_sats, 0.0));
        BDS_C7.assign(num_epochs, std::vector<double>(num_sats, 0.0));
        BDS_L2.assign(num_epochs, std::vector<double>(num_sats, 0.0));
        BDS_L6.assign(num_epochs, std::vector<double>(num_sats, 0.0));
        BDS_L7.assign(num_epochs, std::vector<double>(num_sats, 0.0));

        std::vector<t_gtime> all_epochs = gobs->epochs(station);

        BDS_sats.clear();
        for (int s = 1; s <= 46; ++s) {
            std::ostringstream oss;
            oss << "C" << std::setw(2) << std::setfill('0') << s;
            BDS_sats.push_back(oss.str());
        }

        std::vector<std::string> C1_priority = { "C2I", "C2" }; 
        std::vector<std::string> C2_priority = { "C6I", "C6I", "C7I", "C7" };
        std::vector<std::string> L1_priority = { "L2I", "L2" }; 
        std::vector<std::string> L2_priority = { "L6I", "L6", "L7I", "L7" };

        std::string chosen_C1 = "", chosen_C2 = "", chosen_L1 = "", chosen_L2 = "";

        for (size_t i = 0; i < std::min<size_t>(all_epochs.size(), num_epochs - 1); ++i) {
            std::vector<t_gsatdata> obsvec = gobs->obs(station, all_epochs[i]);
            for (const auto& obs : obsvec) {
                std::string prn = obs.sat();
                if (prn[0] != 'C') continue;

                const auto& GFst = obs.obs();

                if (chosen_C1.empty()) {
                    for (const auto& type : C1_priority) {
                        if (std::find(GFst.begin(), GFst.end(), str2gobs(type)) != GFst.end()) {
                            chosen_C1 = type;
                            break;
                        }
                    }
                }

                if (chosen_C2.empty()) {
                    for (const auto& type : C2_priority) {
                        if (std::find(GFst.begin(), GFst.end(), str2gobs(type)) != GFst.end()) {
                            chosen_C2 = type;
                            break;
                        }
                    }
                }
  
                if (chosen_L1.empty()) {
                    for (const auto& type : L1_priority) {
                        if (std::find(GFst.begin(), GFst.end(), str2gobs(type)) != GFst.end()) {
                            chosen_L1 = type;
                            break;
                        }
                    }
                }

                if (chosen_L2.empty()) {
                    for (const auto& type : L2_priority) {
                        if (std::find(GFst.begin(), GFst.end(), str2gobs(type)) != GFst.end()) {
                            chosen_L2 = type;
                            break;
                        }
                    }
                }
                if (!chosen_C1.empty() && !chosen_C2.empty() && !chosen_L1.empty() && !chosen_L2.empty())
                    break;
            }
            if (!chosen_C1.empty() && !chosen_C2.empty() && !chosen_L1.empty() && !chosen_L2.empty())
                break;
        }

        logger->info("[{}] Selected BDS types: C1={}, C2={}, L1={}, L2={}",
            station, chosen_C1, chosen_C2, chosen_L1, chosen_L2);

        if (chosen_C1.empty()) logger->warn("[{}] No C1 observation type found!", station);
        if (chosen_C2.empty()) logger->warn("[{}] No C2 observation type found!", station);
        if (chosen_L1.empty()) logger->warn("[{}] No L1 observation type found!", station);
        if (chosen_L2.empty()) logger->warn("[{}] No L2 observation type found!", station);

        for (size_t i = 0; i < std::min<size_t>(all_epochs.size(), num_epochs - 1); ++i) {
            const t_gtime& epo = all_epochs[i];
            std::vector<t_gsatdata> obsvec = gobs->obs(station, epo);

            for (const auto& obs : obsvec) {
                std::string prn = obs.sat();
                if (prn[0] != 'C') continue;

                auto it = std::find(BDS_sats.begin(), BDS_sats.end(), prn);
                if (it == BDS_sats.end()) continue;
                size_t j = std::distance(BDS_sats.begin(), it) + 1;
                if (j >= num_sats) continue;

                const auto& GFst = obs.obs();

                if (!chosen_C1.empty() && std::find(GFst.begin(), GFst.end(), str2gobs(chosen_C1)) != GFst.end())
                    BDS_C2[i + 1][j] = obs.getobs(str2gobs(chosen_C1));
  
                if (!chosen_C2.empty() && std::find(GFst.begin(), GFst.end(), str2gobs(chosen_C2)) != GFst.end())
                    (chosen_C2 == "C6I" ? BDS_C6[i + 1][j] : BDS_C7[i + 1][j]) = obs.getobs(str2gobs(chosen_C2));

                if (!chosen_L1.empty() && std::find(GFst.begin(), GFst.end(), str2gobs(chosen_L1)) != GFst.end())
                    BDS_L2[i + 1][j] = obs.getobs(str2gobs(chosen_L1));
    
                if (!chosen_L2.empty() && std::find(GFst.begin(), GFst.end(), str2gobs(chosen_L2)) != GFst.end())
                    (chosen_L2 == "L6I" ? BDS_L6[i + 1][j] : BDS_L7[i + 1][j]) = obs.getobs(str2gobs(chosen_L2));
            }
            epochs.push_back(epo);
        }

        for (int i = 1; i <= 46; ++i) {
            for (int j = 1; j <= 2880; ++j) {
                OBS.C1[i][j] = BDS_C2[j][i];
                OBS.L1[i][j] = BDS_L2[j][i];
                if (!chosen_C2.empty()) {
                    if (chosen_C2 == "C6I") OBS.C2[i][j] = BDS_C6[j][i];
                    else if (chosen_C2 == "C7I") OBS.C2[i][j] = BDS_C7[j][i];
                }
                if (!chosen_L2.empty()) {
                    if (chosen_L2 == "L6I") OBS.L2[i][j] = BDS_L6[j][i];
                    else if (chosen_L2 == "L7I") OBS.L2[i][j] = BDS_L7[j][i];
                }
            }
        }
    }
    void extract_GLO_obs(
        t_gallobs* gobs, const std::string& station,
        std::vector<std::vector<double>>& GLO_C1,
        std::vector<std::vector<double>>& GLO_C2,
        std::vector<std::vector<double>>& GLO_L1,
        std::vector<std::vector<double>>& GLO_L2,
        std::vector<t_gtime>& epochs,
        std::vector<std::string>& GLO_sats, obs& OBS,
        std::shared_ptr<spdlog::logger> logger)
    {
        std::vector<std::string> GLO_C1_type = { "C1P", "P1" };
        std::vector<std::string> GLO_C2_type = { "C2P", "P2" };
        std::vector<std::string> GLO_L1_type = { "L1P", "L1" };
        std::vector<std::string> GLO_L2_type = { "L2P", "L2" };

        const size_t num_epochs = 2881;
        const size_t num_sats = 25;

        GLO_C1.assign(num_epochs, std::vector<double>(num_sats, 0.0));
        GLO_C2.assign(num_epochs, std::vector<double>(num_sats, 0.0));
        GLO_L1.assign(num_epochs, std::vector<double>(num_sats, 0.0));
        GLO_L2.assign(num_epochs, std::vector<double>(num_sats, 0.0));

        std::vector<t_gtime> all_epochs = gobs->epochs(station);

        GLO_sats.clear();
        for (int r = 1; r <= 24; ++r) {
            std::ostringstream oss;
            oss << "R" << std::setw(2) << std::setfill('0') << r;
            GLO_sats.push_back(oss.str());
        }

        std::string chosen_C1 = "", chosen_C2 = "", chosen_L1 = "", chosen_L2 = "";

        for (const auto& type : GLO_C1_type) {
            bool found = false;
            for (size_t i = 0; i < all_epochs.size(); ++i) {
                auto obsvec = gobs->obs(station, all_epochs[i]);
                for (const auto& obs : obsvec) {
                    if (obs.sat()[0] != 'R') continue;
                    auto GFst = obs.obs();
                    if (std::find(GFst.begin(), GFst.end(), str2gobs(type)) != GFst.end()) {
                        chosen_C1 = type;
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
            if (!chosen_C1.empty()) break;
        }
        for (const auto& type : GLO_C2_type) {
            bool found = false;
            for (size_t i = 0; i < all_epochs.size(); ++i) {
                auto obsvec = gobs->obs(station, all_epochs[i]);
                for (const auto& obs : obsvec) {
                    if (obs.sat()[0] != 'R') continue;
                    auto GFst = obs.obs();
                    if (std::find(GFst.begin(), GFst.end(), str2gobs(type)) != GFst.end()) {
                        chosen_C2 = type;
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
            if (!chosen_C2.empty()) break;
        }

        for (const auto& type : GLO_L1_type) {
            bool found = false;
            for (size_t i = 0; i < all_epochs.size(); ++i) {
                auto obsvec = gobs->obs(station, all_epochs[i]);
                for (const auto& obs : obsvec) {
                    if (obs.sat()[0] != 'R') continue;
                    auto GFst = obs.obs();
                    if (std::find(GFst.begin(), GFst.end(), str2gobs(type)) != GFst.end()) {
                        chosen_L1 = type;
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
            if (!chosen_L1.empty()) break;
        }

        for (const auto& type : GLO_L2_type) {
            bool found = false;
            for (size_t i = 0; i < all_epochs.size(); ++i) {
                auto obsvec = gobs->obs(station, all_epochs[i]);
                for (const auto& obs : obsvec) {
                    if (obs.sat()[0] != 'R') continue;
                    auto GFst = obs.obs();
                    if (std::find(GFst.begin(), GFst.end(), str2gobs(type)) != GFst.end()) {
                        chosen_L2 = type;
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
            if (!chosen_L2.empty()) break;
        }

        logger->info("[{}] GLO Selected types: C1={}, C2={}, L1={}, L2={}", station, chosen_C1, chosen_C2, chosen_L1, chosen_L2);

        if (chosen_C1.empty()) logger->warn("[{}] No C1 observation type found!", station);
        if (chosen_C2.empty()) logger->warn("[{}] No C2 observation type found!", station);
        if (chosen_L1.empty()) logger->warn("[{}] No L1 observation type found!", station);
        if (chosen_L2.empty()) logger->warn("[{}] No L2 observation type found!", station);

        for (size_t i = 0; i < std::min<size_t>(all_epochs.size(), num_epochs - 1); ++i) {
            const t_gtime& epo = all_epochs[i];
            std::vector<t_gsatdata> obsvec = gobs->obs(station, epo);

            for (const auto& obs : obsvec) {
                std::string prn = obs.sat();
                if (prn[0] != 'R') continue;

                auto it = std::find(GLO_sats.begin(), GLO_sats.end(), prn);
                if (it == GLO_sats.end()) continue;
                size_t j = std::distance(GLO_sats.begin(), it) + 1;
                if (j >= num_sats) continue;

                auto GFst = obs.obs();

                if (!chosen_C1.empty() && std::find(GFst.begin(), GFst.end(), str2gobs(chosen_C1)) != GFst.end())
                    GLO_C1[i + 1][j] = obs.getobs(str2gobs(chosen_C1));
                if (!chosen_C2.empty() && std::find(GFst.begin(), GFst.end(), str2gobs(chosen_C2)) != GFst.end())
                    GLO_C2[i + 1][j] = obs.getobs(str2gobs(chosen_C2));
                if (!chosen_L1.empty() && std::find(GFst.begin(), GFst.end(), str2gobs(chosen_L1)) != GFst.end())
                    GLO_L1[i + 1][j] = obs.getobs(str2gobs(chosen_L1));
                if (!chosen_L2.empty() && std::find(GFst.begin(), GFst.end(), str2gobs(chosen_L2)) != GFst.end())
                    GLO_L2[i + 1][j] = obs.getobs(str2gobs(chosen_L2));
            }
            epochs.push_back(epo);
        }

        for (int i = 1; i <= 24; i++) {
            for (int j = 1; j <= 2880; j++) {
                OBS.C1[i][j] = GLO_C1[j][i];
                OBS.C2[i][j] = GLO_C2[j][i];
                OBS.L1[i][j] = GLO_L1[j][i];
                OBS.L2[i][j] = GLO_L2[j][i];
            }
        }
    }
    void extract_GAL_obs(t_gallobs* gobs, const std::string& station,
        std::vector<std::vector<double>>& GAL_C1,
        std::vector<std::vector<double>>& GAL_C5,
        std::vector<std::vector<double>>& GAL_L1,
        std::vector<std::vector<double>>& GAL_L5,
        std::vector<t_gtime>& epochs,
        std::vector<std::string>& GAL_sats,
        obs& OBS,
        std::shared_ptr<spdlog::logger> logger,
        bool& isC1QAllZero 
    ) {
        const size_t num_epochs = 2881;
        const size_t num_sats = 37;

        GAL_C1.assign(num_epochs, std::vector<double>(num_sats, 0.0));
        GAL_C5.assign(num_epochs, std::vector<double>(num_sats, 0.0));
        GAL_L1.assign(num_epochs, std::vector<double>(num_sats, 0.0));
        GAL_L5.assign(num_epochs, std::vector<double>(num_sats, 0.0));

        std::vector<t_gtime> all_epochs = gobs->epochs(station);

        GAL_sats.clear();
        for (int e = 1; e <= 36; ++e) {
            std::ostringstream oss;
            oss << "E" << std::setw(2) << std::setfill('0') << e;
            GAL_sats.push_back(oss.str());
        }

        std::vector<std::string> C1_priority = { "C1C", "C1X", "C1" };
        std::vector<std::string> C5_priority = { "C5Q", "C5X", "C5" };
        std::vector<std::string> L1_priority = { "L1C", "L1X", "L1" };
        std::vector<std::string> L5_priority = { "L5Q", "L5X", "L5" };

        auto find_first_exist = [&](const std::vector<std::string>& types) -> std::string {
            for (const auto& type : types) {
                for (size_t i = 0; i < all_epochs.size(); ++i) {
                    std::vector<t_gsatdata> obsvec = gobs->obs(station, all_epochs[i]);
                    for (const auto& obs : obsvec) {
                        if (obs.sat()[0] != 'E') continue;
                        const auto& GFst = obs.obs();
                        if (std::find(GFst.begin(), GFst.end(), str2gobs(type)) != GFst.end())
                            return type;
                    }
                }
            }
            return "";
            };

        std::string chosen_C1 = find_first_exist(C1_priority);
        std::string chosen_C2 = find_first_exist(C5_priority);
        std::string chosen_L1 = find_first_exist(L1_priority);
        std::string chosen_L2 = find_first_exist(L5_priority);

        isC1QAllZero = (chosen_C1 != "C1C" || chosen_C2 != "C5Q" || chosen_L1 != "L1C" || chosen_L2 != "L5Q");

        for (size_t i = 0; i < std::min<size_t>(all_epochs.size(), num_epochs - 1); ++i) {
            const t_gtime& epo = all_epochs[i];
            std::vector<t_gsatdata> obsvec = gobs->obs(station, epo);
            for (const auto& obs : obsvec) {
                std::string prn = obs.sat();
                if (prn[0] != 'E') continue;
                auto it = std::find(GAL_sats.begin(), GAL_sats.end(), prn);
                if (it == GAL_sats.end()) continue;
                size_t j = std::distance(GAL_sats.begin(), it) + 1;
                if (j >= num_sats) continue;

                auto has = [](const t_gsatdata& o, const GOBS& g) {
                    const auto& GFst = o.obs();
                    return std::find(GFst.begin(), GFst.end(), g) != GFst.end();
                    };
                GOBS g;
                if (!chosen_C1.empty()) { g = str2gobs(chosen_C1); if (has(obs, g)) GAL_C1[i + 1][j] = obs.getobs(g); }
                if (!chosen_C2.empty()) { g = str2gobs(chosen_C2); if (has(obs, g)) GAL_C5[i + 1][j] = obs.getobs(g); }
                if (!chosen_L1.empty()) { g = str2gobs(chosen_L1); if (has(obs, g)) GAL_L1[i + 1][j] = obs.getobs(g); }
                if (!chosen_L2.empty()) { g = str2gobs(chosen_L2); if (has(obs, g)) GAL_L5[i + 1][j] = obs.getobs(g); }
            }
            epochs.push_back(epo);
        }

        logger->info("[{}] GAL Selected types: C1={}, C2={}, L1={}, L2={}", station, chosen_C1, chosen_C2, chosen_L1, chosen_L2);

        if (chosen_C1.empty()) logger->warn("[{}] No C1 observation type found!", station);
        if (chosen_C2.empty()) logger->warn("[{}] No C2 observation type found!", station);
        if (chosen_L1.empty()) logger->warn("[{}] No L1 observation type found!", station);
        if (chosen_L2.empty()) logger->warn("[{}] No L2 observation type found!", station);

        for (int i = 1; i <= 36; i++) {
            for (int j = 1; j <= 2880; j++) {
                OBS.C1[i][j] = GAL_C1[j][i];
                OBS.C2[i][j] = GAL_C5[j][i];
                OBS.L1[i][j] = GAL_L1[j][i];
                OBS.L2[i][j] = GAL_L5[j][i];
            }
        }
    }





    //// Extract SNR (Signal-to-Noise Ratio) values for GPS S1 and S2, with auto priority selection
    //// Save results to GPS_S1, GPS_S2, and map to OBS struct; also output S1.txt/S2.txt files
    //void extract_GPS_SNR(
    //    t_gallobs* gobs,
    //    const std::string& station,
    //    std::vector<std::vector<double>>& GPS_S1,
    //    std::vector<std::vector<double>>& GPS_S2,
    //    std::vector<t_gtime>& epochs,
    //    std::vector<std::string>& GPS_sats,
    //    std::shared_ptr<spdlog::logger> logger,
    //    obs& OBS)
    //{
    //    // S1 priority: S1C > S1W > S1
    //    std::vector<std::string> S1_priority = { "S1C", "S1W", "S1" };
    //    // S2 priority: S2W > S2X > S2P > S2C > S2
    //    std::vector<std::string> S2_priority = { "S2W", "S2X", "S2P", "S2C", "S2" };

    //    const size_t num_epochs = 2881;
    //    const size_t num_sats = 33;

    //    // Initialize output arrays with zeros
    //    GPS_S1.assign(num_epochs, std::vector<double>(num_sats, 0.0));
    //    GPS_S2.assign(num_epochs, std::vector<double>(num_sats, 0.0));

    //    std::vector<t_gtime> all_epochs = gobs->epochs(station);

    //    // Fill PRN list "G01"..."G32"
    //    GPS_sats.clear();
    //    for (int s = 1; s <= 32; ++s) {
    //        std::ostringstream oss;
    //        oss << "G" << std::setw(2) << std::setfill('0') << s;
    //        GPS_sats.push_back(oss.str());
    //    }

    //    // Auto-select S1/S2 observation code based on availability and priority
    //    std::string chosen_S1 = "";
    //    std::string chosen_S2 = "";

    //    for (size_t i = 0; i < all_epochs.size(); ++i) {
    //        std::vector<t_gsatdata> obsvec = gobs->obs(station, all_epochs[i]);
    //        for (const auto& obs : obsvec) {
    //            if (obs.sat()[0] != 'G') continue;
    //            const auto& GFst = obs.obs();

    //            if (chosen_S1.empty()) {
    //                for (const auto& type : S1_priority) {
    //                    if (std::find(GFst.begin(), GFst.end(), str2gobs(type)) != GFst.end()) {
    //                        chosen_S1 = type;
    //                        break;
    //                    }
    //                }
    //            }
    //            if (chosen_S2.empty()) {
    //                for (const auto& type : S2_priority) {
    //                    if (std::find(GFst.begin(), GFst.end(), str2gobs(type)) != GFst.end()) {
    //                        chosen_S2 = type;
    //                        break;
    //                    }
    //                }
    //            }
    //            if (!chosen_S1.empty() && !chosen_S2.empty()) break;
    //        }
    //        if (!chosen_S1.empty() && !chosen_S2.empty()) break;
    //    }

    //    // Log the selected SNR observation types
    //    logger->info("[{}] Selected SNR types: S1={}, S2={}", station, chosen_S1, chosen_S2);

    //    if (chosen_S1.empty()) logger->warn("[{}] No S1 SNR observation type found!", station);
    //    if (chosen_S2.empty()) logger->warn("[{}] No S2 SNR observation type found!", station);

    //    // Fill the 2D SNR arrays for each epoch and satellite
    //    for (size_t i = 0; i < std::min<size_t>(all_epochs.size(), num_epochs - 1); ++i) {
    //        const t_gtime& epo = all_epochs[i];
    //        std::vector<t_gsatdata> obsvec = gobs->obs(station, epo);

    //        for (const auto& obs : obsvec) {
    //            std::string prn = obs.sat();
    //            if (prn[0] != 'G') continue;

    //            auto it = std::find(GPS_sats.begin(), GPS_sats.end(), prn);
    //            if (it == GPS_sats.end()) continue;
    //            size_t j = std::distance(GPS_sats.begin(), it) + 1;
    //            if (j >= num_sats) continue;

    //            const auto& GFst = obs.obs();

    //            if (!chosen_S1.empty() && std::find(GFst.begin(), GFst.end(), str2gobs(chosen_S1)) != GFst.end())
    //                GPS_S1[i + 1][j] = obs.getobs(str2gobs(chosen_S1));

    //            if (!chosen_S2.empty() && std::find(GFst.begin(), GFst.end(), str2gobs(chosen_S2)) != GFst.end())
    //                GPS_S2[i + 1][j] = obs.getobs(str2gobs(chosen_S2));
    //        }
    //        epochs.push_back(epo);
    //    }

    //    // Optionally: map to OBS struct (indices: [sat][epoch])
    //    for (int i = 1; i <= 32; i++) {
    //        for (int j = 1; j <= 2880; j++) {
    //            OBS.S1[i][j] = GPS_S1[j][i];
    //            OBS.S2[i][j] = GPS_S2[j][i];
    //        }
    //    }

    //    // Output S1.txt file: rows for epochs, columns for PRNs
    //    std::ofstream fout_s1("S1.txt");
    //    for (int j = 1; j <= 2880; j++) {
    //        for (int i = 1; i <= 32; i++) {
    //            fout_s1 << OBS.S1[i][j];
    //            if (i < 32) fout_s1 << "\t";
    //        }
    //        fout_s1 << std::endl;
    //    }
    //    fout_s1.close();

    //    // Output S2.txt file: rows for epochs, columns for PRNs
    //    std::ofstream fout_s2("S2.txt");
    //    for (int j = 1; j <= 2880; j++) {
    //        for (int i = 1; i <= 32; i++) {
    //            fout_s2 << OBS.S2[i][j];
    //            if (i < 32) fout_s2 << "\t";
    //        }
    //        fout_s2 << std::endl;
    //    }
    //    fout_s2.close();

    //}    

    
    
}