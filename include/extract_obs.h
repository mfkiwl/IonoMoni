#pragma once

#include <vector>
#include <string>
#include "gall/gallobs.h"
#include "gset/gsetbase.h"
#include "gset/gsetgen.h"
#include "gutils/gobs.h"
#include "gutils/gtime.h"
#include <memory>          
#include <spdlog/logger.h> 
#include "obs.h"
namespace gnut {

    // 写矩阵到txt（模板）
    template<typename T>
    void write_matrix_to_txt_raw(const T& filename, const std::vector<std::vector<double>>& matrix)
    {
        std::ofstream fout(filename);
        if (!fout.is_open())
        {
            std::cerr << "Cannot open output file: " << filename << std::endl;
            return;
        }

        for (size_t i = 1; i < matrix.size(); ++i)
        {
            for (size_t j = 1; j < matrix[i].size(); ++j)
            {
                fout << std::fixed << std::setprecision(3) << matrix[i][j];
                if (j < matrix[i].size() - 1)
                    fout << "\t";
            }
            fout << "\n";
        }

        fout.close();
    }

    // GPS


    void extract_GPS_obs(t_gallobs* gobs, const std::string& station,
        std::vector<std::vector<double>>& GPS_C1,
        std::vector<std::vector<double>>& GPS_C2,
        std::vector<std::vector<double>>& GPS_L1,
        std::vector<std::vector<double>>& GPS_L2,
        std::vector<t_gtime>& epochs,
        std::vector<std::string>& GPS_sats,
        bool& isC1WAllZero,
        std::shared_ptr<spdlog::logger> logger,obs &OBS);  

    /*void extract_GPS_SNR(t_gallobs* gobs, const std::string& station,
        std::vector<std::vector<double>>& GPS_S1,
        std::vector<std::vector<double>>& GPS_S2,
        std::vector<t_gtime>& epochs,
        std::vector<std::string>& GPS_sats,
        std::shared_ptr<spdlog::logger> logger, 
        obs& OBS);*/



    // BDS
    void extract_BDS_obs(t_gallobs* gobs, const std::string& station,
        std::vector<std::vector<double>>& BDS_C2,
        std::vector<std::vector<double>>& BDS_C6,
        std::vector<std::vector<double>>& BDS_C7,
        std::vector<std::vector<double>>& BDS_L2,
        std::vector<std::vector<double>>& BDS_MW,
        std::vector<std::vector<double>>& BDS_L7,
        std::vector<t_gtime>& epochs,
        std::vector<std::string>& BDS_sats,
        obs& OBS,
        std::shared_ptr<spdlog::logger> logger,
        bool& isC7IAllZero);  

    // GLO
    void extract_GLO_obs(t_gallobs* gobs, const std::string& station,
        std::vector<std::vector<double>>& GLO_C1,
        std::vector<std::vector<double>>& GLO_C2,
        std::vector<std::vector<double>>& GLO_L1,
        std::vector<std::vector<double>>& GLO_L2,
        std::vector<t_gtime>& epochs,
        std::vector<std::string>& GLO_sats,
        obs& OBS,
        std::shared_ptr<spdlog::logger> logger);

    void extract_GAL_obs(
        t_gallobs* gobs,
        const std::string& station,
        std::vector<std::vector<double>>& GAL_C1,
        std::vector<std::vector<double>>& GAL_C5,
        std::vector<std::vector<double>>& GAL_L1,
        std::vector<std::vector<double>>& GAL_L5,
        std::vector<t_gtime>& epochs,
        std::vector<std::string>& GAL_sats,
        obs& OBS,
        std::shared_ptr<spdlog::logger> logger,
        bool& isC1QAllZero    
    );


} // namespace gnut
