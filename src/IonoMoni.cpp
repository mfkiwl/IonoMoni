#define _HAS_STD_BYTE 0    
#define NOMINMAX 
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <memory>
#include <signal.h>
#include <windows.h>
#include "extract_obs.h"
#include "src/libGREAT/gcfg_ppp.h"
#include "gcoders/rinexo.h"
#include "gio/gfile.h"
#include "gall/gallobs.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "obs.h"
#include "decode_obs.h"
#include "process_sp3.h"
#include "read_sp3.h"
#include "LagrangeInterpolation.h"
#include "read_Time_sp3.h"
#include "sp3.h"
#include "get_sp3_paths.h"
#include "stringToDouble.h"
#include "calc_elevation.h"
#include "extract_TEC.h"
#include "get_dcb_path.h"
#include <sys/stat.h>
#include <direct.h>
#include <io.h>   
#include "gall/gallproc.h"
#include "gdata/gobj.h"
#include "gproc/gpvtflt.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "calc_AATR.h"
#include "calc_ROTI.h"
#include "constants.h"
#include "ecef2elli.h"
#include "calc_ipp.h"
#include "mapping_function.h"
#include "read_DCB_1.0.h"
#include "ppp+.h"
#include <cmath>

double h = 350000.0;  // Default value to avoid linker errors
namespace fs = std::filesystem;
using namespace std;
using namespace gnut;


// Main processing entry points for each function
void run_ccl_stec(t_gcfg_ppp& gset, shared_ptr<spdlog::logger> my_logger, const string& output_dir);
void run_ppp_stec(t_gcfg_ppp& gset, shared_ptr<spdlog::logger> my_logger);
void run_roti(t_gcfg_ppp& gset, std::shared_ptr<spdlog::logger> my_logger);
void run_aatr(t_gcfg_ppp& gset, std::shared_ptr<spdlog::logger> my_logger);


int main(int argc, char* argv[])
{
    std::string output_dir = "result";
    // Support the -x parameter; if not provided, use the default XML file.
    std::string xmlfile = "xml/IonoMoni.xml";
    bool xml_provided = false;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if ((arg == "-x" || arg == "--xml") && i + 1 < argc) {
            xmlfile = argv[i + 1];
            xml_provided = true;
            break;
        }
    }

    if (!xml_provided) {
        std::cout << "[Notice] No XML specified via -x; using default path: " << xmlfile << std::endl;
    }
    else {
        std::cout << "[Notice] Using provided XML configuration file: " << xmlfile << std::endl;
    }

    struct _stat info;
    if (_stat(output_dir.c_str(), &info) != 0 || !(info.st_mode & _S_IFDIR))
    {
        _mkdir(output_dir.c_str());
    }

    // Check if the XML file exists
    if (_access(xmlfile.c_str(), 0) != 0) {
        std::cerr << "[Error] XML file does not exist: " << xmlfile << std::endl;
        return -1;
    }

    std::cout << "Attempting to read XML file: " << xmlfile << std::endl;

    // Read XML
    t_gcfg_ppp gset;
    try {
        if (gset.read(xmlfile) != 0) {
            std::cerr << "Failed to read XML!" << std::endl;
            return -2;
        }
        gset.check();
    }
    catch (const std::exception& e) {
        std::cerr << "Exception occurred while reading XML: " << e.what() << std::endl;
        return -3;
    }


    h = gset.ion_height();

    // Extract the function field
    std::string function_mode = gset.comment("function");
    if (function_mode.empty()) {
        std::cerr << "Function field not found or empty in XML!" << std::endl;
        return -4;
    }

    // Create logger
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    std::string log_filename = output_dir + "/" + function_mode + ".log";
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_filename, true);
    std::vector<spdlog::sink_ptr> sinks{ console_sink, file_sink };
    auto my_logger = std::make_shared<spdlog::logger>("multi_logger", sinks.begin(), sinks.end());
    spdlog::register_logger(my_logger);
    my_logger->set_level(spdlog::level::info);

    // === Logging starts ===
    my_logger->info("IonoMoni program started...");
    my_logger->info("Successfully read XML: {}", xmlfile);
    my_logger->info("Successfully retrieved function field: {}", function_mode);

    // Run the corresponding function
    if (function_mode == "CCL_STEC")
        run_ccl_stec(gset, my_logger, output_dir);
    else if (function_mode == "PPP_STEC")
        run_ppp_stec(gset, my_logger);
    else if (function_mode == "ROTI")
        run_roti(gset, my_logger);
    else if (function_mode == "AATR")
        run_aatr(gset, my_logger);
    else {
        my_logger->error("Unrecognized function type: {}", function_mode);
        return -5;
    }

    my_logger->info("IonoMoni program finished.");
    return 0;
}

void apply_elevation_mask(obs& OBS, const sp3& SP3_data, double elevCutoffDeg, int maxSat, char sys) {
    for (int j = 1; j <= 2880; j++) {
        for (int i = 1; i <= maxSat; i++) {
            if (OBS.C1[i][j] == 0 || OBS.C2[i][j] == 0 || OBS.L1[i][j] == 0 || OBS.L2[i][j] == 0) {
                OBS.C1[i][j] = OBS.C2[i][j] = OBS.L1[i][j] = OBS.L2[i][j] = 0;
            }
            double E, A;
            switch (sys) {
            case 'G': calc_elevation(SP3_data.X[i][j] * 1000, SP3_data.Y[i][j] * 1000, SP3_data.Z[i][j] * 1000, OBS.X, OBS.Y, OBS.Z, E, A); break;
            case 'C': calc_elevation(SP3_data.CX[i][j] * 1000, SP3_data.CY[i][j] * 1000, SP3_data.CZ[i][j] * 1000, OBS.X, OBS.Y, OBS.Z, E, A); break;
            case 'R': calc_elevation(SP3_data.RX[i][j] * 1000, SP3_data.RY[i][j] * 1000, SP3_data.RZ[i][j] * 1000, OBS.X, OBS.Y, OBS.Z, E, A); break;
            case 'E': calc_elevation(SP3_data.EX[i][j] * 1000, SP3_data.EY[i][j] * 1000, SP3_data.EZ[i][j] * 1000, OBS.X, OBS.Y, OBS.Z, E, A); break;
            }
            if (E < elevCutoffDeg) {
                OBS.C1[i][j] = OBS.C2[i][j] = OBS.L1[i][j] = OBS.L2[i][j] = 0;
               //OBS.S1[i][j] = OBS.S2[i][j] = 0;
            }
        }
    }
}

// CCL_STEC function

void run_ccl_stec(t_gcfg_ppp& gset, shared_ptr<spdlog::logger> my_logger, const string& output_dir)
{
    my_logger->info("Starting SMOOTH_STEC function...");
    std::string raw_txt_path = gset.outputs("txt");
    if (raw_txt_path.find("file://") == 0) raw_txt_path = raw_txt_path.substr(7);
    if (raw_txt_path.empty()) raw_txt_path = "result/CCL_STEC_result/$(rec)-$(function).txt";

    std::string function_name = gset.outputs("function");
    if (function_name.empty()) function_name = "CCL_STEC";

    auto replace_all = [](std::string& str, const std::string& from, const std::string& to) {
        size_t pos = 0;
        while ((pos = str.find(from, pos)) != std::string::npos) {
            str.replace(pos, from.length(), to);
            pos += to.length();
        }
        };
    int mappingFuncType = gset.mapping_function();
    double elevCutoffDeg = gset.minimum_elev();
    std::set<std::string> sys_set = gset.sys();
    std::string sys_str;
    for (const auto& s : sys_set) {
        if (!sys_str.empty()) sys_str += " ";
        sys_str += s;
    }

    t_gallobs* gobs = decode_obs(gset, my_logger);
    set<string> stations = dynamic_cast<t_gsetgen*>(&gset)->recs();
 
    for (const auto& station : stations) {
        auto runepoch = t_gtime::current_time(t_gtime::GPS);
        obs OBS;
        std::string stationName = station;

        my_logger->info("Start processing station: {}", station);
        std::string txt_output_path = raw_txt_path;
        replace_all(txt_output_path, "$(rec)", station);
        replace_all(txt_output_path, "$(function)", function_name);
        std::filesystem::create_directories(std::filesystem::path(txt_output_path).parent_path());

        std::vector<t_gtime> epochs;
        std::vector<double> xyz = gobs->get_station_coordinates(station);
        OBS.X = xyz[0], OBS.Y = xyz[1], OBS.Z = xyz[2];

        std::string dcbFilePath = getDcbFilePath(gset);
        std::vector<std::string> filenamesSp3 = getSortedSp3Paths(gset);

        sp3 SP3[3];
        for (int i = 0; i < 3; i++) {
            std::ifstream file(filenamesSp3[i]);
            if (file.is_open()) {
                read_Time_sp3(file, filenamesSp3[i], SP3[i]);
                getSp3Data(file, filenamesSp3[i], SP3[i]);
            }
        }
        processSP3(SP3);

        // Use a separate try-catch block for each system.
        if (sys_str.find("GPS") != std::string::npos) {
            try {
                OBS.reset();
                std::vector<std::vector<double>> GPS_C1, GPS_C2, GPS_L1, GPS_L2;
                std::vector<std::string> GPS_sats;
                bool isC1WAllZero = true;
                extract_GPS_obs(gobs, station, GPS_C1, GPS_C2, GPS_L1, GPS_L2, epochs, GPS_sats, isC1WAllZero, my_logger, OBS);
                apply_elevation_mask(OBS, SP3[1], elevCutoffDeg, 32, 'G');
                G_prepro(OBS, stationName, isC1WAllZero, dcbFilePath, txt_output_path, gset, SP3, mappingFuncType, my_logger);
                std::vector<std::vector<IPP_Point>> IPP;
                calc_ipp_matrix(OBS.X, OBS.Y, OBS.Z,SP3[1].X, SP3[1].Y, SP3[1].Z,32, 2880,"GPS",txt_output_path);
            }
            catch (const std::exception& ex) {
                my_logger->error("Station {} [GPS] skipped: {}", station, ex.what());
            }
        }

        if (sys_str.find("BDS") != std::string::npos) {
            try {
                OBS.reset();
                std::vector<std::vector<double>> BDS_C2, BDS_C6, BDS_C7, BDS_L2, BDS_L6, BDS_L7;
                std::vector<std::string> BDS_sats;
                bool isC7IAllZero = true;
                extract_BDS_obs(gobs, station, BDS_C2, BDS_C6, BDS_C7, BDS_L2, BDS_L6, BDS_L7, epochs, BDS_sats, OBS, my_logger, isC7IAllZero);
                apply_elevation_mask(OBS, SP3[1], elevCutoffDeg, 46, 'C');
                C_prepro(OBS, stationName, isC7IAllZero, dcbFilePath, txt_output_path, gset, SP3, mappingFuncType, my_logger);
                std::vector<std::vector<IPP_Point>> IPP;
                calc_ipp_matrix(OBS.X, OBS.Y, OBS.Z,SP3[1].CX, SP3[1].CY, SP3[1].CZ,46, 2880,"BDS",txt_output_path);
            }
            catch (const std::exception& ex) {
                my_logger->error("Station {} [BDS] skipped: {}", station, ex.what());
            }
        }

        // GLO
        if (sys_str.find("GLO") != std::string::npos) {
            try {
                OBS.reset();
                std::vector<std::vector<double>> GLO_C1, GLO_C2, GLO_L1, GLO_L2;
                std::vector<std::string> GLO_sats;
                extract_GLO_obs(gobs, station, GLO_C1, GLO_C2, GLO_L1, GLO_L2, epochs, GLO_sats, OBS, my_logger);
                apply_elevation_mask(OBS, SP3[1], elevCutoffDeg, 24, 'R');
                R_prepro(OBS, stationName, dcbFilePath, txt_output_path, gset, SP3, mappingFuncType, my_logger);
                std::vector<std::vector<IPP_Point>> IPP;
                calc_ipp_matrix(OBS.X, OBS.Y, OBS.Z,SP3[1].RX, SP3[1].RY, SP3[1].RZ,24, 2880,"GLO",txt_output_path);
            }
            catch (const std::exception& ex) {
                my_logger->error("Station {} [GLO] skipped: {}", station, ex.what());
            }
        }

        // GAL
        if (sys_str.find("GAL") != std::string::npos) {
            try {
                OBS.reset();
                std::vector<std::vector<double>> GAL_C1, GAL_C5, GAL_L1, GAL_L5;
                std::vector<std::string> GAL_sats;
                bool isC1QAllZero = false; 
                extract_GAL_obs(gobs, station, GAL_C1, GAL_C5, GAL_L1, GAL_L5, epochs, GAL_sats, OBS, my_logger, isC1QAllZero);
                apply_elevation_mask(OBS, SP3[1], elevCutoffDeg, 36, 'E');
                E_prepro(OBS, stationName, dcbFilePath, txt_output_path, gset, SP3, mappingFuncType, my_logger);
                std::vector<std::vector<IPP_Point>> IPP;
                calc_ipp_matrix(OBS.X, OBS.Y, OBS.Z,SP3[1].EX, SP3[1].EY, SP3[1].EZ,36, 2880,"GAL",txt_output_path);
            }
            catch (const std::exception& ex) {
                my_logger->error("Station {} [GAL] skipped: {}", station, ex.what());
            }
        }

        auto lstepoch = t_gtime::current_time(t_gtime::GPS);
        my_logger->info("Station {} output completed", station);
        my_logger->info("Station {} CCL processing completed in {:.3f} seconds", station, lstepoch.diff(runepoch));
    }

    delete gobs;
}

// PPP_STEC Function
const int GLO_FCN[24] = { 1, -4, 5, 6, 1, -4, 5, 6, -2, -7, 0, -1,
               -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2 };
void run_ppp_stec(t_gcfg_ppp& gset, shared_ptr<spdlog::logger> my_logger)
{
    std::vector<std::vector<double>> GPS_PPP_STEC(2881, std::vector<double>(33, 0.0));
    std::vector<std::vector<double>> GLO_PPP_STEC(2881, std::vector<double>(25, 0.0));
    std::vector<std::vector<double>> GAL_PPP_STEC(2881, std::vector<double>(37, 0.0));
    std::vector<std::vector<double>> BDS_PPP_STEC(2881, std::vector<double>(47, 0.0));

    my_logger->info("Starting PPP_STEC function...");
    
    // Initialize input objects
    t_gdata* gdata = nullptr;
    t_gnavde* gde = new t_gnavde;
    t_gpoleut1* gerp = new t_gpoleut1;
    t_gallobs* gobs = new t_gallobs();
    gobs->spdlog(my_logger);
    gobs->gset(&gset);

    t_gallprec* gorb = new t_gallprec();
    gorb->spdlog(my_logger);

    t_gallpcv* gpcv = nullptr;
    if (gset.input_size("atx") > 0) gpcv = new t_gallpcv;
    if (gpcv) gpcv->spdlog(my_logger);

    t_gallotl* gotl = nullptr;
    if (gset.input_size("blq") > 0) gotl = new t_gallotl;
    if (gotl) gotl->spdlog(my_logger);

    t_gallbias* gbia = nullptr;
    if (gset.input_size("biasinex") > 0 || gset.input_size("bias") > 0) gbia = new t_gallbias;
    if (gbia) gbia->spdlog(my_logger);

    t_gallobj* gobj = new t_gallobj(my_logger, gpcv, gotl);
    gobj->spdlog(my_logger);

    t_gupd* gupd = nullptr;
    if (gset.input_size("upd") > 0) gupd = new t_gupd;
    if (gupd) gupd->spdlog(my_logger);

    t_gifcb* gifcb = nullptr;
    if (gset.input_size("ifcb") > 0) gifcb = new t_gifcb;
    if (gifcb) gifcb->spdlog(my_logger);

    multimap<IFMT, string> inp = gset.inputs_all();
    for (auto& itINP : inp)
    {
        IFMT ifmt(itINP.first);
        string path(itINP.second);
        t_gio* tgio = nullptr;
        t_gcoder* tgcoder = nullptr;

        if (ifmt == IFMT::RINEXO_INP) { gdata = gobs; tgcoder = new t_rinexo(&gset, "", 4096); }
        else if (ifmt == IFMT::SP3_INP) { gdata = gorb; tgcoder = new t_sp3(&gset, "", 8192); }
        else if (ifmt == IFMT::RINEXC_INP) { gdata = gorb; tgcoder = new t_rinexc(&gset, "", 4096); }
        else if (ifmt == IFMT::RINEXN_INP) { gdata = gorb; tgcoder = new t_rinexn(&gset, "", 4096); }
        else if (ifmt == IFMT::ATX_INP) { gdata = gpcv; tgcoder = new t_atx(&gset, "", 4096); }
        else if (ifmt == IFMT::BLQ_INP) { gdata = gotl; tgcoder = new t_blq(&gset, "", 4096); }
        else if (ifmt == IFMT::UPD_INP) { gdata = gupd; tgcoder = new t_upd(&gset, "", 4096); }
        else if (ifmt == IFMT::BIASINEX_INP) { gdata = gbia; tgcoder = new t_biasinex(&gset, "", 20480); }
        else if (ifmt == IFMT::BIAS_INP) { gdata = gbia; tgcoder = new t_biabernese(&gset, "", 20480); }
        else if (ifmt == IFMT::DE_INP) { gdata = gde; tgcoder = new t_dvpteph405(&gset, "", 4096); }
        else if (ifmt == IFMT::EOP_INP) { gdata = gerp; tgcoder = new t_poleut1(&gset, "", 4096); }
        else if (ifmt == IFMT::IFCB_INP) { gdata = gifcb; tgcoder = new t_ifcb(&gset, "", 4096); }
        else {
            my_logger->warn("Unrecognized input format, skipped: {}", int(ifmt));
            continue;
        }

        if (path.substr(0, 7) == "file://") {
            tgio = new t_gfile(my_logger);
            tgio->spdlog(my_logger);
            tgio->path(path);
        }

        if (tgcoder && tgio) {
            tgcoder->clear();
            tgcoder->path(path);
            tgcoder->spdlog(my_logger);
            tgcoder->add_data("ID" + std::to_string(&itINP - &(*inp.begin())), gdata);
            tgcoder->add_data("OBJ", gobj);
            tgio->coder(tgcoder);
            tgio->run_read();
            delete tgio;
            delete tgcoder;
        }
    }

    t_gtime beg = dynamic_cast<t_gsetgen*>(&gset)->beg();
    gobj->read_satinfo(beg);
    gobj->sync_pcvs();

    t_gallproc* data = new t_gallproc();
    data->Add_Data(t_gdata::type2str(gobs->id_type()), gobs);
    data->Add_Data(t_gdata::type2str(gorb->id_type()), gorb);
    data->Add_Data(t_gdata::type2str(gobj->id_type()), gobj);
    if (gbia) data->Add_Data(t_gdata::type2str(gbia->id_type()), gbia);
    if (gotl) data->Add_Data(t_gdata::type2str(gotl->id_type()), gotl);
    if (gde)  data->Add_Data(t_gdata::type2str(gde->id_type()), gde);
    if (gerp) data->Add_Data(t_gdata::type2str(gerp->id_type()), gerp);
    if (gupd) data->Add_Data(t_gdata::type2str(gupd->id_type()), gupd);
    if (gifcb) data->Add_Data(t_gdata::type2str(gifcb->id_type()), gifcb);

    std::string sys_str;
    for (const auto& sys : gset.sys()) sys_str += sys + " ";

    int mappingFuncType = gset.mapping_function();

    std::string dcb_file = getDcbFilePath(gset);
    std::map<std::string, double> satelliteDCBMap;
    std::map<std::string, double> receiverDCBMap;

    satelliteDCBMap.clear();
    receiverDCBMap.clear();
    read_dcb_to_map(dcb_file, satelliteDCBMap, receiverDCBMap);

    set<string> stations = dynamic_cast<t_gsetgen*>(&gset)->recs();

    for (const auto& station : stations)
    {
        if (!gobs->isSite(station))
        {
            my_logger->warn("No valid observation data for station {}, skipping...", station);
            continue;
        }

        my_logger->info("Start processing station: {}", station);

        std::fill(GPS_PPP_STEC.begin(), GPS_PPP_STEC.end(), std::vector<double>(33, 0.0));
        std::fill(GLO_PPP_STEC.begin(), GLO_PPP_STEC.end(), std::vector<double>(25, 0.0));
        std::fill(GAL_PPP_STEC.begin(), GAL_PPP_STEC.end(), std::vector<double>(37, 0.0));
        std::fill(BDS_PPP_STEC.begin(), BDS_PPP_STEC.end(), std::vector<double>(47, 0.0));

        obs OBS;
        std::vector<double> xyz = gobs->get_station_coordinates(station);
        OBS.X = xyz[0], OBS.Y = xyz[1], OBS.Z = xyz[2];

        std::vector<std::string> filenamesSp3 = getSortedSp3Paths(gset);

        sp3 SP3[3];
        for (int i = 0; i < 3; i++) {
            std::ifstream file(filenamesSp3[i]);
            if (file.is_open()) getSp3Data(file, filenamesSp3[i], SP3[i]);
        }
        processSP3(SP3);
        std::string raw_txt_path = gset.outputs("txt");

        if (raw_txt_path.find("file://") == 0) {
            raw_txt_path = raw_txt_path.substr(7);
        }

        if (raw_txt_path.empty())
            raw_txt_path = "result/PPP_STEC_result/$(rec)-$(function).txt";  // Default output path

        std::string function_name = gset.outputs("function");  
        if (function_name.empty()) function_name = "PPP_STEC";
        auto replace_all = [](std::string& str, const std::string& from, const std::string& to) {
            size_t pos = 0;
            while ((pos = str.find(from, pos)) != std::string::npos) {
                str.replace(pos, from.length(), to);
                pos += to.length();
            }
            };

        replace_all(raw_txt_path, "$(rec)", station);
        replace_all(raw_txt_path, "$(function)", function_name);

        std::string txt_output_path = raw_txt_path;
        std::filesystem::create_directories(std::filesystem::path(txt_output_path).parent_path());


        t_gpvtflt* pppflt = new t_gpvtflt(station, "", &gset, my_logger, data);

        t_gtime beg = dynamic_cast<t_gsetgen*>(&gset)->beg();
        t_gtime end = dynamic_cast<t_gsetgen*>(&gset)->end();

        auto runepoch = t_gtime::current_time(t_gtime::GPS);

        pppflt->processBatch(beg, end, true);

        auto lstepoch = t_gtime::current_time(t_gtime::GPS);

        my_logger->info("Station {} PPP processing completed in {:.3f} seconds", station, lstepoch.diff(runepoch));

        std::set<std::string> reportedRecKeys, reportedSatKeys;

        {
            for (const auto& entry : pppflt->get_sion_memory()) {
                const t_gtime& epo = std::get<0>(entry);
                const std::string& prn = std::get<1>(entry);
                const par_type& type = std::get<2>(entry);
                double value = std::get<3>(entry);
                const std::string& code_type = std::get<4>(entry);   

                auto map_code = [](const std::string& sys, const std::string& code) -> std::string {
                    if (sys == "G") {
                        if (code == "P1") return "C1W";
                        if (code == "C1") return "C1C";
                        if (code == "P2") return "C2W";
                        if (code == "C2") return "C2C";
                    }
                    else if (sys == "R") {
                        if (code == "P1") return "C1P";
                        if (code == "C1") return "C1C";
                        if (code == "P2") return "C2P";
                        if (code == "C2") return "C2C";
                    }
                    else if (sys == "E") {
                        if (code == "C1") return "C1X";
                        if (code == "C5") return "C5X";
                    }
                    return code;
                };

                std::string sys = prn.substr(0, 1);       
                std::string prn_num = prn.substr(1);        

                // Split code_type into two observation types
                std::string orig_code1, orig_code2;
                size_t comma_pos = code_type.find(',');

                if (comma_pos != std::string::npos) {
                    orig_code1 = code_type.substr(0, comma_pos);
                    orig_code2 = code_type.substr(comma_pos + 1);
                }
                //else {
                //    orig_code1 = code_type;
                //    orig_code2 = "";
                //}
                
                // Remove spaces
                auto trim = [](std::string s) {
                    s.erase(s.find_last_not_of(" \n\r\t") + 1);
                    s.erase(0, s.find_first_not_of(" \n\r\t"));
                    return s;
                    };
                orig_code1 = trim(orig_code1);
                orig_code2 = trim(orig_code2);


                std::string dcb_code1 = map_code(sys, orig_code1);
                std::string dcb_code2 = map_code(sys, orig_code2);
                std::string dcb_code_pair = dcb_code1 + dcb_code2;
                std::string satKey;
                if (sys == "G") {
                    satKey = sys + prn_num + "_C1WC2W";
                }
                else if (sys == "R") {
                    satKey = sys + prn_num + "_C1PC2P";
                }
             
                std::string recKey = sys + "_" + station + "_" + dcb_code_pair;        // eg: G_ABPO_C1WC2W

                bool recDCB_found = (receiverDCBMap.find(recKey) != receiverDCBMap.end());
                bool satDCB_found = (satelliteDCBMap.find(satKey) != satelliteDCBMap.end());

                if (!recDCB_found) {
                    if (reportedRecKeys.find(recKey) == reportedRecKeys.end()) {
                        my_logger->error("Receiver DCB not found for key: {}.", recKey);
                        reportedRecKeys.insert(recKey);
                    }
                    value = 0.0;
                }
                else if (!satDCB_found) {                                          
                    if (reportedSatKeys.find(satKey) == reportedSatKeys.end()) {
                        my_logger->warn("Satellite DCB not found for key: {}.", satKey);
                        reportedSatKeys.insert(satKey);
                    }
                    value = 0.0;
                }
                else {

                    double f1 = 0, f2 = 0;
                    if (sys == "G") {         // GPS
                        f1 = 1575.42e6; f2 = 1227.60e6;
                    }
                    else if (sys == "E") {    // Galileo
                        f1 = 1575.42e6; f2 = 1176.45e6;
                    }
                    else if (sys == "C") {    // BDS
                        f1 = 1561.098e6;  // BDS B1I
                        if (dcb_code_pair.find('6') != std::string::npos) {
                            f2 = 1268.52e6;   // BDS B3I
                        }
                        else if (dcb_code_pair.find('7') != std::string::npos) {
                            f2 = 1207.14e6;   // BDS B2I
                        }
                        else {
                            throw std::runtime_error("In current PPP_STEC mode, only B1I/B2I/B3I observation types are supported for BDS. dcb_code_pair: " + dcb_code_pair);
                        }

                    }

                    else if (sys == "R") {    // GLONASS
                        int sat_idx = std::stoi(prn_num) - 1;
                        f1 = (1602.0e6 + GLO_FCN[sat_idx] * 0.5625e6);
                        f2 = (1246.0e6 + GLO_FCN[sat_idx] * 0.4375e6); 
                    }

                    if (f1 > 0 && f2 > 0) {
                        double satDCB = satelliteDCBMap.at(satKey);
                        double recDCB = receiverDCBMap.at(recKey);
                        double dcb_total_ns = recDCB;  
                        double dcb_total_m = dcb_total_ns * 1e-9 * c;
                        value = value + ((f2 * f2) / (f1 * f1 - f2 * f2)) * dcb_total_m;          
                        //my_logger->info("Epoch {} PRN {} code {}: satDCB={} ns, recDCB={} ns, dcb_total_m={} m, corrected_value={}",
                        //    epo.str_ymdhms(), prn, code_type, satDCB, recDCB, dcb_total_m, value);
                    }

                }

                if (type == par_type::SION) {
                    if (prn.size() >= 3) {
                        char sys = prn[0];
                        int sat_num = stoi(prn.substr(1));
                        int seconds_of_day = static_cast<int>(epo.sod());
                        int epoch_index = seconds_of_day / 30 + 1;

                        if (epoch_index >= 1 && epoch_index <= 2880) {
                            if (sys == 'G' && sat_num > 0 && sat_num < 33) {
                                GPS_PPP_STEC[epoch_index][sat_num] = value;
                            }
                            else if (sys == 'R' && sat_num > 0 && sat_num < 25) {
                                GLO_PPP_STEC[epoch_index][sat_num] = value;
                            }
                            else if (sys == 'E' && sat_num > 0 && sat_num < 37) {
                                GAL_PPP_STEC[epoch_index][sat_num] = value;
                            }
                            else if (sys == 'C' && sat_num > 0 && sat_num < 47) {
                                BDS_PPP_STEC[epoch_index][sat_num] = value;
                            }
                        }
                    }
                }
            }


            auto save_stec_matrix = [&](const std::vector<std::vector<double>>& mat, const std::string& filename, const std::string& sys_prefix) {
                std::ofstream fout(filename);
                if (!fout.is_open()) {
                    my_logger->warn("Unable to open matrix file: {}", filename);
                    return;
                }

                int col_num = static_cast<int>(mat[0].size()) - 1;

                fout << std::setw(12) << "Epoch \\ PRN";
                for (int i = 1; i <= col_num; ++i) {
                    char prn_buf[12];
                    sprintf(prn_buf, "%s%02d", sys_prefix.c_str(), i); 
                    fout << std::setw(11) << prn_buf;
                }
                fout << "\n";

                for (size_t i = 1; i < mat.size(); ++i) {
                    char epoch_buf[20];
                    sprintf(epoch_buf, "Epoch %04d:", static_cast<int>(i));
                    fout << std::setw(12) << epoch_buf;

                    for (size_t j = 1; j < mat[i].size(); ++j) {
                        double val = mat[i][j];
                        if (fabs(val) < 1e-8) val = 0.0;
                        fout << std::setw(11) << std::fixed << std::setprecision(5) << val;
                    }
                    fout << "\n";
                }

                fout.close();
                my_logger->info("Matrix saved successfully: {}", filename);
                };


            //Convert SION to STEC before saving
            const double GPS_F1 = 1575.42e6;// Hz
            const double GAL_F1 = 1575.42e6;// Hz
            const double BDS_F1 = 1561.098e6;// Hz
           
            // GLONASS L1 frequency (Hz) for each satellite, 24 satellites
            double GLONASS_F1[24];

            for (int i = 0; i < 24; ++i) {
                GLONASS_F1[i] = (1602.0 + GLO_FCN[i] * 0.5625) * 1e6;
            }

            auto convert_sion_to_stec = [&](std::vector<std::vector<double>>& mat, const std::string& sys) {
                for (size_t i = 1; i < mat.size(); ++i) {       
                    for (size_t j = 1; j < mat[i].size(); ++j) { 
                        if (mat[i][j] != 0) {  
                            if (sys == "GPS") {
                                mat[i][j] = mat[i][j] * (GPS_F1 * GPS_F1) / (IONO_COEFF * 1e16);
                            }
                            else if (sys == "GAL") {
                                mat[i][j] = mat[i][j] * (GAL_F1 * GAL_F1) / (IONO_COEFF * 1e16);
                            }
                            else if (sys == "BDS") {
                                mat[i][j] = mat[i][j] * (BDS_F1 * BDS_F1) / (IONO_COEFF * 1e16);
                            }
                            else if (sys == "GLO") {
                                int sat_idx = j - 1;  
                                if (sat_idx >= 0 && sat_idx < 24) {
                                    double f1 = GLONASS_F1[sat_idx];
                                    mat[i][j] = mat[i][j] * (f1 * f1) / (IONO_COEFF * 1e16);
                                }
                            }
                        }
                    }
                }
            };

            convert_sion_to_stec(GPS_PPP_STEC, "GPS");
            convert_sion_to_stec(GLO_PPP_STEC, "GLO");
            convert_sion_to_stec(GAL_PPP_STEC, "GAL");
            convert_sion_to_stec(BDS_PPP_STEC, "BDS");

            int arc_edge_remove_len = 1;  
            int arc_min_length = gset.arc_min_length();

            remove_arc_edges_by_column(GPS_PPP_STEC, arc_edge_remove_len);
            remove_short_arcs_by_column(GPS_PPP_STEC, arc_min_length);

            remove_arc_edges_by_column(GLO_PPP_STEC, arc_edge_remove_len);
            remove_short_arcs_by_column(GLO_PPP_STEC, arc_min_length);

            remove_arc_edges_by_column(GAL_PPP_STEC, arc_edge_remove_len);
            remove_short_arcs_by_column(GAL_PPP_STEC, arc_min_length);

            remove_arc_edges_by_column(BDS_PPP_STEC, arc_edge_remove_len);
            remove_short_arcs_by_column(BDS_PPP_STEC, arc_min_length);


            if (sys_str.find("GPS") != std::string::npos) {
                save_stec_matrix(GPS_PPP_STEC, txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GPS.txt", "G");
                calc_ipp_matrix(OBS.X, OBS.Y, OBS.Z, SP3[1].X, SP3[1].Y, SP3[1].Z, 32, 2880, "GPS", txt_output_path);
                output_ppp_vtec_txt(txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GPS_VTEC.txt", "G", 32, GPS_PPP_STEC, OBS, &SP3[1], mappingFuncType,
                    [](const sp3* SP3ptr, int prn, int epoch, double& x, double& y, double& z) {
                        x = SP3ptr->X[prn][epoch]; y = SP3ptr->Y[prn][epoch]; z = SP3ptr->Z[prn][epoch];
                    }
                );
            }
            if (sys_str.find("GLO") != std::string::npos) {
                save_stec_matrix(GLO_PPP_STEC, txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GLO.txt", "R");
                calc_ipp_matrix(OBS.X, OBS.Y, OBS.Z, SP3[1].RX, SP3[1].RY, SP3[1].RZ, 24, 2880, "GLO", txt_output_path);
                output_ppp_vtec_txt(txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GLO_VTEC.txt", "R", 24, GLO_PPP_STEC, OBS, &SP3[1], mappingFuncType,
                    [](const sp3* SP3ptr, int prn, int epoch, double& x, double& y, double& z) {
                        x = SP3ptr->RX[prn][epoch]; y = SP3ptr->RY[prn][epoch]; z = SP3ptr->RZ[prn][epoch];
                    }
                );
            }
            if (sys_str.find("GAL") != std::string::npos) {
                save_stec_matrix(GAL_PPP_STEC, txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GAL.txt", "E");
                calc_ipp_matrix(OBS.X, OBS.Y, OBS.Z, SP3[1].EX, SP3[1].EY, SP3[1].EZ, 36, 2880, "GAL", txt_output_path);
                output_ppp_vtec_txt(txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_GAL_VTEC.txt", "E", 36, GAL_PPP_STEC, OBS, &SP3[1], mappingFuncType,
                    [](const sp3* SP3ptr, int prn, int epoch, double& x, double& y, double& z) {
                        x = SP3ptr->EX[prn][epoch]; y = SP3ptr->EY[prn][epoch]; z = SP3ptr->EZ[prn][epoch];
                    }
                );
            }
            if (sys_str.find("BDS") != std::string::npos) {
                save_stec_matrix(BDS_PPP_STEC, txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_BDS.txt", "C");
                calc_ipp_matrix(OBS.X, OBS.Y, OBS.Z, SP3[1].CX, SP3[1].CY, SP3[1].CZ, 46, 2880, "BDS", txt_output_path);
                output_ppp_vtec_txt(txt_output_path.substr(0, txt_output_path.find_last_of(".")) + "_BDS_VTEC.txt", "C", 46, BDS_PPP_STEC, OBS, &SP3[1], mappingFuncType,
                    [](const sp3* SP3ptr, int prn, int epoch, double& x, double& y, double& z) {
                        x = SP3ptr->CX[prn][epoch]; y = SP3ptr->CY[prn][epoch]; z = SP3ptr->CZ[prn][epoch];
                    }
                );
            }

        }
        delete pppflt;
    }

    delete gobs;
    delete gorb;
    if (gpcv) delete gpcv;
    if (gotl) delete gotl;
    if (gobj) delete gobj;
    if (gbia) delete gbia;
    if (gde) delete gde;
    if (gerp) delete gerp;
    if (gupd) delete gupd;
    if (gifcb) delete gifcb;
    if (data) delete data;

    my_logger->info("PPP_STEC function completed.");
}

void run_roti(t_gcfg_ppp& gset, std::shared_ptr<spdlog::logger> my_logger)
{
    my_logger->info("Starting ROTI function...");

    std::string raw_txt_path = gset.outputs("txt");
    if (raw_txt_path.find("file://") == 0)
        raw_txt_path = raw_txt_path.substr(7);
    if (raw_txt_path.empty())
        raw_txt_path = "result/ROTI_result/$(rec)_$(function).txt";

    std::string function_name = gset.outputs("function");
    if (function_name.empty()) function_name = "ROTI";

    auto replace_all = [](std::string& str, const std::string& from, const std::string& to) {
        size_t pos = 0;
        while ((pos = str.find(from, pos)) != std::string::npos) {
            str.replace(pos, from.length(), to);
            pos += to.length();
        }
    };

    double elevCutoffDeg = gset.minimum_elev();
    std::string sys_str;
    for (const auto& sys : gset.sys()) sys_str += sys + " ";

    t_gallobs* gobs = decode_obs(gset, my_logger);
    set<string> stations = dynamic_cast<t_gsetgen*>(&gset)->recs();

    for (const auto& station : stations)
    {
        auto runepoch = t_gtime::current_time(t_gtime::GPS); 
        obs OBS;
        string stationName = station;
        std::vector<t_gtime> epochs;

        my_logger->info("Start processing station: {}", station);

        std::string txt_output_path = raw_txt_path;
        replace_all(txt_output_path, "$(rec)", station);
        replace_all(txt_output_path, "$(function)", function_name);
        std::filesystem::create_directories(std::filesystem::path(txt_output_path).parent_path());

        std::vector<double> xyz = gobs->get_station_coordinates(station);
        OBS.X = xyz[0], OBS.Y = xyz[1], OBS.Z = xyz[2];

        std::vector<std::string> filenamesSp3 = getSortedSp3Paths(gset);
        sp3 SP3[3];
        for (int i = 0; i < 3; i++) {
            std::ifstream file(filenamesSp3[i]);
            if (file.is_open()) getSp3Data(file, filenamesSp3[i], SP3[i]);
        }
        processSP3(SP3);

        if (sys_str.find("GPS") != std::string::npos) {
            OBS.reset();
            std::vector<std::vector<double>> C1, C2, L1, L2,GPS_S1,GPS_S2;
            std::vector<std::string> sats;
            std::vector<std::string> GPS_sats;
            bool isC1WAllZero = true;
            extract_GPS_obs(gobs, station, C1, C2, L1, L2, epochs, sats, isC1WAllZero, my_logger, OBS);
            //extract_GPS_SNR(gobs, station, GPS_S1, GPS_S2, epochs, GPS_sats, my_logger, OBS);
            apply_elevation_mask(OBS, SP3[1], elevCutoffDeg, 32, 'G');
            calc_roti_GPS(OBS, stationName, SP3[1], isC1WAllZero, txt_output_path, gset);
            //calc_S4C(OBS, 32, 2880, 10, "S4C_S1.txt", "S4C_S2.txt");
        }

        if (sys_str.find("BDS") != std::string::npos) {
            OBS.reset();
            std::vector<std::vector<double>> C2, C6, C7, L2, L6, L7;
            std::vector<std::string> sats;
            bool isC7IAllZero = true;
            extract_BDS_obs(gobs, station, C2, C6, C7, L2, L6, L7, epochs, sats, OBS, my_logger, isC7IAllZero);
            apply_elevation_mask(OBS, SP3[1], elevCutoffDeg, 46, 'C');
            calc_roti_BDS(OBS, stationName, SP3[1], isC7IAllZero, txt_output_path, gset);
        }

        if (sys_str.find("GLO") != std::string::npos) {
            OBS.reset();
            std::vector<std::vector<double>> C1, C2, L1, L2;
            std::vector<std::string> sats;
            extract_GLO_obs(gobs, station, C1, C2, L1, L2, epochs, sats, OBS, my_logger);
            apply_elevation_mask(OBS, SP3[1], elevCutoffDeg, 24, 'R');
            calc_roti_GLO(OBS, stationName, SP3[1], txt_output_path, gset);
        }

        if (sys_str.find("GAL") != std::string::npos) {
            OBS.reset();
            std::vector<std::vector<double>> C1, C5, L1, L5;
            std::vector<std::string> sats;
            bool isC1QAllZero = false; 
            extract_GAL_obs(gobs, station, C1, C5, L1, L5, epochs, sats, OBS, my_logger, isC1QAllZero);
            apply_elevation_mask(OBS, SP3[1], elevCutoffDeg, 36, 'E');
            calc_roti_GAL(OBS, stationName, SP3[1], txt_output_path, gset);
        }

        
        auto lstepoch = t_gtime::current_time(t_gtime::GPS); 
        my_logger->info("Station {} ROTI output completed", station);
        my_logger->info("Station {} ROTI processing completed in {:.3f} seconds", station, lstepoch.diff(runepoch));

    }

    delete gobs;
}


void run_aatr(t_gcfg_ppp& gset, std::shared_ptr<spdlog::logger> my_logger)
{
    my_logger->info("Starting AATR function...");

    std::string raw_txt_path = gset.outputs("txt");
    if (raw_txt_path.find("file://") == 0) raw_txt_path = raw_txt_path.substr(7);
    if (raw_txt_path.empty()) raw_txt_path = "result/AATR_result/$(rec)_$(function).txt";

    std::string function_name = gset.outputs("function");
    if (function_name.empty()) function_name = "AATR";

    auto replace_all = [](std::string& str, const std::string& from, const std::string& to) {
        size_t pos = 0;
        while ((pos = str.find(from, pos)) != std::string::npos) {
            str.replace(pos, from.length(), to);
            pos += to.length();
        }
    };

    double elevCutoffDeg = gset.minimum_elev();
    std::string sys_str;
    for (const auto& sys : gset.sys()) sys_str += sys + " ";

    t_gallobs* gobs = decode_obs(gset, my_logger);
    set<string> stations = dynamic_cast<t_gsetgen*>(&gset)->recs();

    std::vector<std::string> filenamesSp3 = getSortedSp3Paths(gset);
    sp3 SP3[3];
    for (int i = 0; i < 3; i++) {
        std::ifstream file(filenamesSp3[i]);
        if (file.is_open()) {
            read_Time_sp3(file, filenamesSp3[i], SP3[i]);
            getSp3Data(file, filenamesSp3[i], SP3[i]);
        }
    }
    processSP3(SP3);

    for (const auto& station : stations) {
        auto runepoch = t_gtime::current_time(t_gtime::GPS);
        obs OBS;
        std::string stationName = station;

        std::vector<t_gtime> epochs;
        std::vector<std::string> sats;

        std::string txt_output_path = raw_txt_path;
        replace_all(txt_output_path, "$(rec)", station);
        replace_all(txt_output_path, "$(function)", function_name);
        std::filesystem::create_directories(std::filesystem::path(txt_output_path).parent_path());

        // Set receiver coordinates
        std::vector<double> xyz = gobs->get_station_coordinates(station);
        OBS.X = xyz[0], OBS.Y = xyz[1], OBS.Z = xyz[2];

        if (sys_str.find("GPS") != std::string::npos) {
            std::vector<std::vector<double>> C1, C2, L1, L2;
            bool isC1WAllZero = true;
            extract_GPS_obs(gobs, station, C1, C2, L1, L2, epochs, sats, isC1WAllZero, my_logger, OBS);
            apply_elevation_mask(OBS, SP3[1], elevCutoffDeg, 32, 'G');
            calc_aatr_GPS(OBS, station, SP3, txt_output_path, gset);
        }

        if (sys_str.find("BDS") != std::string::npos) {
            OBS.reset(); OBS.X = xyz[0], OBS.Y = xyz[1], OBS.Z = xyz[2];
            std::vector<std::vector<double>> C2, C6, C7, L2, L6, L7;
            bool isC7IAllZero = true;
            extract_BDS_obs(gobs, station, C2, C6, C7, L2, L6, L7, epochs, sats, OBS, my_logger, isC7IAllZero);
            apply_elevation_mask(OBS, SP3[1], elevCutoffDeg, 46, 'C');
            calc_aatr_BDS(OBS, station, isC7IAllZero, SP3, txt_output_path, gset);
        }

        if (sys_str.find("GLO") != std::string::npos) {
            OBS.reset(); OBS.X = xyz[0], OBS.Y = xyz[1], OBS.Z = xyz[2];
            std::vector<std::vector<double>> C1, C2, L1, L2;
            extract_GLO_obs(gobs, station, C1, C2, L1, L2, epochs, sats, OBS, my_logger);
            apply_elevation_mask(OBS, SP3[1], elevCutoffDeg, 24, 'R');
            calc_aatr_GLO(OBS, station, SP3, txt_output_path, gset);
        }

        if (sys_str.find("GAL") != std::string::npos) {
            OBS.reset(); OBS.X = xyz[0], OBS.Y = xyz[1], OBS.Z = xyz[2];
            std::vector<std::vector<double>> C1, C5, L1, L5;
            bool isC1QAllZero = false;
            extract_GAL_obs(gobs, station, C1, C5, L1, L5, epochs, sats, OBS, my_logger, isC1QAllZero);
            apply_elevation_mask(OBS, SP3[1], elevCutoffDeg, 36, 'E');
            calc_aatr_GAL(OBS, station, SP3, txt_output_path, gset);
        }

        auto lstepoch = t_gtime::current_time(t_gtime::GPS);
        my_logger->info("Station {} AATR output completed", station);
        my_logger->info("Station {} AATR processing completed in {:.3f} seconds", station, lstepoch.diff(runepoch));
    }

    delete gobs;
}

