#include "decode_obs.h"
#include "gcoders/rinexo.h"
#include "gio/gfile.h"
#include <iostream>
using namespace std;
using namespace gnut;

namespace gnut {

    // Function to decode all input GNSS observation files and store results in a t_gallobs object
    t_gallobs* decode_obs(t_gcfg_ppp& gset, shared_ptr<spdlog::logger> my_logger)
    {
        // Log the start of the decoding process
        my_logger->info("Start decoding observation files...");

        // Allocate and initialize the observation container
        t_gallobs* gobs = new t_gallobs();
        gobs->gset(&gset);
        gobs->spdlog(my_logger);

        // Retrieve all input file paths of RINEXO_INP (RINEX observation) type
        multimap<IFMT, string> inputs = gset.inputs_all();
        int obs_count = 0;  // Counter for successfully decoded files

        // Iterate through all input files and decode only observation files
        for (auto it = inputs.begin(); it != inputs.end(); ++it)
        {
            IFMT fmt = it->first;
            string filename = it->second;

            // Skip non-observation files
            if (fmt != IFMT::RINEXO_INP) continue;

            // Log the file being processed
            my_logger->info("Decoding observation file: {}", filename);

            // Create decoder and set file path
            t_rinexo* decoder = new t_rinexo(&gset);
            decoder->path(filename);
            decoder->add_data("ID0", gobs);      // Associate output container
            decoder->spdlog(my_logger);          // Pass logger for consistent logging

            // Create file reader, assign decoder, and run decoding
            t_gfile* file_reader = new t_gfile(my_logger);
            file_reader->path(filename);
            file_reader->coder(decoder);
            file_reader->run_read();             // Perform file reading and decoding

            // Free resources
            delete decoder;
            delete file_reader;

            obs_count++;
        }

        // Log total number of decoded observation files
        my_logger->info("Number of successfully decoded observation files: {}", obs_count);

        // Return pointer to filled t_gallobs object
        return gobs;
    }

}
