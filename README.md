# IonoMoni: Ionospheric Parameter Monitoring and Analysis Software 

## Overview

Ionospheric delay is one of the primary error sources affecting Global Navigation Satellite Systems (GNSS). In addition, ionospheric parameters such as Total Electron Content (TEC) are widely applied in various fields, including space weather monitoring, earthquake detection, and radio signal propagation studies. Accurate extraction and analysis of these ionospheric parameters are essential for both scientific research and practical applications. To meet this need, we developed the IonoMoni ionospheric parameter computation software, based on the secondary development of the GREAT-PVT software within the Visual Studio 2022 environment. The IonoMoni program offers multiple functionalities, including:

## Key Features

- Supports STEC extraction using the dual-frequency carrier-to-code leveling (DFCCL) method
- Supports STEC extraction based on the undifferenced and uncombined precise point positioning (UCPPP) method
- Supports the calculation of the Rate of TEC Index (ROTI), a widely used indicator for detecting ionospheric irregularities and scintillation phenomena
- Supports the estimation of the Along Arc TEC Rate (AATR), an effective metric for monitoring ionospheric disturbances, especially during geomagnetic storms or in equatorial and polar regions
- Supports processing of multi-GNSS observational data, including GPS, BDS, Galileo, and GLONASS constellations
- Compatible with both Rinex 2.x and Rinex 3.x file formats
- Includes batch processing and plotting capabilities for efficient extraction and visualization of ionospheric parameters across multiple stations and days

## Directory Structure

```
./bin                   The executable binary applications for Windows *
./src                   IonoMoni source programs *
./include               IonoMoni header files *
./lib                   Static libraries and export files *
./build_resources       Dynamic link libraries and sample data *
./batch_process         Batch processing Python scripts *
./plot                  Plotting Python scripts *
./poleut1               Earth Orientation Parameters (EOP) generation program *
./doc                   Document files *
  IonoMoni_User_Manual.pdf    User manual *
  sample_config.xml           Sample XML configuration file *
./CMakeLists.txt         CMake build configuration file *
./CMakePresets.json      CMake build preset settings *
```



## Installation and Usage



For detailed installation and usage instructions, please refer to the **IonoMoni_manual_\<ver\>.pdf** document included in the `./doc` directory.
<br>
QQ交流群：1049702653
