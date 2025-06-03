# IonoMoni: Ionospheric monitoring based on dual-frequency data from single GNSS station v1.0 by Henan University (last update: 2025-06-03)

## Overview

The ionosphere has a significant impact on Global Navigation Satellite System (GNSS) signals, introducing refraction and diffraction effects that can severely degrade positioning performance. With the current peak of Solar Cycle 25, increasingly prominent ionospheric space weather events and disturbances have further intensified adverse effects on GNSS applications. On the other hand, dual-frequency GNSS observations enable effective monitoring of ionospheric conditions, serving as a crucial foundation for the mitigation of ionospheric effects on GNSS positioning and the advancement of research into the physical mechanisms of ionospheric variability. The extraction and analysis of ionospheric parameters are essential for both scientific research and practical applications. In this context, we developed the IonoMoni ionospheric parameter monitoring software, based on the Microsoft Visual Studio 2022 integrated development environment and GREAT-PVT platform, using the C++ programming language and following the C++17 standard. The project adopts CMake as the build system, which is primarily used to organize the source code structure and manage the linking of external libraries. CMake integrates seamlessly with Visual Studio, enabling automated compilation and debugging, and also facilitates future migration to Linux or other platforms. In addition, IonoMoni uses XML files for parameter configuration, which not only allows flexible control over operation modes, input paths, and functional modules, but also enables seamless integration with batch processing scripts for efficient automated workflows.

## Key Features

- Supports STEC extraction using the dual-frequency carrier-to-code leveling (CCL) method
- Supports STEC extraction based on the undifferenced and uncombined precise point positioning (UCPPP) method
- Supports VTEC conversion based on ionospheric mapping function and STEC
- Supports the calculation of the Rate of TEC Index (ROTI), a widely used indicator for detecting ionospheric irregularities of ionospheric diffractive effects
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

## Environmental Requirements
### Developer Environment
The executable CUI AP for Windows in the package was built using Microsoft Visual Studio 2022 on Windows 11 (64-bit), with the ReleaseWithXml configuration. All required dynamic link libraries are included in the directory.
### License
IonoMoni is released as open-source software under the GNU General Public License, version 3 (GPLv3). This license permits anyone to use, modify, and redistribute the software, provided that the same license terms are preserved and the source code remains accessible.
For detailed information, please refer to the full license text at: https://www.gnu.org/licenses/gpl-3.0.html.
### Python version
IonoMoni includes several auxiliary Python scripts that require a Python environment and a set of dependencies listed in the requirements.txt file. For optimal compatibility, Python version 3.13 is recommended. It is further advised to use a virtual environment manager, such as Conda or venv, to ensure reproducibility and isolate dependencies.

## Installation and Usage

For detailed installation and usage instructions, please refer to the **IonoMoni_manual_\<ver\>.pdf** document included in the `./doc` directory.
<br>

## Visualization
The Python plotting scripts provided with IonoMoni support multi-station plotting for a single day. The directory structure under the plot folder is as follows:
```
./data_aatr      Input folder for AATR plot data
./data_ccl       Input folder for CCL_STEC plot data
./data_ppp       Input folder for PPP_STEC plot data
./data_roti      Input folder for ROTI plot data
./output         Output folder for plot results
./PythonScripts  Plot script folder
```

## Contributing
### Developer
Yuzhi Zhang, Yihong Zhang, Qi Liu, Manuel Hernández-Pajares, Enric Monte-Moreno
<br>
Research group of Henan University, UPC-Ionsat, Zhengzhou University.

### Third-Party Libraries
GREAT-PVT library

G-Nut Library (http://www.pecny.cz)
Copyright (C) 2011-2016 GOP - Geodetic Observatory Pecny, RIGTC.

pugixml Library (http://pugixml.org)
Copyright (C) 2006-2014 Arseny Kapoulkine.

Newmat Library (http://www.robertnz.net/nm_intro.htm)
Copyright (C) 2008: R B Davies.

spdlog Library (https://github.com/gabime/spdlog)
Copyright (c) 2015-present, Gabi Melman & spdlog contributors.

Eigen Library (https://eigen.tuxfamily.org)
Copyright (C) 2008-2011 Gael Guennebaud.

FAST Source Code (https://github.com/ChangChuntao/FAST)
Copyright (C) The GNSS Center, Wuhan University & Chinese Academy of Surveying and Mapping.

## Getting in touch
You can contact us for bug reports and comments by sending an email or leaving a message on our website:
<br>
Email: pride@whu.edu.cn
<br>
For Chinese users, we provide Tencent QQ Group service.
<br>
QQ group：1049702653
