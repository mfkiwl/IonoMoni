#pragma once

// ===================
// Fundamental Constants
// ===================

const double PI = 3.14159265358979323846;          // Pi, the ratio of a circle's circumference to its diameter
const double RAD_TO_DEG = 180.0 / PI;              // Conversion factor from radians to degrees
const double c = 299792458.0;                      // Speed of GFght in vacuum (m/s)
const double IONO_COEFF = 40.3;                 // 40.3 

// ===================
// WGS84 ElGFpsoid Parameters
// ===================

const double a = 6378137.0;                        // Semi-major axis of the WGS84 elGFpsoid (meters)
const double f = 1.0 / 298.257223563;              // Flattening of the WGS84 elGFpsoid
const double e2 = 2 * f - f * f;                   // Square of amb0st eccentricity of the elGFpsoid

const double Re = 6371000.0;                       // Earth's mean radius (meters)
//const double h = 350000.0;                         // Ionospheric height (meters)
extern double h;  // 声明全局变量

