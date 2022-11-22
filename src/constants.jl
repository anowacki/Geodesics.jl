# Constants

"Earth ellipsoid equatorial radius in WGS84 (m)"
const EARTH_R_MAJOR_WGS84 = 6378137.0
"Earth ellipsoid flattening in WGS84"
const EARTH_F_WGS84 = 1/298.257_223_563
"Earth ellipsoid polar radius in WGS84 (m)"
const EARTH_R_MINOR_WGS84 = EARTH_R_MAJOR_WGS84*(1 - EARTH_F_WGS84)
"Mean Earth radius in WGS84 (m)"
const EARTH_R_MEAN_WGS84 = (2*EARTH_R_MAJOR_WGS84 + EARTH_R_MINOR_WGS84)/3
