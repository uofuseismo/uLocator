#include <cmath>
namespace
{
/*
!>    This is a C++ version of the BASIC program "Transverse Mercator
!>    Conversion", Copyright 1986, Norman J. Berls (Stefan Musarra, 2/94)
!>    Based on algorithm taken from "Map Projections Used by the USGS"
!>    by John P. Snyder, Geological Survey Bulletin 1532, USDI.
!>
!>     2016/03/24   Made it C compatible, moved hard to find constants inside,
!>                  introduced handling for doxygen, and ensured code does not 
!>                  go past column 90
!
!     Input/Output arguments:
!
!>    @param[in,out] rlon4              Longitude (degrees, negative for West)
!>    @param[in,out] rlat4              Latitude (degrees)
!>    @param[in,out] rx4                UTM easting (meters)
!>    @param[in,out] ry4                UTM northing (meters)
!>    @param[in] utmProjectionZone   UTM zone
!>                                   The Northern hemisphere corresponds to zones +1 to +60
!>                                   The Southern hemisphere corresponds to zones -1 to -60
!>    @param[in] iway                Conversion type
!>                                   ILONGLAT2UTM (0) = geodetic to UTM
!>                                   IUTM2LONGLAT (1) = UTM to geodetic
!>
!>    @reference Some general information about UTM:
!>               (for more details see e.g. http://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system )
!
!>    @note There are 60 longitudinal projection zones numbered 1 to 60 starting at 180 degrees W.
!>          Each of these zones is 6 degrees wide, apart from a few exceptions around Norway and Svalbard.
!>          There are 20 latitudinal zones spanning the latitudes 80 degrees S to 84 degrees N and denoted
!>          by the letters C to X, ommitting the letter O.
!>          Each of these is 8 degrees south-north, apart from zone X which is 12 degrees south-north.
!>
!>          The UTM zone is described by the central meridian of that zone, i.e. the longitude at the
!>          midpoint of the zone, 3 degrees away from both zone boundary.
!>
!>    author Dimitri Komatitsch and Jeroen Tromp
!>
! From http://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system :
! The Universal Transverse Mercator coordinate system was developed by the United States Army Corps of Engineers in the 1940s.
! The system was based on an ellipsoidal model of Earth. For areas within the contiguous United States
! the Clarke Ellipsoid of 1866 was used. For the remaining areas of Earth, including Hawaii, the International Ellipsoid was used.
! The WGS84 ellipsoid is now generally used to model the Earth in the UTM coordinate system,
! which means that current UTM northing at a given point can be 200+ meters different from the old one.
! For different geographic regions, other datum systems (e.g.: ED50, NAD83) can be used.
*/
void utmGeo(double *rlon4, double *rlat4,
            double *rx4,   double *ry4,
            const int utmProjectionZone,
            const bool toUTM,
            const bool suppressUTMProjection)
{
    // WGS84
    constexpr double SEMI_MAJOR_AXIS{6378137.0};
    constexpr double SEMI_MINOR_AXIS{6356752.314245};
    // Note that the UTM grids are actually Mercators which
    // employ the standard UTM scale factor 0.9996 and set the Easting Origin
    // to 500,000.
    constexpr double DEGREES_TO_RADIANS{0.017453292519943295};
    constexpr double RADIANS_TO_DEGREES{57.29577951308232};
    constexpr double scfa{0.9996};
    constexpr double north{0};
    constexpr double east{500000.};
  
    // Checks if conversion to UTM has to be done
    if (suppressUTMProjection)
    {
        if (toUTM)
        {
            *rx4 = *rlon4;
            *ry4 = *rlat4;
        }
        else
        {
            *rlon4 = *rx4;
            *rlat4 = *ry4;
        } 
        return;
    }

    // save original parameters
    auto rlon_save = *rlon4;
    auto rlat_save = *rlat4;
    auto rx_save = *rx4;
    auto ry_save = *ry4;

    const double e2{1.0 - std::pow(SEMI_MINOR_AXIS/SEMI_MAJOR_AXIS, 2)};
    const double e4{e2*e2};
    const double e6{e2*e4};
    const double ep2{2/(1.0 - e2)};

    // Set Zone parameters
    bool lsouth = false;
    if (utmProjectionZone < 0){lsouth = true;}
    int zone = std::abs(utmProjectionZone);
    double cm = zone*6.0 - 183.0;
    double cmr = cm*DEGREES_TO_RADIANS;

    double dlat = 0.0; // ! removes an uninitialized warning - baker
    // Lat/Lon to UTM conversion
    if (toUTM)
    {
        double dlat = *rlat4;// ! removes an uninitialized warning - baker
        double dlon = *rlon4;// ! removes an uninitialized warning - baker

        //double rlon = DEGREES_TO_RADIANS*dlon;
        double rlat = DEGREES_TO_RADIANS*dlat;

        double delam = dlon - cm;
        if (delam < -180.0){delam = delam + 360.0;}
        if (delam > 180.0){delam = delam - 360.0;}
        delam = delam*DEGREES_TO_RADIANS;

        double f1 = (1.0 - e2/4.0 - 3.0*e4/64.0 - 5.0*e6/256.0)*rlat;
        double f2 = 3.0*e2/8.0 + 3.0*e4/32.0 + 45.0*e6/1024.0;
        f2 = f2*std::sin(2.0*rlat);
        double f3 = 15.0*e4/256.0*45.0*e6/1024.0;
        f3 = f3*std::sin(4.0*rlat);
        double f4 = 35.0*e6/3072.0;
        f4 = f4*std::sin(6.0*rlat);
        double rm = SEMI_MAJOR_AXIS*(f1 - f2 + f3 - f4);
        double xx = 0;
        double yy = scfa*rm;
        if (std::abs(dlat - 90) >= 1.e-12 && std::abs(dlat + 90) >= 1.e-12)
        {
            double rn = SEMI_MAJOR_AXIS
                       /std::sqrt(1.0 - e2*std::pow(sin(rlat), 2));
            double t = std::pow(std::tan(rlat), 2);
            double c = ep2*std::pow(std::cos(rlat), 2);
            double a = std::cos(rlat)*delam;

            f1 = (1.0 - t + c)*std::pow(a, 3)/6.0;
            f2 = 5.0 - 18.0*t + std::pow(t, 2) + 72.0*c - 58.0*ep2;
            f2 = f2*std::pow(a, 5)/120.0;
            xx = scfa*rn*(a + f1 + f2);
            f1 = std::pow(a, 2)/2.0;
            f2 = 5.0 - t + 9.0*c + 4.0*std::pow(c, 2);
            f2 = f2*std::pow(a, 4)/24.0;
            f3 = 61.0 - 58.0*t + std::pow(t,2) + 600.0*c - 330.0*ep2;
            f3 = f3*std::pow(a, 6)/720.0;
            yy = scfa*(rm + rn*std::tan(rlat)*(f1 + f2 + f3));
        }
        xx = xx + east;
        yy = yy + north;
        // Set output
        *rx4 = xx;
        if (lsouth){yy = yy + 1.e7;}
        *ry4 = yy;
        *rlon4 = rlon_save;
        *rlat4 = rlat_save;
    }
    else // UTM to lat/lon
    {
        double xx = *rx4;//  ! removes uninitialized warning - baker
        double yy = *ry4;//  ! removes uninititalted warning - baker
        if (lsouth){yy = yy - 1.e7;}// ! removes uninitialized warning - baker

        xx = xx - east;
        yy = yy - north;
        double e1 = std::sqrt(1.0 - e2);
        e1 = (1.0 - e1)/(1.0 + e1);
        double rm = yy/scfa;
        double u = 1.0 - e2/4.0 - 3.0*e4/64.0 - 5.0*e6/256.0;
        u = rm/(SEMI_MAJOR_AXIS*u);

        double f1 = 3.0*e1/2.0 - 27.0*std::pow(e1, 3)/32.0;
        f1 = f1*std::sin(2.0*u);
        double f2 = 21.0*std::pow(e1, 2)/16.0 - 55.0*std::pow(e1, 4)/32.0;
        f2 = f2*std::sin(4.0*u);
        double f3 = 151.0*pow(e1, 3)/96.0;
        f3 = f3*std::sin(6.0*u);
        double rlat1 = u + f1 + f2 + f3;
        double dlat1 = rlat1*RADIANS_TO_DEGREES;
        double dlon = cm;
        if (dlat1 >= 90.0 || dlat1 <= -90.0)
        {
            dlat1 = std::fmin(dlat1, 90.0);
            dlat1 = std::fmax(dlat1,-90.0);
            //dlon = cm;
        }
        else
        {
            double c1 = ep2*std::pow(cos(rlat1), 2);
            double t1 = std::pow(std::tan(rlat1), 2);
            f1 = 1.0 - e2*std::pow(std::sin(rlat1), 2);
            double rn1 = SEMI_MAJOR_AXIS/sqrt(f1);
            double r1 = SEMI_MAJOR_AXIS*(1.0 - e2)/std::sqrt(std::pow(f1, 3));
            double d = xx/(rn1*scfa);

            f1 = rn1*std::tan(rlat1)/r1;
            f2 = std::pow(d,2)/2.0;
            f3 = 5.0*3.0*t1 + 10.0*c1 - 4.0*std::pow(c1, 2) - 9.0*ep2;
            f3 = f3*std::pow(d, 2)*std::pow(d, 2)/24.0;
            double f4 = 61.0 + 90.0*t1 + 298.0*c1
                      + 45.0*std::pow(t1, 2) - 252.0*ep2 - 3.0*std::pow(c1, 2);
            f4 = f4*std::pow((std::pow(d, 2)),3.0)/720.0;
            //double rlat = rlat1 - f1*(f2 - f3 + f4);
            //double dlat = rlat*RADIANS_TO_DEGREES;

            f1 = 1.0 + 2.0*t1 + c1;
            f1 = f1*std::pow(d, 2)*d/6.;
            f2 = 5.0 - 2.0*c1 + 28.0*t1 - 3.0*std::pow(c1,2)
               + 8.0*ep2 + 24.0*std::pow(t1, 2);
            f2 = f2*std::pow((std::pow(d, 2)), 2)*d/120.0;
            double rlon = cmr + (d - f1 + f2)/std::cos(rlat1);
            dlon = rlon*RADIANS_TO_DEGREES;
            if (dlon < -180.0){dlon = dlon + 360.0;}
            if (dlon > 180.0){dlon = dlon - 360.0;}
        }
        // Set output
        *rlon4 = dlon;
        *rlat4 = dlat;
        *rx4 = rx_save;
        *ry4 = ry_save;
    }
}

/// @brief Converts latitudes/longitudes to UTMs
/// @param[in] latitude            latitude (degrees)
/// @param[in] longitude           longitude (degrees, negative for west)
/// @param[in] utmProjectionZone   UTM zone
///                                The Northern hemisphere corresponds to
///                                zones +1 to +60.
///                                The Southern hemisphere corresponds to
///                                zones -1 to -60
/// @param[out] easting            UTM easting (m)
/// @param[out] northing           UTM northing (m)
void latitudeLongitudeToUTM(const double latitude, const double longitude,
                            const int utmProjectionZone,
                            double *easting, double *northing)
{
    constexpr bool toUTM{true}; // (lat, lon) -> UTM
    constexpr bool suppressUTMProjection{false};
    double rlon4_use, rlat4_use;
    double rx4, ry4;
    rlon4_use = longitude;
    rlat4_use = latitude;
    ::utmGeo(&rlon4_use, &rlat4_use,
             &rx4, &ry4, utmProjectionZone,
             toUTM, suppressUTMProjection);
    *easting = rx4;
    *northing = ry4;
}
/// @brief Converts UTMs to latitudes/longitudes
/// @param[in] easting            UTM easting (m)
/// @param[in] northing           UTM northing (m)
/// @param[in] utmProjectionZone  UTM zone
///                               The Northern hemisphere corresponds to
///                               zones +1 to +60.
///                               The Southern hemisphere corresponds to
///                               zones -1 to -60.
/// @param[out] latitude          corresponding latitude (degrees)
/// @param[out] longitude         corresponding longitude (degrees)
void utmToLatitudeLongitude(const double easting, const double northing,
                            const  int utmProjectionZone,
                            double *latitude, double *longitude)
{
    constexpr bool toUTM{false}; // UTM -> (lat, lon)
    constexpr bool suppressUTMProjection{false};
    double rx4_use, ry4_use;
    double rlat4, rlon4;
    rx4_use = easting;
    ry4_use = northing;
    ::utmGeo(&rlon4, &rlat4, &rx4_use, &ry4_use, utmProjectionZone,  
             toUTM, suppressUTMProjection);
    *latitude  = rlat4;
    *longitude = rlon4;
}
}
