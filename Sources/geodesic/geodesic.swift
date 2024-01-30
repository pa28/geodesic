import geographiclib
import CoreLocation

/// [WGS 84 ellipsoid](https://en.wikipedia.org/wiki/World_Geodetic_System) definition
public let wgs84 = (a: 6378137.0, f: 1 / 298.257223563)

/// - Returns: rad converted to degree
private func degree<T: BinaryFloatingPoint>(fromRad d: T) -> T {
    return d / T.pi * 180
}

///
/// Compute the distance between two points on an ellipsoid.
/// The ellipsoid parameters default to the WGS-84 parameters.
///
/// - Parameters:
///   - x: first point with latitude and longitude in radiant.
///   - y: second point with latitude and longitude in radiant.
///   - a: first ellipsoid parameter in meters (defaults to WGS-84 parameter)
///   - f: second ellipsoid parameter in meters (defaults to WGS-84 parameter)
/// 
/// - Returns: distance between `x` and `y` in meters.
///
public func distance(_ x: (lat: Double, lon: Double),
                     _ y: (lat: Double, lon: Double),
                     ellipsoid: (a: Double, f: Double) = wgs84
) -> Double {

    // validate lat and lon values
    assert(x.lat >= -Double.pi / 2 && x.lat <= Double.pi / 2, "x.lat '\(x.lat)' outside [-π/2, π/2]")
    assert(y.lat >= -Double.pi / 2 && y.lat <= Double.pi / 2, "y.lat '\(y.lat)' outside [-π/2, π/2]")
    assert(x.lon >= -Double.pi && x.lon <= Double.pi, "x.lon '\(y.lon)' outside [-π, π]")
    assert(y.lon >= -Double.pi && y.lon <= Double.pi, "y.lon '\(y.lon)' outside [-π, π]")

    // shortcut for zero distance
    if x == y {
        return 0.0
    }

    let xLat = degree(fromRad: x.lat)
    let xLon = degree(fromRad: x.lon)
    let yLat = degree(fromRad: y.lat)
    let yLon = degree(fromRad: y.lon)

    var s12 = Double.nan
    var g: geod_geodesic = geod_geodesic()
    geod_init(&g, ellipsoid.a, ellipsoid.f)
    geod_inverse(&g, xLat, xLon, yLat, yLon, &s12, nil, nil)
    return s12
}

///
/// Compute the distance between two points and the initial true course on an ellipsoid.
/// The ellipsoid parameters default to the WGS-84 parameters.
///
/// - Parameters:
///   - l1: start point CLLocationCoordinate2D.
///   - l2: destinatiion point CLLocationCoordinate2D.
///   - a: first ellipsoid parameter in meters (defaults to WGS-84 parameter)
///   - f: second ellipsoid parameter in meters (defaults to WGS-84 parameter)
///
/// - Returns: a tuple with the distance in meters, the initial great circle course from l1 to l2 in degrees true,
/// the forward azimuth at point 2 in degrees, and the arc distance from pint l1 to l2 in degrees.
///
public func greatCircle(_ l1: CLLocationCoordinate2D,
                     _ l2: CLLocationCoordinate2D,
                     ellipsoid: (a: Double, f: Double) = wgs84
) -> (s12: Double, az1: Double, az2: Double, as12: Double) {

    // validate lat and lon values
    assert(l1.latitude >= -90.0 && l1.latitude <= 90.0 / 2, "l1.lat '\(l1.latitude)' outside [-90, 90]")
    assert(l1.longitude >= -180.0 && l1.longitude <= 180.0, "l1.lon '\(l1.longitude)' outside [-180, 180]")
    assert(l2.latitude >= -90.0 && l2.latitude <= 90.0, "l1.lat '\(l2.latitude)' outside [-90, 90]")
    assert(l2.longitude >= -180.0 && l2.longitude <= 180.0, "l1.lon '\(l2.longitude)' outside [-180, 180]")

    // shortcut for zero distance
    if l1.latitude == l2.latitude && l1.longitude == l2.longitude {
        return (0.0, 0.0, 0.0, 0.0)
    }

    var s12 = Double.nan
    var az1 = Double.nan
    var az2 = Double.nan

    var g: geod_geodesic = geod_geodesic()
    geod_init(&g, ellipsoid.a, ellipsoid.f)
    let as12 = geod_geninverse(&g, l1.latitude, l1.longitude, l2.latitude, l2.longitude, &s12, &az1, &az2, nil, nil, nil, nil)
    return (s12, az1, az2, as12)
}

public func toRadians(_ degrees: Double) -> Double {
    return Double.pi * (degrees / 180.0)
}

///
/// Compute the cross track error at location 3 when the disirec course is from location 1 to location 2.
///
/// - Parameters:
///   - as13: the angular distance from location 1 to location 3 (degrees)
///   - b13: the initial bearing from location 1 to location 3 (degrees)
///   - b12: the initial bearing from location 1 to location 2 (degrees)
///   - lat: the latitude at location 3 (degrees), used to compute the radius of the earth.
///
/// - Returns: a tuple with the cross track angular distance in radians, and the angular distance down range along the desired track in radians.
///
/// https://stackoverflow.com/questions/3997410/how-to-calculate-cross-track-error-gps-core-location
/// https://gis.stackexchange.com/questions/20200/how-do-you-compute-the-earths-radius-at-a-given-geodetic-latitude
///

public func crossTrackError(as13: Double, b13: Double, b12: Double, lat: Double) -> (xt: Double, at: Double) {
    let ellipsoid: (a: Double, f: Double) = wgs84
    let a = ellipsoid.a
    let b = a * (1.0 - ellipsoid.f)
    let f = toRadians(lat)
    let R = sqrt(( pow((a*a) * cos(f), 2) + pow(b*b * sin(f), 2)) / (pow(a * cos(f), 2) + pow(b * sin(f),2)))
    let asxt = asin(sin(toRadians(as13)) * sin(toRadians(b13 - b12))) * R
    let asat = acos(cos(toRadians(as13))/cos(asxt/R)) * R
    return (asxt, asat)
}

///
/// Compute a new point point on an ellipsoid give a a starting point, a true azimuth and distance in meteres.
/// The ellipsoid parameters default to the WGS-84 parameters.
///
/// - Parameters:
///   - l1: first point as a tuple with latitude (lat) and longitude (lon) in degrees.
///   - l2: second tuple with azimuth in degrees clockwise from True North, and distance in meters.
///   - a: first ellipsoid parameter in meters (defaults to WGS-84 parameter)
///   - f: second ellipsoid parameter in meters (defaults to WGS-84 parameter)
///
/// - Returns: the second point as a tuple with latitude (lat) and longitude (lon).
///
public func project(_ l1: (lat: Double, lon: Double),
                    _ l2: (azi: Double, dis: Double),
                    ellipsoid: (a: Double, f: Double) = wgs84
) -> (lat: Double, lon: Double) {
    assert(l1.lat >= -90.0 && l1.lat <= 90.0, "l1.lat '\(l1.lat)' outside [-90, 90]")
    assert(l1.lon >= -180.0 && l1.lon <= 180.0, "l1.lon '\(l1.lon)' outside [-180, 180]")
    assert(l2.azi >= 0.0 && l2.azi <= 360, "l2.azi '\(l2.azi))' outside [0, 360]")

    if l2.dis <= 0.0 {
        return l1
    }

    let lat1 = l1.lat
    let lon1 = l1.lon
    let azi1 = l2.azi
    let s12 = l2.dis

    var lat2 = Double.nan
    var lon2 = Double.nan
    var g: geod_geodesic = geod_geodesic()
    geod_init(&g, ellipsoid.a, ellipsoid.f)
    geod_direct(&g, lat1, lon1, azi1, s12, &lat2, &lon2, nil)
    return (lat2, lon2)
}

///
/// Compute a new point point on an ellipsoid give a a starting point, a true azimuth and distance in meteres.
/// The ellipsoid parameters default to the WGS-84 parameters.
///
/// - Parameters:
///   - loc: first point as a CLLocation.
///   - az: azimuth in degrees clockwise from True North
///   - meteres: distance in meters.
///
/// - Returns: the second point as a CLLocation.
///
public func project(loc: CLLocation, az: Double, meters: Double) -> CLLocation {
    let dest = project((lat: loc.coordinate.latitude, lon: loc.coordinate.longitude),
                   (azi: az, dis: meters))
    return CLLocation(latitude: dest.lat, longitude: dest.lon)
}

///
/// Compute a new point point on an ellipsoid give a a starting point, a true azimuth and distance in meteres.
/// The ellipsoid parameters default to the WGS-84 parameters.
///
/// - Parameters:
///   - loc: first point as a CLLocationCoordinate2D.
///   - az: azimuth in degrees clockwise from True North
///   - meteres: distance in meters.
///
/// - Returns: the second point as a CLLocationCoordinate2D.
///
public func project(loc: CLLocationCoordinate2D, az: Double, meters: Double) -> CLLocationCoordinate2D {
    let dest = project((lat: loc.latitude, lon: loc.longitude), (azi: az, dis: meters))
    return CLLocationCoordinate2D(latitude: dest.lat, longitude: dest.lon)
}
