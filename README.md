# geodesic

[![Swift Version](https://img.shields.io/badge/swift-5.9-blue.svg)](https://swift.org)
![Platform](https://img.shields.io/badge/platform-macOS|linux--64-lightgray.svg)
![Build](https://github.com/dastrobu/geodesic/actions/workflows/ci.yaml/badge.svg)
[![GeographicLib Version](https://img.shields.io/badge/GeographicLib-2.1-blue.svg)](https://github.com/geographiclib/geographiclib-c/releases/tag/v2.1)

Solver for a number of geodesic, UTM and MGRS problems in Swift. Forked from 
[dastrobu/geodesic](https://github.com/dastrobu/geodesic)

The inverse geodesic problem must be solved to compute the distance between two points on an oblate spheroid, or
ellipsoid in general. The generalization to ellipsoids, which are not oblate spheroids is not further considered here,
hence the term ellipsoid will be used synonymous with oblate spheroid.

The distance between two points is also known as the
[Vincenty distance](https://en.wikipedia.org/wiki/Vincenty's_formulae).

Here is an example to compute the distance between two points (the poles in this case) on the
[WGS 84 ellipsoid](https://en.wikipedia.org/wiki/World_Geodetic_System).

    import geodesic
    let d = distance((lat: Double.pi / 2,lon: 0), (lat: -Double.pi / 2, lon: 0))

Here is an example to project a new point from an existing point on the
[WGS 84 ellipsoid](https://en.wikipedia.org/wiki/World_Geodetic_System).

    import geodesic
    let dest = project(CLLocationCoordinate2D(latitude: 40.7128, longitude: -74.0060), az: 180.0, meters: 5000.0)

and that's it so far. I will be adding conversion to and from UTM and MGRS coordinates using
[hobuinc/mgrs](https://github.com/hobuinc/mgrs.git)

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
## Table of Contents

- [Installation](#installation)
  - [Swift Package Manager](#swift-package-manager)
- [Implementation Details](#implementation-details)
- [Convergence and Tolerance](#convergence-and-tolerance)
- [WGS 84 and other Ellipsoids](#wgs-84-and-other-ellipsoids)
- [Known Issues](#known-issues)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Installation

At least `clang-3.6` is required. On linux one might need to install it explicitly. There are no dependencies on macOS.

### Swift Package Manager

```swift
let package = Package(
    dependencies: [
        .package(url: "https://github.com/dastrobu/geodesic.git", from: "1.4.0"),
    ]
)
```

## Implementation Details

This Swift package is a wrapper for the
[C implementation of the geodesic routines in GeographicLib](https://github.com/geographiclib/geographiclib-c).
The goal of this Swift package is to make some algorithms from
GeographicLib available to the Swift world. Alternatively one can employ the
package [vincenty](https://github.com/dastrobu/vincenty)
which is a much simpler solver for the inverse geodesic problem, completely written in Swift. Vincenty's formulae does,
however, have some convergence problems in rare cases and may not give the same accuracy as Karney's algorithm.

## Convergence and Tolerance

The computation does always converge and is said to compute up to machine precision. See documentation
of [GeographicLib](https://geographiclib.sourceforge.io/) for details.

## WGS 84 and other Ellipsoids

By default the
[WGS 84 ellipsoid](https://en.wikipedia.org/wiki/World_Geodetic_System)
is employed, but different parameters can be specified, e.g. for the
[GRS 80 ellipsoid](https://en.wikipedia.org/wiki/GRS_80).

    distance((lat: Double.pi / 2, lon: 0), (lat: -Double.pi / 2, lon: 0), 
             ellipsoid (a: 6378137.0, f: 1/298.257222100882711))

## Known Issues

* Compilation with gcc on Linux does not work. `swift build` fails. No problems with clang on Linux. 
 

