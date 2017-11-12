//iosluna.swift
//A swift 4.0 class to calculate where the moon is in the sky, as well as moon rise and set times
//Adapted from original code in QBASIC by Keith Burnett
//Acessed from http://www.stargazing.net/kepler/moon2.html
//Further adapted with algorithms/equations from page "How to Compute Planetary Positions" by Paul Schlyter http://www.stjarnhimlen.se/ppcomp.html#20

import Foundation

func degrees(rads: Double) -> Double {
    return (rads/Double.pi * 180.0)
}

func radians(degrees: Double) -> Double {
    return (degrees/180.0 * Double.pi)
}

func fnday(y: Double, m:Double, d: Double, h: Double) -> Double {
    var a: Double
    a = 367 * y - floor(7 * (y + floor((m + 9) / 12))/4) + floor(275 * m / 9) + d - 730530 + h / 24
    return a
}

func fnipart(x: Double) -> Double {
    return copysign(1.0,x) * floor(abs(x))
}

func fnrange(angle: Double) -> Double {
    var newangle: Double
    newangle = angle.truncatingRemainder(dividingBy: 2*Double.pi)
    if(newangle < 0){
        newangle = newangle + 2*Double.pi
    }
    return newangle
}

//these check out so far compared to python code

func fnatan2(y: Double, x: Double) -> Double {
    var a: Double
    a = atan2(y, x)
    if(a < 0){
        a = a + 2 * Double.pi
    }
    return a
}

//This one works too!

func timeadjust(x: Double, utc: Double) -> Double {
    var a: Double
    a = x + utc
    while(a < 0){
        a = a + 24
    }
    if(a > 24){
        a = a - 24
    }
    return a
}

func decimaltominsecs(x: Double) -> [Double] {
    var a: Double
    var mins: Double
    var secs: Double
    a = x.truncatingRemainder(dividingBy: 1.0)
    mins = floor(a * 60)
    a = (a * 60).truncatingRemainder(dividingBy: 1.0)
    secs = floor(a * 60)
    return [mins, secs]
}

func longituderange(phi: Double) -> Double {
    return (phi + 180.0).truncatingRemainder(dividingBy: 360.0) - 180
}

func latituderange(theta: Double) -> Double {
    return abs(((theta - 90.0).truncatingRemainder(dividingBy: 360.0)) - 180) - 90
}

class Luna{
    var localLat = 45.0
    var localLong = -122.0
    
    func calculate(y: Double, m: Double, d: Double, h: Double, mins: Double, second: Double) -> [Double] {
        var h2:Double
        h2 = h + mins / 60.0 + second / 3600.0
        var d2:Double
        d2 = fnday(y: y, m: m, d: d, h: h2)
        
        //These are the moon's orbital elements:
        
        var Nm:Double, im:Double, wm:Double, am:Double, ecm:Double, Mm:Double
        Nm = fnrange(angle: radians(degrees: (125.1228 - 0.0529538083 * d2)))
        im = radians(degrees: 5.1454)
        wm = fnrange(angle: radians(degrees: (318.0634 + 0.1643573223 * d2)))
        am = 60.2666 //apparently the mean radius of the moon in earth radii?
        ecm = 0.0549
        Mm = fnrange(angle: radians(degrees: (115.3654 + 13.0649929509 * d2)))
        
        //Sun orbital elements
        var Ns:Double, isun:Double, ws:Double, asun:Double, ecs:Double, Ms:Double
        Ns = 0.0
        isun = 0.0
        ws = fnrange(angle: radians(degrees: (282.9404 + 4.70935E-05 * d2)))
        asun = 1.0 //Astronomical units
        ecs = 0.016709 - 1.151E-09 * d2
        Ms = fnrange(angle: radians(degrees: (356.047 + 0.9856002585 * d2)))
        
        //Sun position
        var Es:Double, xv:Double, yv:Double, vs:Double, rs:Double, lonsun:Double, xs:Double, ys:Double, xes:Double, yes:Double, zes:Double, ra:Double, dec:Double
        Es = Ms + ecs * sin(Ms) * (1.0 + ecs * cos(Ms))
        xv = cos(Es) - ecs
        yv = sqrt(1 - ecs * ecs) * sin(Es)
        vs = atan2(yv, xv)
        rs = sqrt(xv * xv + yv * yv)
        lonsun = ws + vs
        xs = rs * cos(lonsun)
        ys = rs * sin(lonsun)
        //Since the sun is always in the ecliptic plane, zs is 0.
        xes = xs
        yes = ys * cos(ecs)
        zes = ys * sin(ecs)
        ra = atan2(yes, xes)
        dec = atan2(zes, sqrt(xes * xes + yes * yes))
        
        //Moon position
        var xh:Double, yh:Double, zh:Double, Em:Double, vm:Double, rm:Double
        Em = Mm + ecm * sin(Mm) * (1 + ecm * cos(Mm))
        xv = am * (cos(Em) - ecm)
        yv = am * (sqrt(1 - ecm * ecm) * sin(Em))
        vm = fnatan2(y: yv, x: xv)
        rm = sqrt(xv * xv + yv * yv)
        xh = rm * (cos(Nm) * cos(vm + wm) - sin(Nm) * sin(vm + wm) * cos(im))
        yh = rm * (sin(Nm) * cos(vm + wm) + cos(Nm) * sin(vm + wm) * cos(im))
        zh = rm * (sin(vm + wm) * sin(im))
        
        //Moon geometric longitude and latitude
        var lon:Double, lat:Double
        lon = fnatan2(y: yh, x: xh)
        lat = fnatan2(y: zh, x: (sqrt(xh * xh + yh * yh)))
        
        //Perturbations are calculated next in radians.
        //Ms, Mm - Mean anomaly of the sun and moon
        //Nm - longitutde of the Moon's node
        //ws, wm - Argument of the perihelion of the Sun and Moon
        
        var Ls:Double, Lm:Double, dm:Double, F:Double
        Ls = Ms + ws
        Lm = Mm + wm + Nm
        dm = Lm - Ls
        F = Lm - Nm
        
        //Then add the following terms to the longitude
        //Note the amplitudes are in degrees, convert at end
        var dlon:Double
        dlon = -1.274 * sin(Mm - 2.0 * dm) //The evection
        dlon = dlon + 0.658 * sin(2.0 * dm) //The variation
        dlon = dlon - 0.186 * sin(Ms) //The yearly equation
        dlon = dlon - 0.059 * sin(2.0 * Mm - 2.0 * dm)
        dlon = dlon - 0.057 * sin(Mm - 2.0 * dm + Ms)
        dlon = dlon + 0.053 * sin(Mm + 2.0 * dm)
        dlon = dlon + 0.046 * sin(2.0 * dm - Ms)
        dlon = dlon + 0.041 * sin(Mm - Ms)
        dlon = dlon - 0.035 * sin(dm)
        dlon = dlon - 0.031 * sin(Mm + Ms)
        dlon = dlon - 0.015 * sin(2.0 * F - 2.0 * dm)
        dlon = dlon + 0.011 * sin(Mm - 4.0 * dm)
        lon = radians(degrees: dlon) + lon
        
        //Then these terms adjust the latitude
        var dlat:Double
        dlat = -0.173 * sin(F - 2 * dm)
        dlat = dlat - 0.055 * sin(Mm - F - 2 * dm)
        dlat = dlat - 0.046 * sin(Mm + F - 2 * dm)
        dlat = dlat + 0.033 * sin(F + 2 * dm)
        dlat = dlat + 0.017 * sin(2 * Mm + F)
        lat = radians(degrees: dlat) + lat
        
        //Need to correct so latitude is between pi/2 and -pi/2
        //(+90 or -90 degrees) The above formula had been giving
        //latitudes like 355.679..."
        //There's some math involved that doesn't make much sense
        //unless you plot desired output vs input
        
        while(lat < 0){
            lat = lat + 2 * Double.pi
        }
        if((lat > Double.pi/2.0) && (lat <= 3.0 * Double.pi/2.0)){
            lat = Double.pi - lat
        } else if((lat>3.0 * Double.pi / 2.0) && (lat <= 5.0 * Double.pi/2.0)){
            lat = lat - 2.0 * Double.pi
        }
        
        //Distance terms earth radii
        rm = rm - 0.58 * cos(Mm - 2.0 * dm)
        rm = rm - 0.46 * cos(2.0 * dm)
        
        //Next to find the cartesian coordinates of the geocentric lunar position
        var xg:Double, yg:Double, zg:Double
        xg = rm * cos(lon) * cos(lat)
        yg = rm * sin(lon) * cos(lat)
        zg = rm * sin(lat)
        
        //now rotate the equatorial coordinates
        //obliquity of ecliptic date
        var ecl:Double, xe:Double, ye:Double, ze:Double
        ecl = radians(degrees: (23.4393 - 3.563E-07 * d2))
        xe = xg
        ye = yg * cos(ecl) - zg * sin(ecl)
        ze = yg * sin(ecl) + zg * cos(ecl)
        
        //geocentric right ascension and declension
        var latdegrees:Double, londegrees:Double
        ra = fnatan2(y: ye, x: xe)
        dec = fnatan2(y: ze, x: sqrt(xe * xe + ye * ye))
        latdegrees = latituderange(theta: degrees(rads: lat))
        londegrees = longituderange(phi: degrees(rads: lon))
        
        var lat1:Double, gmsto:Double, gmst:Double, lst:Double, ha:Double, x1:Double, y1:Double, z1:Double
        lat1 = radians(degrees: localLat)
        gmsto = 12.0 * (Ls + Double.pi) / Double.pi
        gmst = gmsto + h2
        lst = gmst + localLong/15.0
        lst = lst.truncatingRemainder(dividingBy: 24.0)
        ha = lst - degrees(rads: ra) / 15.0
        while((ha > 24.0) || (ha < -24.0)){
            if(ha > 12.0){
                ha = ha - 24.0
            } else if(ha < -12.0){
                ha = ha + 24.0
            }
        }
        ha = radians(degrees: (ha * 15.0))
        x1 = cos(ha) * cos(dec)
        y1 = sin(ha) * cos(dec)
        z1 = sin(dec)
        
        var xhor:Double, yhor:Double, zhor:Double
        xhor = x1 * sin(lat1) - z1 * cos(lat1)
        yhor = y1
        zhor = x1 * cos(lat1) + z1 * sin(lat1)
        
        var az:Double, azdegrees:Double, alt:Double
        az = atan2(yhor, xhor) + Double.pi
        azdegrees = degrees(rads: az)
        alt = asin(zhor)
        
        //Correct for topocentricity
        //see www.stjarnhimlen.se/comp/ppcomp.html#12b for more info
        var mpar:Double, msd:Double, alttopoc:Double, altdegrees:Double
        mpar = asin(1.0 / rm)
        msd = asin(0.2725 / rm)
        alttopoc = alt - mpar * cos(alt)
        altdegrees = degrees(rads: alttopoc)
        
        //Calculate geocentric latitude accounting for flattening of the earth
        var gclat:Double, rho:Double, g:Double, topRA:Double, topDecl:Double
        gclat = lat1 - radians(degrees: 0.1924) * sin(2 * lat1)
        rho = 0.99833 + 0.00167 * cos(2 * lat1)
        g = atan(tan(gclat) / cos(ha)) //g is the auxiliary angle
        topRA = ra - mpar * rho * cos(gclat) * sin(ha) / cos(dec)
        topDecl = 0.0
        if((abs(abs(dec) - Double.pi/2.0) < 0.001) || (abs(gclat) < 0.001)){
            topDecl = dec - mpar * rho * sin(-dec) * cos(ha)
        } else {
            topDecl = dec - mpar * rho * sin(gclat) * sin(g - dec) / sin(g)
        }
        
        //Moon phase
        var slon:Double, elong:Double, FV:Double, phase:Double
        slon = lonsun
        elong = acos(cos(slon - lon) * cos(lat))
        FV = Double.pi - elong
        phase = (1.0 + cos(FV))/2.0
        
        return  [azdegrees, altdegrees, ra, dec, phase, topRA, topDecl, mpar, msd, londegrees, latdegrees, gmsto, d]
    }
    func now() -> [Double]{
        var y:Double, m:Double, d:Double, h:Double, mins:Double, second:Double
        var current = Date()
        var formatter = DateFormatter()
        formatter.timeZone = TimeZone(abbreviation: "GMT")
        formatter.dateFormat = "YYYY"
        y = Double(formatter.string(from: current))!
        formatter.dateFormat = "MM"
        m = Double(formatter.string(from: current))!
        formatter.dateFormat = "dd"
        d = Double(formatter.string(from: current))!
        formatter.dateFormat = "hh"
        h = Double(formatter.string(from: current))!
        formatter.dateFormat = "mm"
        mins = Double(formatter.string(from: current))!
        formatter.dateFormat = "ss"
        second = Double(formatter.string(from: current))!
        return calculate(y: y, m: m, d: d, h: h, mins: mins, second: second)
    }
    func riseset(dis: Double, utcdis: Double) -> [Double]{
        var y:Double, m:Double, d:Double, h:Double, mins:Double, second:Double
        var current = Date()
        var formatter = DateFormatter()
        formatter.timeZone = TimeZone(abbreviation: "GMT")
        formatter.dateFormat = "YYYY"
        y = Double(formatter.string(from: current))!
        formatter.dateFormat = "MM"
        m = Double(formatter.string(from: current))!
        formatter.dateFormat = "dd"
        d = Double(formatter.string(from: current))! + dis
        formatter.dateFormat = "hh"
        h = Double(formatter.string(from: current))!
        formatter.dateFormat = "mm"
        mins = Double(formatter.string(from: current))!
        formatter.dateFormat = "ss"
        second = Double(formatter.string(from: current))!
        var moonArray = calculate(y: y, m: m, d: d, h: h, mins: mins, second: second)
        var msd:Double, mra:Double, mdecl:Double, gmsto:Double, hmm:Double, mlha:Double, utcmoon:Double
        msd = moonArray[8]
        mra = moonArray[5]
        mdecl = moonArray[6]
        gmsto = moonArray[11]
        hmm = radians(degrees: -0.583) - msd
        mlha = (sin(hmm) - sin(radians(degrees: localLat)) * sin(mdecl))/(cos(radians(degrees: localLat)) * cos(mdecl))
        mlha = acos(mlha) * 12.0 * 15.0 / (15.04107 * Double.pi)
        utcmoon = (degrees(rads: mra) - gmsto * 15.0 - localLong) / 15.0
        
        //Now for a loop of successive approximations
        var a: Double, c: Double, ff: Double, b: Double
        a = mlha
        c = 1
        ff = 0
        while (abs(c) > 0.000001){
            formatter.timeZone = TimeZone(abbreviation: "GMT")
            formatter.dateFormat = "YYYY"
            y = Double(formatter.string(from: current))!
            formatter.dateFormat = "MM"
            m = Double(formatter.string(from: current))!
            formatter.dateFormat = "dd"
            d = Double(formatter.string(from: current))! + dis
            h = utcmoon + a
            moonArray = calculate(y: y, m: m, d: d, h: h, mins: 0, second: 0)
            msd = moonArray[8]
            mra = moonArray[5]
            mdecl = moonArray[6]
            gmsto = moonArray[11]
            hmm = radians(degrees: -0.583) - msd
            mlha = (sin(hmm) - sin(radians(degrees: localLat)) * sin(mdecl))/(cos(radians(degrees: localLat)) * cos(mdecl))
            mlha = acos(mlha) * 12.0 * 15.0 / (15.04107 * Double.pi)
            utcmoon = (degrees(rads: mra) - gmsto * 15.0 - localLong) / 15.0
            b = mlha
            c = abs(b - a)
            a = mlha
            ff += 1
            if(ff >= 100){
                break
            }
        }
        var moonset = mlha + utcmoon
        c = 1
        while (abs(c) > 0.000001){
            formatter.timeZone = TimeZone(abbreviation: "GMT")
            formatter.dateFormat = "YYYY"
            y = Double(formatter.string(from: current))!
            formatter.dateFormat = "MM"
            m = Double(formatter.string(from: current))!
            formatter.dateFormat = "dd"
            d = Double(formatter.string(from: current))! + dis
            h = utcmoon - a
            moonArray = calculate(y: y, m: m, d: d, h: h, mins: 0, second: 0)
            msd = moonArray[8]
            mra = moonArray[5]
            mdecl = moonArray[6]
            gmsto = moonArray[11]
            hmm = radians(degrees: -0.583) - msd
            mlha = (sin(hmm) - sin(radians(degrees: localLat)) * sin(mdecl))/(cos(radians(degrees: localLat)) * cos(mdecl))
            mlha = acos(mlha) * 12.0 * 15.0 / (15.04107 * Double.pi)
            utcmoon = (degrees(rads: mra) - gmsto * 15.0 - localLong) / 15.0
            b = mlha
            c = abs(b - a)
            a = mlha
            ff += 1
            if(ff >= 200){
                break
            }
        }
        var moonrise = utcmoon - mlha
        moonrise = timeadjust(x: moonrise, utc: utcdis)
        var moonrise2 = decimaltominsecs(x: moonrise)
        moonset = timeadjust(x: moonset, utc: utcdis)
        var moonset2 = decimaltominsecs(x: moonset)
        
        return [floor(moonrise), moonrise2[0], moonrise2[1], floor(moonset), moonset2[0], moonset2[1]]
    }
}

//here's an example of how to use this class

var portland = Luna()
portland.localLat = 45.6
portland.localLong = -122.5
portland.calculate(y:2017, m:11, d:12, h:6, mins:13, second:0)
portland.now()
portland.riseset(dis: 1, utcdis: -8)
portland.riseset(dis: 2, utcdis: -8)
