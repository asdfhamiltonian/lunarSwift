//
//  iosluna.swift
//  moon info
//
// Moon 4 arcminute calculation
//Adapted from original code in QBASIC by Keith Burnett
//Acessed from http://www.stargazing.net/kepler/moon2.html
//Further adapted with algorithms/equations from page "How to Compute Planetary Positions"
//by Paul Schlyter http://www.stjarnhimlen.se/ppcomp.html#20

import Foundation

func degrees(rads: Double) -> Double {
    return (rads/M_PI) * 180
}

func radians(degrees: Double) -> Double {
    return (degrees/180) * M_PI
}

func fnday(y: Double, m: Double, d: Double, h:Double) -> Double {
    var a = 367 * y
    a -= 7 * floor((y + floor((m+9) / 12))/4)
    a += floor(275 * m / 9)
    a += d - 730530
    a += h / 24
    return a
}

//Ok this next function may be a problem depending on how swift interprets the copysign function compared to python
func fnipart(x: Double)->Double {
    return copysign(1, x) * floor(abs(x))
}

func fnrange(angle: Double) -> Double {
    var newangle = angle % (2*M_PI)
    if newangle < 0 {
        newangle = newangle + 2*M_PI
    }
    return newangle
}

func fnatan2(y: Double, x: Double) -> Double {
    var a = atan2(y, x)
    if a < 0 {
        a = a + 2 * M_PI
    }
    return a
}

func timeadjust(x: Double, utc: Double) -> Double {
    var xnew = x + utc
    while (xnew < 0) {
        xnew = xnew + 24.0
    }
    if (xnew > 24) {
        xnew = xnew - 24
    }
    return xnew
}

func decimaltominsecs(x: Double) -> [Double] {
    var xnew = x%1
    let mins = floor(xnew * 60)
    xnew = (xnew * 60) % 1
    let secs = floor(xnew * 60)
    return [mins, secs]
}

func calculate(y: Double, m: Double, d: Double, h:Double, mins: Double, second: Double, llat: Double, llon: Double) -> [Double]{
    let h = h + mins / 60 + second / 3600
    let d = fnday(y,m: m,d: d,h: h)
    
    // Moon orbital elements
    let Nm = fnrange(radians(125.1228 - 0.0529538083 * d))
    let im = radians(5.1454)
    let wm = fnrange(radians(318.0634 + 0.1643573223 * d))
    let am = 60.2666 //mean radius in earth radii?
    let ecm = 0.0549
    let Mm = fnrange(radians(115.3654 + 13.0649929509 * d))
    
    //Sun elements next:
    //let Ns = 0
    //let isun = 0
    let ws = fnrange(radians(282.9404 + 4.70935E-05 * d))
    //let asun = 1.0 //Astronomical units
    let ecs = 0.016709 - 1.151E-09 * d
    let Ms = fnrange(radians(356.047 + 0.9856002585 * d))
    
    //Sun Position (For Aziumuth Info)
    let Es = Ms + ecs * sin(Ms) * (1.0 + ecs * cos(Ms))
    let xv = cos(Es) - ecs
    let yv = sqrt(1 - ecs * ecs) * sin(Es)
    let vs = atan2(yv, xv)
    let rs = sqrt(xv*xv + yv*yv )
    let lonsun = ws + vs
    let xs = rs * cos(lonsun)
    let ys = rs * sin(lonsun)
    //since the Sun always is in the ecliptic plane, zs is of course zero??
    let xes = xs
    let yes = ys * cos(ecs)
    let zes = ys * sin(ecs)
    var ra = atan2(yes, xes)
    var dec = atan2(zes, sqrt(xes*xes + yes*yes))
    
    
    //Position of Moon
    let Em = Mm + ecm * sin(Mm) * (1 + ecm * cos(Mm))
    let xvm = am * (cos(Em) - ecm)
    let yvm = am * (sqrt(1 - ecm * ecm) * sin(Em))
    let vm = fnatan2(yvm, x: xvm)
    var rm = sqrt(xvm * xvm + yvm * yvm)
    let xh = rm * (cos(Nm) * cos(vm + wm) - sin(Nm) * sin(vm + wm) * cos(im))
    let yh = rm * (sin(Nm) * cos(vm + wm) + cos(Nm) * sin(vm + wm) * cos(im))
    let zh = rm * (sin(vm + wm) * sin(im))
    
    //So looks like I recycled the variables named xv and yv in the original python code? Woops!
    
    //moon's geometric longitude and latitude
    var lon = fnatan2(yh, x: xh)
    var lat = fnatan2(zh, x: sqrt(xh * xh + yh * yh))
    if (degrees(lat) > 90) {
        lat = lat - 2*M_PI
    }
    
    //perturbations
    //first need to calculate the args below, which are in radians
    //MS, Mm - Mean anomaly of the sun and moon
    //Nm - longitude of the Moon's node
    //ws, wm - Argument of the perihelion for the Sun and the Moon
    let Ls = Ms + ws // Mean longitude of the Sun (Ns = 0)
    let Lm = Mm + wm + Nm //Mean longitude of the Moon
    let dm = Lm - Ls //Mean elongation of the Moon
    let F = Lm - Nm //Argument of latitude for the Moon
    
    //Then add the following terms to the longitude
    // Note amplitudes are in degrees, convert at end
    var dlon = -1.274 * sin(Mm - 2 * dm) //the Evection?
    dlon = dlon + 0.658 * sin(2 * dm) // the Variation?
    dlon = dlon - 0.186 * sin(Ms) // the Yearly Equation?
    dlon = dlon - 0.059 * sin(2 * Mm - 2 * dm)
    dlon = dlon - 0.057 * sin(Mm - 2*dm + Ms)
    dlon = dlon + 0.053 * sin(Mm + 2 * dm)
    dlon = dlon + 0.046 * sin(2 * dm - Ms)
    dlon = dlon + 0.041 * sin(Mm - Ms)
    dlon = dlon - 0.035 * sin(dm)   //The Parallactic Equation??
    dlon = dlon - 0.031 * sin(Mm + Ms)
    dlon = dlon - 0.015 * sin(2 * F - 2 * dm)
    dlon = dlon + 0.011 * sin(Mm - 4 * dm)
    lon = radians(dlon) + lon
    
    //latitude terms
    var dlat = -0.173 * sin(F - 2 * dm)
    dlat = dlat - 0.055 * sin(Mm - F - 2 * dm)
    dlat = dlat - 0.046 * sin(Mm + F - 2 * dm)
    dlat = dlat + 0.033 * sin(F + 2 * dm)
    dlat = dlat + 0.017 * sin(2 * Mm + F)
    lat = radians(dlat) + lat
    
    //distance terms earth radii
    rm = rm - 0.58 * cos(Mm - 2 * dm)
    rm = rm - 0.46 * cos(2 * dm)
    
    //Next find the cartesian coordinates of the geocentric lunar position
    let xg = rm * cos(lon) * cos(lat)
    let yg = rm * sin(lon) * cos(lat)
    let zg = rm * sin(lat)
    //rotate the equatorial coords
    //obliquity of ecliptic of date
    let ecl = radians(23.4393 - 3.563E-07 * d)
    let xe = xg
    let ye = yg * cos(ecl) - zg * sin(ecl)
    let ze = yg * sin(ecl) + zg * cos(ecl)
    
    //geocentric RA and Dec
    ra = fnatan2(ye, x: xe)
    dec = fnatan2(ze, x: sqrt(xe*xe + ye*ye))
    //var lstring = "Lat: " + degrees(lat) commenting these out for now
    //var lonstring = "Lon: " + degrees(lon)
    
    //next sidereal time:
    //print(Ls)
    //print(vs+ws)
    //localLong = float(input("longitude? "))
    //localLat = float(input("latitude? "))
    
    let lat1 = radians(llat)
    let gmsto = 12 * (Ls + M_PI)/(M_PI)
    let gmst = gmsto + h
    var lst = gmst + llon/15
    lst = lst%24
    var ha = lst - degrees(ra)/15
    while ((ha > 24) || (ha < -24)){ //This part could get weird with swift not having elif
        if (ha > 12) {
            ha = ha - 24 }
        if (ha < -12) {
            ha = ha + 24 }
    }
    ha = radians(ha*15)
    let x = cos(ha) * cos(dec)
    let y = sin(ha) * cos(dec)
    let z = sin(dec)
    
    let xhor = x * sin(lat1) - z * cos(lat1)
    let yhor = y
    let zhor = x * cos(lat1) + z * sin(lat1)
    
    let az = atan2( yhor, xhor) + M_PI
    let alt = asin(zhor)
    
    //#print("Azimuth is {} degrees".format(str(degrees(az))))
    //#print("Alt is {} degrees".format(str(degrees(alt))))
    
    //#now to correct for topocentricity?
    //#see http://www.stjarnhimlen.se/comp/ppcomp.html#12b for more info about this
    let mpar = asin(1 / rm)
    let msd = asin(0.2725 / rm)
    let alttopoc = alt - mpar * cos(alt)

    //Next calculating "geocentric lattitude" accounting for flattening of earth:
    let gclat = lat1 - radians(0.1924) * sin(2 * lat1)
    let rho = 0.99833 + 0.00167 * cos(2 * lat1)
    let g = atan( tan(gclat) / cos(ha)) //g is the auxiliary angle?
    let topRA = ra - mpar * rho * cos(gclat) * sin(ha) / cos(dec)
    var topDecl = 0.0
    if (abs(abs(dec) - M_PI/2) < 0.001) || (abs(gclat) < 0.001) {
        topDecl = dec - mpar * rho * sin(-dec) * cos(ha)} else {
    topDecl = dec - mpar * rho * sin(gclat) * sin(g - dec) / sin(g)
    }
    
    //moon phase calcs
    let slon = lonsun
    let elong = acos( cos(slon - lon) * cos(lat) )
    let FV = M_PI - elong
    let phase = (1 + cos(FV))/2
    // phasestr = "The current phase of the moon is {}.".format(str(phase))
    
    /*self.d = d
    self.ra = ra
    self.dec = dec
    self.lstring = lstring
    self.lonstring = lonstring
    self.az = az
    self.alttopoc = alttopoc
    self.phasestr = phasestr
    self.mra = topRA
    self.mdecl = topDecl
    self.mpar = mpar
    self.msd = msd
    self.gmsto = gmsto*/
    
    return [d, ra, dec, lat, lon, az, alttopoc, phase, topRA, topDecl, mpar, msd, gmsto]
}

func now(llat: Double, llon: Double, utcDis: Double) -> [Double] {
    //NSDate code? from https://ios8programminginswift.wordpress.com/2014/08/16/get-current-date-and-time-quickcode/
    var date = NSDate()
    date = date.dateByAddingTimeInterval(-3600 * utcDis)
    let calendar = NSCalendar.currentCalendar()
    let components = calendar.components([.Hour, .Minute, .Second, .Month, .Year, .Day], fromDate: date)
    let h = Double(components.hour)
    let mins = Double(components.minute)
    let second = Double(components.second)
    let m = Double(components.month)
    let y = Double(components.year)
    let d = Double(components.day)
    
    let nowarray = calculate(y, m: m, d: d, h: h, mins: mins, second: second, llat: llat, llon: llon)
    return nowarray
}

func riseset(llat: Double, llon: Double, day:Double, utcDis: Double) -> [Double] {
    
    var date = NSDate()
    date = date.dateByAddingTimeInterval(-3600 * utcDis)
    date = date.dateByAddingTimeInterval(3600*24*day)
    var calendar = NSCalendar.currentCalendar()
    var components = calendar.components([.Hour, .Minute, .Second, .Month, .Year, .Day], fromDate: date)
    var h = Double(components.hour)
    let mins = Double(components.minute)
    let second = Double(components.second)
    var m = Double(components.month)
    var y = Double(components.year)
    var d = Double(components.day)
    var array = calculate(y, m: m, d: d, h: h, mins: mins, second: second, llat: llat, llon: llon)
    // return [d, ra, dec, lat, lon, az, alttopoc, phase, topRA, topDecl, mpar, msd, gmsto]
    
    d = array[0]
    /* var ra = array[1]
    var dec = array[2]
    var lat = array[3]
    var lon = array[4]
    var az = array[5]
    var alttopoc = array[6]
    var phase = array[7] */
    var topRA = array[8]
    var topDecl = array[9]
    //var mpar = array[10]
    var msd = array[11]
    var gmsto = array[12]
    
    var hmm = radians(-0.583) - msd
    var mlha = (sin(hmm) - sin(radians(llat)) * sin(topDecl))/(cos(radians(llat)) * cos(topDecl))
    mlha = acos(mlha) * 12 * 15.0 / (15.04107 * M_PI)
    var utcmoon = (degrees(topRA) - gmsto * 15 - llon) / 15
    var a = mlha
    var b = 0.0
    var c = 1.0
    var ff = 0.0
    
    while (abs(c) > 0.000001) {
        date = NSDate()
        date = date.dateByAddingTimeInterval(-3600 * utcDis)
        date = date.dateByAddingTimeInterval(3600*24*day)
        calendar = NSCalendar.currentCalendar()
        components = calendar.components([.Hour, .Minute, .Second, .Month, .Year, .Day], fromDate: date)
        y = Double(components.year)
        m = Double(components.month)
        d = Double(components.day)
        h = utcmoon + a
        array = calculate(y, m: m, d: d, h: h, mins: 0, second: 0, llat: llat, llon: llon)
        
        d = array[0]
        /* ra = array[1]
        dec = array[2]
        lat = array[3]
        lon = array[4]
        az = array[5]
        alttopoc = array[6]
        phase = array[7] */
        topRA = array[8]
        topDecl = array[9]
        //mpar = array[10]
        msd = array[11]
        gmsto = array[12]
            
        hmm = radians(-0.583) - msd
        mlha = (sin(hmm) - sin(radians(llat)) * sin(topDecl))/(cos(radians(llat)) * cos(topDecl))
        mlha = acos(mlha) * 12 * 15.0 / (15.04107 * M_PI)
        utcmoon = (degrees(topRA) - gmsto * 15 - llon) / 15
        b = mlha
        c = abs(b - a)
        a = mlha
        ff += 1
        if (ff >= 100){
            break;
        }
    }
    
    var moonset = mlha + utcmoon
    c = 1
    
    while (abs(c) > 0.000001) {
        date = NSDate()
        date = date.dateByAddingTimeInterval(-3600 * utcDis)
        date = date.dateByAddingTimeInterval(3600*24*day)
        calendar = NSCalendar.currentCalendar()
        components = calendar.components([.Hour, .Minute, .Second, .Month, .Year, .Day], fromDate: date)
        y = Double(components.year)
        m = Double(components.month)
        d = Double(components.day)
        h = utcmoon - a
        array = calculate(y, m: m, d: d, h: h, mins: 0, second: 0, llat: llat, llon: llon)
        
        d = array[0]
        /* ra = array[1]
        dec = array[2]
        lat = array[3]
        lon = array[4]
        az = array[5]
        alttopoc = array[6]
        phase = array[7] */
        topRA = array[8]
        topDecl = array[9]
        //mpar = array[10]
        msd = array[11]
        gmsto = array[12]
        
        hmm = radians(-0.583) - msd
        mlha = (sin(hmm) - sin(radians(llat)) * sin(topDecl))/(cos(radians(llat)) * cos(topDecl))
        mlha = acos(mlha) * 12 * 15.0 / (15.04107 * M_PI)
        utcmoon = (degrees(topRA) - gmsto * 15 - llon) / 15
        b = mlha
        c = b - a
        a = mlha
        ff += 1
        if (ff >= 200){
            break;
        }
    }
    
    var moonrise = utcmoon - mlha
    moonrise = timeadjust(moonrise, utc: utcDis)
    var moonrise2 = decimaltominsecs(moonrise)
    moonset = timeadjust(moonset, utc: utcDis)
    var moonset2 = decimaltominsecs(moonset)
    
    return [floor(moonrise), moonrise2[0], moonrise2[1], floor(moonset), moonset2[0], moonset2[1]]
}

func isUp(llat: Double, llon: Double, utcDis: Double) -> Bool {
    //reports is moon is currently up at a given location
    var nowArray = now(llat, llon: llon, utcDis: utcDis)
    //copying some of the calculations from the riseset function
    //Uses dummy value of "0" for UTC dis. This is important for displaying moon
    //parameters in local time but doesn't affect the calculation.
    let alttopoc = nowArray[6]
    let msd = nowArray[11]
    let hmm = radians(-0.583) - msd
    if alttopoc > hmm {
        return true
    } else {
        return false
    }
}

func upAt(llat: Double, llon: Double, y: Double, m: Double, d: Double, h:Double, mins: Double, second: Double) -> Bool {
    //reports if moon is up at a given location at a specific UTC time
    var array = calculate(y, m: m, d: d, h: h, mins: mins, second: second, llat: llat, llon: llon)
    //copying some of the calculations from the riseset function
    //Uses dummy value of "0" for UTC dis. This is important for displaying moon
    //parameters in local time but doesn't affect the calculation.
    let alttopoc = array[6]
    let msd = array[11]
    let hmm = radians(-0.583) - msd
    if alttopoc > hmm {
        return true
    } else {
        return false
    }
}

func dateNow () -> [Double] {
    //returns current local time
    let date = NSDate()
    let calendar = NSCalendar.currentCalendar()
    let components = calendar.components([.Hour, .Minute, .Second, .Month, .Year, .Day], fromDate: date)
    let h = Double(components.hour)
    let mins = Double(components.minute)
    let second = Double(components.second)
    let m = Double(components.month)
    let y = Double(components.year)
    let d = Double(components.day)
    return [y, m, d, h, mins, second]
}

func utcNow (utcDis: Double) -> [Double] {
    var date = NSDate()
    date = date.dateByAddingTimeInterval(-3600 * utcDis)
    let calendar = NSCalendar.currentCalendar()
    let components = calendar.components([.Hour, .Minute, .Second, .Month, .Year, .Day], fromDate: date)
    let h = Double(components.hour)
    let mins = Double(components.minute)
    let second = Double(components.second)
    let m = Double(components.month)
    let y = Double(components.year)
    let d = Double(components.day)
    return [y, m, d, h, mins, second]
}

func advanceDate(y: Double, m: Double, d: Double, h:Double, mins: Double, second: Double, addDays: Double, utcDis: Double) -> [Double] {
    let asdfDate = NSDateComponents()
    asdfDate.year = Int(y)
    asdfDate.month = Int(m)
    asdfDate.day = Int(d)
    asdfDate.hour = Int(h)
    asdfDate.minute = Int(mins)
    asdfDate.second = Int(second)
    let greg = NSCalendar(calendarIdentifier: NSCalendarIdentifierGregorian)
    var returnDate = greg?.dateFromComponents(asdfDate)!
    returnDate = returnDate?.dateByAddingTimeInterval(-3600*utcDis)
    returnDate = returnDate?.dateByAddingTimeInterval(3600*24*addDays)
    let calendar = NSCalendar.currentCalendar()
    let components = calendar.components([.Hour, .Minute, .Second, .Month, .Year, .Day], fromDate: returnDate!)
    let h = Double(components.hour)
    let mins = Double(components.minute)
    let second = Double(components.second)
    let m = Double(components.month)
    let y = Double(components.year)
    let d = Double(components.day)
    return [y, m, d, h, mins, second]
}

func isNextRiseSet(llat: Double, llon: Double, utcDis: Double, y: Double, m: Double, d: Double, h:Double, mins: Double, second: Double) -> Bool {
    //Reports if a given time is the next rise/set time for the current time
    var y = y
    var m = m
    var d = d
    var h = h
    var mins = mins
    var second = second
    var currentTimeArray = dateNow()
    var currentHour = currentTimeArray[3]
    let currentMins = currentTimeArray[4]
    let currentSecond = currentTimeArray[5]
    let riseSetHour = h + mins/60 + second/3600
    currentHour = currentHour + currentMins/60 + currentSecond/3600
    if riseSetHour < currentHour {
        var riseSetArray = advanceDate(y, m: m, d: d, h: h, mins: mins, second: second, addDays: 1, utcDis: 0)
        y = riseSetArray[0]
        m = riseSetArray[1]
        d = riseSetArray[2]
        h = riseSetArray[3]
        mins = riseSetArray[4]
        second = riseSetArray[5]
    }
    //convert to utc
    var utcriseSetArray = advanceDate(y, m: m, d: d, h: h, mins: mins, second: second, addDays: 0, utcDis: utcDis)
    y = utcriseSetArray[0]
    m = utcriseSetArray[1]
    d = utcriseSetArray[2]
    h = utcriseSetArray[3]
    mins = utcriseSetArray[4]
    second = utcriseSetArray[5]
    
    var array = calculate(y, m: m, d: d, h: h, mins: mins, second: second, llat: llat, llon: llon)
    //copying some of the calculations from the riseset function
    //Uses dummy value of "0" for UTC dis. This is important for displaying moon
    //parameters in local time but doesn't affect the calculation.
    let alttopoc = array[6]
    let msd = array[11]
    let hmm = radians(-0.583) - msd
    //figures out if moon alt is within 1 degree of the rise/set altitude
    //(the calculations for rise/set aren't exact so there needs to be
    //some wiggle room
    if (alttopoc < hmm + radians(0.2)) && (alttopoc > hmm - radians(0.2)) {
        return true
    } else {
        return false
    }
}

func moonTest(moon: [Double], riseOrSet: Int, llat: Double, llon: Double, utcDis: Double, y: Double, m: Double, d: Double) -> Bool {
    //tests whether a particular moonrise or set time is correct
    if riseOrSet == 0 {
        //testing for moonrise time
        let hour = moon[0]
        let mins = moon[1]
        let second = moon[2]
        let a = isNextRiseSet(llat, llon: llon, utcDis: utcDis, y: y, m: m, d: d, h: hour, mins: mins, second: second)
        return a
    } else {
        //testing for moonset time
        let hour = moon[3]
        let mins = moon[4]
        let second = moon[5]
        let a = isNextRiseSet(llat, llon: llon, utcDis: utcDis, y: y, m: m, d: d, h: hour, mins: mins, second: second)
        return a
    }
}

func dplaces(number: Double, sigfigs: Double) -> Double {
    let part1 = floor(number)
    let part2 = Double(round((number%1) * pow(10, sigfigs))) * pow(10, -sigfigs)
    return part1 + part2
}

func nextRiseSetTime(llat: Double, llon: Double, utcDis: Double) -> [Double] {
    var currentTimeArray = dateNow()
    let currentYear = currentTimeArray[0]
    let currentMonth = currentTimeArray[1]
    let currentDay = currentTimeArray[2]
    //let currentHour = currentTimeArray[3]
    //let currentMins = currentTimeArray[4]
    //let currentSecond = currentTimeArray[5]
    
    //rise set arrays:
    var moon1 = riseset(llat, llon: llon, day: 1, utcDis: utcDis)
    var moon2 = riseset(llat, llon: llon, day: 2, utcDis: utcDis)
    var moon3 = riseset(llat, llon: llon, day: 3, utcDis: utcDis)
    
    if isUp(llat, llon: llon, utcDis: utcDis) {
        if moonTest(moon1, riseOrSet: 1, llat: llat, llon: llon, utcDis: utcDis, y: currentYear, m: currentMonth, d: currentDay) {
            return [moon1[3], moon1[4], moon1[5]]
        }
        if moonTest(moon2, riseOrSet: 1, llat: llat, llon: llon, utcDis: utcDis, y: currentYear, m: currentMonth, d: currentDay) {
            return [moon2[3], moon2[4], moon2[5]]
        }
        if moonTest(moon3, riseOrSet: 1, llat: llat, llon: llon, utcDis: utcDis, y: currentYear, m: currentMonth, d: currentDay) {
            return [moon3[3], moon3[4], moon3[5]]
        }
    } else {
        if moonTest(moon1, riseOrSet: 0, llat: llat, llon: llon, utcDis: utcDis, y: currentYear, m: currentMonth, d: currentDay) {
            return [moon1[0], moon1[1], moon1[2]]
        }
        if moonTest(moon2, riseOrSet: 0, llat: llat, llon: llon, utcDis: utcDis, y: currentYear, m: currentMonth, d: currentDay) {
            return [moon2[0], moon2[1], moon2[2]]
        }
        if moonTest(moon3, riseOrSet: 0, llat: llat, llon: llon, utcDis: utcDis, y: currentYear, m: currentMonth, d: currentDay) {
            return [moon3[0], moon3[1], moon3[2]]
        }
    }
    
    return[2.0, 3.0, 4.0]
}
