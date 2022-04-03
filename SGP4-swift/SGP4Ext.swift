//
//  SGP4Ext.swift
//  pxSat3D-SK
//
//  Created by Xav perso on 03/08/2020.
//  Copyright © 2020 P'tit Xav. All rights reserved.
//
/*     ----------------------------------------------------------------
*
*                               adapted from sgp4ext.swift
*
*    this file contains extra routines needed for the main test program for sgp4.
*    these routines are derived from the astro libraries.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                    2007
*                              by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*               7 may 08  david vallado
*                           fix sgn
*    changes :
*               2 apr 07  david vallado
*                           fix jday floor and str lengths
*                           updates for constants
*              14 aug 06  david vallado
*                           original baseline
*       ----------------------------------------------------------------      */
import Foundation


/* -----------------------------------------------------------------------------
*
*                           function newtonnu
*
*  this function solves keplers equation when the true anomaly is known.
*    the mean and eccentric, parabolic, or hyperbolic anomaly is also found.
*    the parabolic limit at 168¯ is arbitrary. the hyperbolic anomaly is also
*    limited. the hyperbolic sine is used because it's not double valued.
*
*  author        : david vallado                  719-573-2600   27 may 2002
*
*  revisions
*    vallado     - fix small                                     24 sep 2002
*
*  inputs          description                    range / units
*    ecc         - eccentricity                   0.0  to
*    nu          - true anomaly                   -2pi to 2pi rad
*
*  outputs       :
*    e0          - eccentric anomaly              0.0  to 2pi rad       153.02 ¯
*    m           - mean anomaly                   0.0  to 2pi rad       151.7425 ¯
*
*  locals        :
*    e1          - eccentric anomaly, next value  rad
*    sine        - sine of e
*    cose        - cosine of e
*    ktr         - index
*
*  coupling      :
*    asinh       - arc hyperbolic sine
*
*  references    :
*    vallado       2007, 85, alg 5
* --------------------------------------------------------------------------- */

func newtonNu(
	_ ecc : Double, _ nu : Double
	) -> (e0 : Double, m : Double)
{
	var e0, m : Double
	var sine, cose : Double

	// ---------------------  implementation   ---------------------
	e0 = infinite
	m = infinite

	// --------------------------- circular ------------------------
	if ( fabs( ecc ) < small  ) {
		m = nu
		e0 = nu
	} else {
		// ---------------------- elliptical -----------------------
		if ( ecc < 1.0-small  ) {
			sine = ( √( 1.0 - ecc*ecc ) * sin(nu) ) / ( 1.0 + ecc*cos(nu) )
			cose = ( ecc + cos(nu) ) / ( 1.0  + ecc*cos(nu) )
			e0  = ArcTan( sine,cose )
			m   = e0 - ecc*sin(e0)
		} else {
			// -------------------- hyperbolic  --------------------
			if ( ecc > 1.0 + small  ) {
				if ((ecc > 1.0 ) && (fabs(nu)+0.00001 < π - acos(1.0 / ecc))) {
					sine = ( √( ecc*ecc-1.0  ) * sin(nu) ) / ( 1.0  + ecc*cos(nu) )
					e0  = asinh( sine )
					m   = ecc*sinh(e0) - e0
				}
			} else {
				// ----------------- parabolic ---------------------
				if fabs(nu) < 168.0 * π/180.0 {
					e0 = tan( nu*0.5  )
					m = e0 + (e0 * e0 * e0)/3.0
				}
			}
		}
	}

	if ( ecc < 1.0  ) {
		m = fmod( m,2.0 * π )
		if ( m < 0.0  ) {
			m = m + 2.0 * π
		}
//		e0 = Modulus( e0,2.0 * π )
		e0 = e0 % (2.0 * π)
	}

	return (e0, m)
}  // end newtonnu


/* -----------------------------------------------------------------------------
*
*                           function rv2coe
*
*  this function finds the classical orbital elements given the geocentric
*    equatorial position and velocity vectors.
*
*  author        : david vallado                  719-573-2600   21 jun 2002
*
*  revisions
*    vallado     - fix special cases                              5 sep 2002
*    vallado     - delete extra check in inclination code        16 oct 2002
*    vallado     - add constant file use                         29 jun 2003
*    vallado     - add mu                                         2 apr 2007
*
*  inputs          description                    range / units
*    r           - ijk position vector            km
*    v           - ijk velocity vector            km / s
*    mu          - gravitational parameter        km3 / s2
*
*  outputs       :
*    p           - semilatus rectum               km
*    a           - semimajor axis                 km
*    ecc         - eccentricity
*    incl        - inclination                    0.0  to pi rad
*    omega       - longitude of ascending node    0.0  to 2pi rad
*    argp        - argument of perigee            0.0  to 2pi rad
*    nu          - true anomaly                   0.0  to 2pi rad
*    m           - mean anomaly                   0.0  to 2pi rad
*    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
*    truelon     - true longitude            (ce) 0.0  to 2pi rad
*    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
*
*  locals        :
*    hbar        - angular momentum h vector      km2 / s
*    ebar        - eccentricity     e vector
*    nbar        - line of nodes    n vector
*    c1          - v**2 - u/r
*    rdotv       - r dot v
*    hk          - hk unit vector
*    sme         - specfic mechanical energy      km2 / s2
*    i           - index
*    e           - eccentric, parabolic,
*                  hyperbolic anomaly             rad
*    temp        - temporary variable
*    typeorbit   - type of orbit                  ee, ei, ce, ci
*
*  coupling      :
*    mag         - magnitude of a vector
*    cross       - cross product of two vectors
*    angle       - find the angle between two vectors
*    newtonnu    - find the mean anomaly
*
*  references    :
*    vallado       2007, 126, alg 9, ex 2-5
* --------------------------------------------------------------------------- */

func rv2coe(_ r : Position, _ v : Velocity, _ mu : Double) ->
	(out_p : Double, out_a : Double, out_ecc : Double, out_incl : Double, out_omega : Double, out_argp : Double,
	out_nu : Double, out_m : Double, out_arglat : Double, out_truelon : Double, out_lonper : Double)
{
	var out_p, out_a, out_ecc, out_incl, out_omega, out_argp,
		out_nu, out_m, out_arglat, out_truelon, out_lonper : Double

	var hbar : Position = kVector3Null
	var nbar : Position = kVector3Null
	var ebar : Position = kVector3Null
	var magr, magv, magn, sme,
		rdotv, temp, c1, hk, magh : Double

	var typeorbit : String

	out_m = undefined

	// -------------------------  implementation   -----------------
	magr = mag( r )
	magv = mag( v )

	// ------------------  find h n and e vectors   ----------------
//	hbar = cross( r, v)
	hbar = r∧v
	magh = mag( hbar )
	if ( magh > small )
	{
		nbar[0] = -hbar[1]
		nbar[1] =  hbar[0]
		nbar[2] =   0.0
		magn = mag( nbar )
		c1 = magv*magv - mu / magr
//		rdotv = dot( r,v )
		rdotv = r⋅v
		for i in 0...2 {
			ebar[i] = (c1*r[i] - rdotv*v[i])/mu
		}
		out_ecc = mag( ebar )

		// ------------  find a e and semi-latus rectum   ----------
		sme = ( magv*magv*0.5  ) - ( mu / magr )
		if ( fabs( sme ) > small ) {
			out_a = -mu  / (2.0 * sme)
		} else {
			out_a = infinite
		}

		out_p = magh*magh/mu

		// -----------------  find inclination   -------------------
		hk = hbar[2]/magh
		out_incl = acos( hk )

		// --------  determine type of orbit for later use  --------
		// ------ elliptical, parabolic, hyperbolic inclined -------
		typeorbit = "ei"
		if (  out_ecc < small ) {
			// ----------------  circular equatorial ---------------
			if  ( out_incl < small) || (fabs( out_incl - π) < small) {
				typeorbit = "ce"
			} else {
				// --------------  circular inclined ---------------
				typeorbit = "ci"
			}
		}
		else
		{
			// - elliptical, parabolic, hyperbolic equatorial --
			if  ( out_incl < small) || (fabs( out_incl - π) < small) {
				typeorbit = "ee"
			}
		}

		// ----------  find longitude of ascending node ------------
		if ( magn > small )
		{
			temp = nbar[0] / magn
			if ( fabs(temp) > 1.0  ) {
				temp = sgn(temp)
			}
			out_omega = acos( temp )
			if ( nbar[1] < 0.0  ) {
				out_omega = Double.twoπ - out_omega
			}
		} else {
			out_omega = undefined
		}

		// ---------------- find argument of perigee ---------------
		if ( typeorbit == "ei" ) {
			out_argp = angle( nbar,ebar)
			if ( ebar[2] < 0.0  ) {
				out_argp = Double.twoπ - out_argp
			}
		} else {
			out_argp = undefined
		}

		// ------------  find true anomaly at epoch    -------------
		if typeorbit.hasPrefix("e") {
			out_nu = angle( ebar,r)
			if ( rdotv < 0.0  ) {
				out_nu = Double.twoπ -  out_nu
			}
		} else {
			out_nu = undefined
		}

		// ----  find argument of latitude - circular inclined -----
		if typeorbit == "ci" {
			out_arglat = angle( nbar,r )
			if ( r[2] < 0.0  ) {
				out_arglat = Double.twoπ - out_arglat
			}
			out_m = out_arglat
		} else {
			out_arglat = undefined
		}

		// -- find longitude of perigee - elliptical equatorial ----
		if  ( out_ecc > small ) && (typeorbit == "ee")
		{
			temp = ebar[0] / out_ecc
			if ( fabs(temp) > 1.0  ) {
				temp = sgn(temp)
			}
			out_lonper = acos( temp )
			if ( ebar[1] < 0.0  ) {
				out_lonper = Double.twoπ -  out_lonper
			}
            if ( out_incl > Double.halfπ ) {
				out_lonper = Double.twoπ -  out_lonper
			}
		} else {
			out_lonper = undefined
		}

		// -------- find true longitude - circular equatorial ------
		if  (( magr>small ) && ( strcmp(typeorbit,"ce") == 0 ))
		{
			temp = r[0]/magr
			if ( fabs(temp) > 1.0  ) {
				temp = sgn(temp)
			}
			out_truelon = acos( temp )
			if ( r[1] < 0.0  ) {
				out_truelon = Double.twoπ -  out_truelon
			}
            if (  out_incl > Double.halfπ ) {
				out_truelon = Double.twoπ -  out_truelon
			}
			out_m =  out_truelon
		}
		else {
			out_truelon = undefined
		}

		// ------------ find mean anomaly for all orbits -----------
		if typeorbit.hasPrefix("e") {
			(_,  out_m) = newtonNu( out_ecc, out_nu)
		}
	}
	else
	{
		out_p       = undefined
		out_a       = undefined
		out_ecc     = undefined
		out_incl    = undefined
		out_omega   = undefined
		out_argp    = undefined
		out_nu      = undefined
		out_m       = undefined
		out_arglat  = undefined
		out_truelon = undefined
		out_lonper  = undefined
	}

	return (out_p, out_a, out_ecc, out_incl, out_omega, out_argp,
		out_nu, out_m, out_arglat, out_truelon, out_lonper)
}  // end rv2coe

/* -----------------------------------------------------------------------------
*
*                           procedure jday
*
*  this procedure finds the julian date given the year, month, day, and time.
*    the julian date is defined by each elapsed day since noon, jan 1, 4713 bc.
*
*  algorithm     : calculate the answer in one step for efficiency
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    year        - year                           1900 .. 2100
*    mon         - month                          1 .. 12
*    day         - day                            1 .. 28,29,30,31
*    hr          - universal time hour            0 .. 23
*    min         - universal time min             0 .. 59
*    sec         - universal time sec             0.0 .. 59.999
*
*  outputs       :
*    jd          - julian date                    days from 4713 bc
*
*  locals        :
*    none.
*
*  coupling      :
*    none.
*
*  references    :
*    vallado       2007, 189, alg 14, ex 3-14
*
* --------------------------------------------------------------------------- */

func jday(_ year : Int, _ mon : Int, _ day : Int, _ hr : Int, _ minute : Int, _ sec : Double) -> Double {
	return 367.0 * Double(year) -
		floor((7 * (Double(year) + floor((Double(mon) + 9) / 12.0))) * 0.25) +
		floor( 275.0 * Double(mon) / 9.0 ) +
		Double(day) + 1721013.5 +
		((sec / 60.0 + Double(minute)) / 60.0 + Double(hr)) / 24.0  // ut in days
}  // end jday

/* -----------------------------------------------------------------------------
*
*                           procedure days2mdhms
*
*  this procedure converts the day of the year, days, to the equivalent month
*    day, hour, minute and second.
*
*  algorithm     : set up array for the number of days per month
*                  find leap year - use 1900 because 2000 is a leap year
*                  loop through a temp value while the value is < the days
*                  perform int conversions to the correct day and month
*                  convert remainder into h m s using type conversions
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    year        - year                           1900 .. 2100
*    days        - julian day of the year         0.0  .. 366.0
*
*  outputs       :
*    mon         - month                          1 .. 12
*    day         - day                            1 .. 28,29,30,31
*    hr          - hour                           0 .. 23
*    min         - minute                         0 .. 59
*    sec         - second                         0.0 .. 59.999
*
*  locals        :
*    dayofyr     - day of year
*    temp        - temporary extended values
*    inttemp     - temporary int value
*    i           - index
*    lmonth[12]  - int array containing the number of days per month
*
*  coupling      :
*    none.
* --------------------------------------------------------------------------- */

func days2mdhmsS(_ year : Int, _ days : Double) -> (out_mon : Int, out_day : Int, out_hr : Int, out_minute : Int, out_sec : Double) {
	var out_mon : Int
	var out_day : Int
	var out_hr : Int
	var out_minute : Int
	var out_sec : Double

	var i, inttemp, dayofyr : Int
	var temp : Double
	var lmonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

	dayofyr = (Int)(floor(days))

	/* ----------------- find month and day of month ---------------- */
	if ( (year % 4) == 0 ) {
		lmonth[1] = 29
	}

	i = 1
	inttemp = 0
	while ((dayofyr > inttemp + lmonth[i-1]) && (i < 12)) {
		inttemp = inttemp + lmonth[i-1]
		i = i + 1
	}
	out_mon = i
	out_day = dayofyr - inttemp

	/* ----------------- find hours minutes and seconds ------------- */
	temp = (days - Double(dayofyr)) * 24.0
	out_hr  = Int(floor(temp))
	temp = (temp - Double(out_hr)) * 60.0
	out_minute  = Int(floor(temp))
	out_sec = (temp - Double(out_minute)) * 60.0
	return (out_mon, out_day, out_hr, out_minute, out_sec)
}  // end days2mdhms

/* -----------------------------------------------------------------------------
*
*                           procedure invjday
*
*  this procedure finds the year, month, day, hour, minute and second
*  given the julian date. tu can be ut1, tdt, tdb, etc.
*
*  algorithm     : set up starting values
*                  find leap year - use 1900 because 2000 is a leap year
*                  find the elapsed days through the year in a loop
*                  call routine to find each individual value
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    jd          - julian date                    days from 4713 bc
*
*  outputs       :
*    year        - year                           1900 .. 2100
*    mon         - month                          1 .. 12
*    day         - day                            1 .. 28,29,30,31
*    hr          - hour                           0 .. 23
*    min         - minute                         0 .. 59
*    sec         - second                         0.0 .. 59.999
*
*  locals        :
*    days        - day of year plus fractional
*                  portion of a day               days
*    tu          - julian centuries from 0 h
*                  jan 0, 1900
*    temp        - temporary double values
*    leapyrs     - number of leap years from 1900
*
*  coupling      :
*    days2mdhms  - finds month, day, hour, minute and second given days and year
*
*  references    :
*    vallado       2007, 208, alg 22, ex 3-13
* --------------------------------------------------------------------------- */

func  invjdayS (_ jd : Double) -> (out_year : Int, out_mon : Int, out_day : Int, out_hr : Int, out_minute : Int, out_sec : Double) {

	var out_year : Int
	var out_mon : Int
	var out_day : Int
	var out_hr : Int
	var out_minute : Int
	var out_sec : Double

	var leapyrs : Int
	var days, tu, temp : Double

	/* --------------- find year and days of the year --------------- */
	temp    = jd - 2415019.5
	tu      = temp / 365.25
	out_year = 1900 + Int(floor(tu))
	leapyrs = Int(floor(( Double(out_year) - 1901) * 0.25))

	// optional nudge by 8.64x10-7 sec to get even outputs
	days    = temp - Double(( out_year - 1900) * 365 + leapyrs) + 0.00000000001

	/* ------------ check for case of beginning of a year ----------- */
	if (days < 1.0) {
		out_year    =  out_year - 1
		leapyrs = Int(floor(( Double(out_year) - 1901) * 0.25))
		days = temp - Double(( out_year - 1900) * 365 + leapyrs)
	}

	/* ----------------- find remaing data  ------------------------- */
	(out_mon, out_day, out_hr, out_minute, out_sec) = days2mdhmsS(out_year, days)
	out_sec =  out_sec - 0.00000086400
	return (out_year, out_mon, out_day, out_hr, out_minute, out_sec)
}  // end invjday





