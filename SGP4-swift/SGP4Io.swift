//
//  SGP4Io.swift
//  pxSat3D-SK
//
//  Created by Xav perso on 05/08/2020.
//  Copyright © 2020 P'tit Xav. All rights reserved.
//

/*     ----------------------------------------------------------------
*
*                               sgp4io.cpp
*
*    this file contains a function to read two line element sets. while
*    not formerly part of the sgp4 mathematical theory, it is
*    required for practical implemenation.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                    2007
*                              by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*               3 sep 08  david vallado
*                           add operationmode for afspc (a) or improved (i)
*    changes :
*               9 may 07  david vallado
*                           fix year correction to 57
*              27 mar 07  david vallado
*                           misc fixes to manual inputs
*              14 aug 06  david vallado
*                           original baseline
*       ----------------------------------------------------------------      */

/* -----------------------------------------------------------------------------
*
*                           function twoline2rv
*
*  this function converts the two line element set character string data to
*    variables and initializes the sgp4 variables. several intermediate varaibles
*    and quantities are determined. note that the result is a structure so multiple
*    satellites can be processed simultaneously without having to reinitialize. the
*    verification mode is an important option that permits quick checks of any
*    changes to the underlying technical theory. this option works using a
*    modified tle file in which the start, stop, and delta time values are
*    included at the end of the second line of data. this only works with the
*    verification mode. the catalog mode simply propagates from -1440 to 1440 min
*    from epoch and is useful when performing entire catalog runs.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs        :
*    longstr1    - first line of the tle
*    longstr2    - second line of the tle
*    typerun     - type of run                    verification 'v', catalog 'c',
*                                                 manual 'm'
*    typeinput   - type of manual input           mfe 'm', epoch 'e', dayofyr 'd'
*    opsmode     - mode of operation afspc or improved 'a', 'i'
*    whichconst  - which set of constants to use  72, 84
*
*  outputs       :
*    satrec      - structure containing all the sgp4 satellite information
*
*  coupling      :
*    getGravityConstants-
*    days2mdhms  - conversion of days to month, day, hour, minute, second
*    jday        - convert day month year hour minute second into julian date
*    sgp4init    - initialize the sgp4 variables
*
*  references    :
*    norad spacetrack report #3
*    vallado, crawford, hujsak, kelso  2006
--------------------------------------------------------------------------- */
import Foundation

func stringToCars(_ chaine : String ) -> [Character]{
	var cars : [Character] = []
	for car in chaine {
		cars.append(car)
	}
	return cars
}

func carsToString(_ cars : [Character]) -> String {
	var str = ""
	for car in cars {
		str.append(car)
	}
	return str
}

func carsToString(_ cars : [Character], from : Int, to : Int) -> String {
	var str = ""
	for i in from...to {
		str.append(cars[i])
	}
	return str
}

func lineCarsToString(_ cars : [Character], from : Int, to : Int) -> String {
	return carsToString(cars, from: from-1, to: to-1)
}

func lineCarsToInt(_ cars : [Character], from : Int, to : Int) -> Int {
	return Int(lineCarsToString(cars, from: from, to: to)) ?? 0
}

func lineCarsToChar(_ cars : [Character], from : Int) -> Character {
	return cars[from-1]
}

func lineCarsToDouble(_ cars : [Character], from : Int, to : Int) -> Double {
	let str = lineCarsToString(cars, from: from, to: to)
	var lcars = stringToCars(str)
	for i in 0..<lcars.count {
		if lcars[i] == " " {
			lcars[i] = "0"
		}
	}
	let lstr = carsToString(lcars)
	let v = Double(lstr)
//	print("lineCarsToDouble \(str) -> \(v ?? -1)")
	return v ?? 0.0
}

func lineCarsToFloat(_ cars : [Character], from : Int, to : Int) -> Float {
	return Float(lineCarsToString(cars, from: from, to: to)) ?? 0.0
}

func nextValueFromCars(_ cars : [Character], from : Int) -> (value: String, to : Int) {
	var end = false
	var to = from

	while !end {
		if to < cars.count {
			if cars[to] == " " {
				to = to + 1
			} else {
				end = true
			}
		} else {
			end = true
			return ("", to)
		}
	}

	var value : [Character] = []

	while !end {
		if to < cars.count {
			if cars[to] != " " {
				value.append(cars[to])
			} else {
				end = true
			}
			to = to + 1
		} else {
			end = true
		}
	}

	return (carsToString(value), to)
}

func nextDoubleFromCars(_ cars : [Character], from : Int) -> (value: Double, to : Int) {
	let (value, to) = nextValueFromCars(cars, from: from)
	return (Double(value) ?? 0.0, to)
}

func nextIntFromCars(_ cars : [Character], from : Int) -> (value: Int, to : Int) {
	let (value, to) = nextValueFromCars(cars, from: from)
	return (Int(value) ?? 0, to)
}

func nextCharFromCars(_ cars : [Character], from : Int) -> (value: Character, to : Int) {
	let (value, to) = nextValueFromCars(cars, from: from)
	return (stringToCars(value)[0], to)
}


func twoline2rv(_ in_longstr1 : String, _ in_longstr2: String,
				_ in_typerun : Character, _ in_typeinput : Character, _ in_opsmode : Character,
				_ in_whichconst : GravityConstantsType) -> (
	satrec : Elsetrec,
	starfme : Double,
	stopmfe : Double,
	deltamin : Double)
{
	var out_satrec : Elsetrec = kNullElsetrec
	var out_startmfe = 0.0, out_stopmfe = 0.0, out_deltamin = 0.0

	var longstr1 = stringToCars(in_longstr1)
	var longstr2 = stringToCars(in_longstr2)

	var sec = 0.0
//	var startsec = 0.0, stopsec = 0.0, startdayofyr = 0.0, stopdayofyr = 0.0, jdstart = 0.0, jdstop = 0.0
//	var startyear, stopyear, startmon, stopmon, startday, stopday,
//	starthr, stophr, startmin, stopmin : Int
	var year = 0
	var mon = 0, day = 0, hr = 0, minute = 0
	var nexp, ibexp : Float
//	var cksum1, cksum2 : Int

	let gravityConstants = GravityConstants(type: in_whichconst)

	out_satrec.error = .success

	// set the implied decimal points since doing a formated read
	// fixes for bad input data values (missing, ...)
	for j in 10...15 {
		if longstr1[j] == " " {
			longstr1[j] = "_"
		}
	}

	if longstr1[44] != " " {
		longstr1[43] = longstr1[44]
	}
	longstr1[44] = "."
	if longstr1[7] == " " {
		longstr1[7] = "U"
	}
	if longstr1[9] == " " {
		longstr1[9] = "."
	}
	for j in 45...49 {
		if longstr1[j] == " " {
			longstr1[j] = "0"
		}
	}
	if longstr1[51] == " " {
		longstr1[51] = "0"
	}
	if longstr1[53] != " " {
		longstr1[52] = longstr1[53]
	}
	longstr1[53] = "."
	longstr2[25] = "."
	for j in 26...32 {
		if longstr2[j] == " "{
			longstr2[j] = "0"
		}
	}
	if longstr1[62] == " " {
		longstr1[62] = "0"
	}
	if longstr1[68] == " " {
		longstr1[68] = "0"
	}

//	sscanf(in_longstr1,
//		   "%2d %5ld %1c %10s %2d %12lf %11lf %7lf %2d %7lf %2d %2d %6ld ",
//		   &cardnumb,
//		   &out_satrec.satnum,
//		   &classification,
//		   intldesg,
//		   &out_satrec.epochyr,
//		   &out_satrec.epochdays,
//		   &out_satrec.ndot,
//		   &out_satrec.nddot,
//		   &nexp,
//		   &out_satrec.bstar,
//		   &ibexp,
//		   &numb,
//		   &elnum )

	//	Extract Contents of Line 1
	//	01 	Line Number of Element Data
	let cardnumb1 = lineCarsToInt(longstr1, from: 1, to: 1)
//	print("Line 1 type \(cardnumb1)")
	if cardnumb1 != 1 {
		print("Invalid card num line 1 \(cardnumb1)")
		out_satrec.error = .invalidTle
		return (out_satrec, 0.0, 0.0, 0.0)
	}

	//	03-07 	Satellite Number
	out_satrec.satnum = lineCarsToInt(longstr1, from: 3, to: 7)

	//	08 	Classification (U=Unclassified)
//	let classification = lineCarsToChar(longstr1, from: 8)
//	print("Line 1 Classification \(classification)")

	//	10-11 	International Designator (Last two digits of launch year)
	//	12-14 	International Designator (Launch number of the year)
	//	15-17 	International Designator (Piece of the launch)
//	let intldesg = lineCarsToString(longstr1, from: 10, to: 17)
//	print("Line 1 International Designator \(intldesg)")

	//	19-20 	Epoch Year (Last two digits of year)
	out_satrec.epochyr = lineCarsToInt(longstr1, from: 19, to: 20)

	//	21-32 	Epoch (Day of the year and fractional portion of the day)
	out_satrec.epochdays = lineCarsToDouble(longstr1, from: 21, to: 32)

	//	34-43 	First Time Derivative of the Mean Motion
	out_satrec.ndot = lineCarsToDouble(longstr1, from: 34, to: 43)

	//	45-52 	Second Time Derivative of Mean Motion (decimal point assumed)
	out_satrec.nddot = lineCarsToDouble(longstr1, from: 45, to: 50)
	nexp = lineCarsToFloat(longstr1, from: 51, to: 52)

	//	54-61 	BSTAR drag term (decimal point assumed)
	out_satrec.bstar = lineCarsToDouble(longstr1, from: 54, to: 59)
	ibexp = lineCarsToFloat(longstr1, from: 60, to: 61)

	//	63 	Ephemeris type
//	let numb = lineCarsToInt(longstr1, from: 63, to: 63)
//	print("Line 1 Ephemeris type \(numb)")

	//	65-68 	Element number
//	let elnum = lineCarsToInt(longstr1, from: 65, to: 68)
//	print("Line 1 Element number \(elnum)")

	//	69 	Checksum (Modulo 10)
//	cksum1 = lineCarsToInt(longstr1, from: 69, to: 69)
//	print("Line 1 Checksum \(cksum1)")

	//	Line 2
	//	Column 	Description
	//	01 	Line Number of Element Data
	let cardnumb2 = lineCarsToInt(longstr2, from: 1, to: 1)
//	print("Line 2 type \(cardnumb2)")
	if cardnumb2 != 2 {
		print("Invalid card num line 2 \(cardnumb2)")
		out_satrec.error = .invalidTle
		return (out_satrec, 0.0, 0.0, 0.0)
	}
	//	03-07 	Satellite Number
	out_satrec.satnum = lineCarsToInt(longstr2, from: 3, to: 7)

	//	09-16 	Inclination [Degrees]
	out_satrec.inclo = lineCarsToDouble(longstr2, from: 9, to:16)

	//	18-25 	Right Ascension of the Ascending Node [Degrees]
	out_satrec.nodeo = lineCarsToDouble(longstr2, from: 18, to: 25)

	//	27-33 	Eccentricity (decimal point assumed)
	out_satrec.ecco = lineCarsToDouble(longstr2, from: 26, to:33)

	//	35-42 	Argument of Perigee [Degrees]
	out_satrec.argpo = lineCarsToDouble(longstr2, from: 35, to: 42)

	//	44-51 	Mean Anomaly [Degrees]
	out_satrec.mo = lineCarsToDouble(longstr2, from: 44, to: 51)

	//	53-63 	Mean Motion [Revs per day]
	out_satrec.no = lineCarsToDouble(longstr2, from: 53, to: 63)

	//	64-68 	Revolution number at epoch [Revs]
//	let revnum = lineCarsToInt(longstr2, from: 64, to: 68)
//	print("Line 2 Revolution number \(revnum)")

	//	69 	Checksum (Modulo 10)
//	cksum2 = lineCarsToInt(longstr2, from: 69, to: 69)
//	print("Line 2 Checksum \(cksum2)")

//	if in_typerun == "v" {  // run for specified times from the file
//		if (longstr2[52] == " ") {
//			sscanf(in_longstr2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %10lf %6ld %lf %lf %lf \n",
//				   &cardnumb,
//				   &out_satrec.satnum,
//				   &out_satrec.inclo,
//				   &out_satrec.nodeo,
//				   &out_satrec.ecco,
//				   &out_satrec.argpo,
//				   &out_satrec.mo,
//				   &out_satrec.no,
//				   &revnum,
//				   &*out_startmfe,
//				   &*out_stopmfe,
//				   &*out_deltamin )
//		} else {
//			sscanf(in_longstr2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld %lf %lf %lf \n",
//				   &cardnumb,
//				   &out_satrec.satnum,
//				   &out_satrec.inclo,
//				   &out_satrec.nodeo,
//				   &out_satrec.ecco,
//				   &out_satrec.argpo,
//				   &out_satrec.mo,
//				   &out_satrec.no,
//				   &revnum,
//				   &*out_startmfe
//				, &*out_stopmfe,
//				  &*out_deltamin )
//		}
//	} else  {// simply run -1 day to +1 day or user input times
//		if longstr2[52] == " " {
//			sscanf(in_longstr2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %10lf %6ld \n",
//				   &cardnumb,
//				   &out_satrec.satnum,
//				   &out_satrec.inclo,
//				   &out_satrec.nodeo,
//				   &out_satrec.ecco,
//				   &out_satrec.argpo,
//				   &out_satrec.mo,
//				   &out_satrec.no,
//				   &revnum )
//		} else {
//			sscanf(in_longstr2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld \n",
//				   &cardnumb,
//				   &out_satrec.satnum,
//				   &out_satrec.inclo,
//				   &out_satrec.nodeo,
//				   &out_satrec.ecco,
//				   &out_satrec.argpo,
//				   &out_satrec.mo,
//				   &out_satrec.no,
//				   &revnum )
//		}
//	}

	// ---- find no, ndot, nddot ----
	out_satrec.no   = out_satrec.no / xpdotp //* rad/min
	out_satrec.nddot = out_satrec.nddot * Double(pow(10.0, nexp))
	out_satrec.bstar = out_satrec.bstar * Double(pow(10.0, ibexp))

	// ---- convert to sgp4 units ----
	out_satrec.a     = pow( out_satrec.no * gravityConstants.tumin , (-2.0/3.0) )
	out_satrec.ndot  = out_satrec.ndot  / (xpdotp * 1440.0)  //* ? * minperday
	out_satrec.nddot = out_satrec.nddot / (xpdotp * 1440.0*1440)

	// ---- find standard orbital elements ----
    out_satrec.inclo = out_satrec.inclo.degToRad
    out_satrec.nodeo = out_satrec.nodeo.degToRad
    out_satrec.argpo = out_satrec.argpo.degToRad
    out_satrec.mo    = out_satrec.mo.degToRad

	out_satrec.alta = out_satrec.a*(1.0 + out_satrec.ecco) - 1.0
	out_satrec.altp = out_satrec.a*(1.0 - out_satrec.ecco) - 1.0

	// ----------------------------------------------------------------
	// find sgp4epoch time of element set
	// remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
	// and minutes from the epoch (time)
	// ----------------------------------------------------------------

	// ---------------- temp fix for years from 1957-2056 -------------------
	// --------- correct fix will occur when year is 4-digit in tle ---------
	if (out_satrec.epochyr < 57) {
		year = out_satrec.epochyr + 2000
	} else {
		year = out_satrec.epochyr + 1900
	}

	(mon, day, hr, minute, sec) = days2mdhmsS(year,out_satrec.epochdays)
	out_satrec.jdsatepoch = jday( year,mon,day,hr,minute,sec)

//	// ---- input start stop times manually
//	if ((in_typerun != "v") && (in_typerun != "c"))
//	{
//		// ------------- enter start/stop ymd hms values --------------------
//		if (in_typeinput == "e")
//		{
//			printf("input start prop year mon day hr min sec \n")
//			// make sure there is no space at the end of the format specifiers in scanf!
//			scanf( "%i %i %i %i %i %lf",&startyear, &startmon, &startday, &starthr, &startmin, &startsec)
//			fflush(stdin)
//			jday( startyear,startmon,startday,starthr,startmin,startsec, &jdstart )
//
//			printf("input stop prop year mon day hr min sec \n")
//			scanf( "%i %i %i %i %i %lf",&stopyear, &stopmon, &stopday, &stophr, &stopmin, &stopsec)
//			fflush(stdin)
//			jday( stopyear,stopmon,stopday,stophr,stopmin,stopsec, &jdstop )
//
//			*out_startmfe = (jdstart - out_satrec.jdsatepoch) * 1440.0
//			*out_stopmfe  = (jdstop - out_satrec.jdsatepoch) * 1440.0
//
//			printf("input time step in minutes \n")
//			scanf( "%lf",&*out_deltamin )
//		}
//		// -------- enter start/stop year and days of year values -----------
//		if (in_typeinput == "d")
//		{
//			printf("input start year dayofyr \n")
//			scanf( "%i %lf",&startyear, &startdayofyr )
//			printf("input stop year dayofyr \n")
//			scanf( "%i %lf",&stopyear, &stopdayofyr )
//
//			days2mdhmsS( startyear,startdayofyr, &mon,&day,&hr,&minute,&sec )
//			jday( startyear,mon,day,hr,minute,sec, &jdstart )
//			days2mdhmsS( stopyear,stopdayofyr, &mon,&day,&hr,&minute,&sec )
//			jday( stopyear,mon,day,hr,minute,sec, &jdstop )
//
//			*out_startmfe = (jdstart - out_satrec.jdsatepoch) * 1440.0
//			*out_stopmfe  = (jdstop - out_satrec.jdsatepoch) * 1440.0
//
//			printf("input time step in minutes \n")
//			scanf( "%lf",&*out_deltamin )
//		}
//		// ------------------ enter start/stop mfe values -------------------
//		if (in_typeinput == "m")
//		{
//			printf("input start min from epoch \n")
//			scanf( "%lf",&*out_startmfe )
//			printf("input stop min from epoch \n")
//			scanf( "%lf",&*out_stopmfe )
//			printf("input time step in minutes \n")
//			scanf( "%lf",&*out_deltamin )
//		}
//	}

	// ------------ perform complete catalog evaluation, -+ 1 day -----------
	if in_typerun == "c" {
		out_startmfe = -1440.0
		out_stopmfe  =  1440.0
		out_deltamin =    10.0
	}

	// ---------------- initialize the orbit at sgp4epoch -------------------
	out_satrec = sgp4init( in_whichconst, in_opsmode, out_satrec.satnum, out_satrec.jdsatepoch - 2433281.5, out_satrec.bstar,
						   out_satrec.ecco, out_satrec.argpo, out_satrec.inclo, out_satrec.mo, out_satrec.no,
						   out_satrec.nodeo,
						   out_satrec)

	return (out_satrec, out_startmfe, out_stopmfe, out_deltamin)
} // end twoline2rv


