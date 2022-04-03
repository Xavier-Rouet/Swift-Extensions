//
//  SGP4Model.swift
//  pxSat3D-SK
//
//  Created by Xav perso on 05/08/2020.
//  Copyright © 2020 P'tit Xav. All rights reserved.
//

import Foundation

class SGP4Model {
	var gravityConstantType = GravityConstantsType.wgs84
	var satRec = kNullElsetrec
	var stationInfo = kVector4Null
	var stationInfoGeodetic = kVector4Null
	var initStatus : SGP4Status = .undefined

	func setupSatrec(gravityConstantType inGravityConstantType : GravityConstantsType,
					 operationMode inOperationMode : Character,
					 line1 inLine1 : String,
					 line2 inLine2 : String) {
		gravityConstantType = inGravityConstantType

		(satRec, _, _, _) = twoline2rv (inLine1, inLine2, "c",  "e", inOperationMode, gravityConstantType)

		/*
		print("twoline2rv ->")
		print("satnum %ld\n", satRec.satnum)
		print("epochyr %d\n", satRec.epochyr)
		print("epochtynumrev %d\n", satRec.epochtynumrev)
		print("error %d\n", satRec.error)
		print("operationmode %c\n", satRec.operationmode)
		print("init %c\n", satRec.inited)
		print("a %le\n", satRec.a)
		print("altp %le\n", satRec.altp)
		print("alta %le\n", satRec.alta)
		print("epochdays %le\n", satRec.epochdays)
		print("jdsatepoch %le\n", satRec.jdsatepoch)
		print("nddot %le\n", satRec.nddot)
		print("ndot %le\n", satRec.ndot)
		print("bstar %le\n", satRec.bstar)
		print("rcse %le\n", satRec.rcse)
		print("inclo %le\n", satRec.inclo)
		print("nodeo %le\n", satRec.nodeo)
		print("ecco %le\n", satRec.ecco)
		print("argpo %le\n", satRec.argpo)
		print("mo %le\n", satRec.mo)
		*/

		satRec =  sgp4init(gravityConstantType, inOperationMode, satRec.satnum, satRec.jdsatepoch-2433281.5, satRec.bstar,
						   satRec.ecco, satRec.argpo, satRec.inclo, satRec.mo, satRec.no,
						   satRec.nodeo,
						   satRec)
		initStatus = satRec.error

		/*
		print("sgp4init ->")
		print("satnum %ld\n", satRec.satnum)
		print("epochyr %d\n", satRec.epochyr)
		print("epochtynumrev %d\n", satRec.epochtynumrev)
		print("error %d\n", satRec.error)
		print("operationmode %c\n", satRec.operationmode)
		print("init %c\n", satRec.inited)
		print("a %le\n", satRec.a)
		print("altp %le\n", satRec.altp)
		print("alta %le\n", satRec.alta)
		print("epochdays %le\n", satRec.epochdays)
		print("jdsatepoch %le\n", satRec.jdsatepoch)
		print("nddot %le\n", satRec.nddot)
		print("ndot %le\n", satRec.ndot)
		print("bstar %le\n", satRec.bstar)
		print("rcse %le\n", satRec.rcse)
		print("inclo %le\n", satRec.inclo)
		print("nodeo %le\n", satRec.nodeo)
		print("ecco %le\n", satRec.ecco)
		print("argpo %le\n", satRec.argpo)
		print("mo %le\n", satRec.mo)
		*/
	}

	func setStationLatitude(_ inLatitude : Double,
							longitude inLongitude : Double,
							altitude inAltitude : Double,
							declinaison inDeclinaison  : Double) {
		stationInfo[0] = inLatitude
		stationInfo[1] = inLongitude
		stationInfo[2] = inAltitude
		stationInfo[3] = inDeclinaison
        stationInfoGeodetic[0] = stationInfo[0].degToRad
		stationInfoGeodetic[1] = stationInfo[1].degToRad
		stationInfoGeodetic[2] = stationInfo[2] / 1000.0
		stationInfoGeodetic[3] = 0.0 // calcule par les procedures appelees en fonction du temps et de la longitude
	}

	init(gravityConstantType inGravityConstantType : GravityConstantsType,
		 operationMode inOperationMode : Character,
		 line1 inLine1 : String,
		 line2 inLine2 : String) {
		self.setupSatrec(gravityConstantType: inGravityConstantType,
						 operationMode: inOperationMode,
						 line1: inLine1,
						 line2: inLine2)
	}

	init(WithGravityConstantType inGravityConstantType : GravityConstantsType,
		 operationMode inOperationMode : Character,
		 line1 inLine1 : String,
		 line2 inLine2 : String,
		 latitude inLatitude : Double,
		 longitude inLongitude : Double,
		 altitude inAltitude : Double,
		 declinaison inDeclinaison : Double) {
		self.setStationLatitude(inLatitude,
								longitude: inLongitude,
								altitude: inAltitude,
								declinaison: inDeclinaison)
		self.setupSatrec(gravityConstantType: inGravityConstantType,
						 operationMode: inOperationMode,
						 line1: inLine1,
						 line2: inLine2)
	}

	func computePosition(atEpoch tsince : Double) -> (sgp4Status : SGP4Status, position : Position, vitesse : Velocity) {
		// tsince in minutes since ephemeris time
		let lEpoch = (tsince - satRec.jdsatepoch) * 1440.0
		let (sgp4Status, position, vitesse) = sgp4(gravityConstantType, satRec, lEpoch)

		if sgp4Status.isError() {
			print("computePosition sgp4Status=\(sgp4Status.description)")
		}

		return (sgp4Status, position, vitesse)
	}

	func computePosition(atEpoch tsince : Double) -> (sgp4Status : SGP4Status, position : Position, vitesse : Velocity, visible : Bool) {
		var position = kVector3Null, vitesse = kVector3Null
		var sgp4Status : SGP4Status
		(sgp4Status, position, vitesse) = self.computePosition(atEpoch:tsince)
		var	lObservatorInfo  = kVector4Null
		let lIsvisible = eciToObservator(gravityConstantType, position, vitesse, &stationInfoGeodetic, tsince, &lObservatorInfo)
		return (sgp4Status, position, vitesse, lIsvisible)
	}


	func computePosition(_ outPosition : inout Position,
						 vitesse outVitesse : inout Velocity,
						 atEpoch tsince : Double) -> SGP4Status {
		// tsince in minutes since ephemeris time
		let lEpoch = (tsince - satRec.jdsatepoch) * 1440.0
		var sgp4Status : SGP4Status
		(sgp4Status, outPosition, outVitesse) = sgp4(gravityConstantType, satRec, lEpoch)

		if sgp4Status.isError() {
			print("computePosition sgp4Status=\(sgp4Status.description)")
		}


		return sgp4Status
	}

	func computePosition(_ outPosition : inout Position,
						 vitesse outVitesse : inout Velocity,
						 visible outVisible : inout Bool,
						 atEpoch tsince : Double) -> SGP4Status {

		let sgp4Status = self.computePosition(&outPosition, vitesse:&outVitesse, atEpoch:tsince)
		var	lObservatorInfo  = kVector4Null
		let lIsvisible = eciToObservator(gravityConstantType, outPosition, outVitesse, &stationInfoGeodetic, tsince, &lObservatorInfo)
		outVisible = lIsvisible
		return sgp4Status

	}

	func calculateLatitude(atTime inTime : Double,
						   forPosition inPosition : Position) -> (latitude : Double, longutude : Double, altitude : Double) {
		let lla = eciToGeodetic(gravityConstantType, inPosition, inTime)
        return (lla[0].radToDeg, lla[1].radToDeg, lla[2])
	}

	static func julianDayFromDate(_ inDateComponents : DateComponents) -> Double {
		return jday(inDateComponents.year!, inDateComponents.month!, inDateComponents.day!,
					inDateComponents.hour!, inDateComponents.minute!, Double(inDateComponents.second!)
		)
	}

	func julianDayFromDate(_ inDateComponents : DateComponents) -> Double {
		return SGP4Model.julianDayFromDate(inDateComponents)
	}

	var tleEpochDate : Double {
		return satRec.jdsatepoch
	}

	static func dateFromJulianDay(_ inJulianDay : Double) -> DateComponents {
		var year = 0,
		mon = 0,
		day = 0,
		hr = 0,
		minute = 0
		var sec = 0.0
		var dateComponents : DateComponents

		(year, mon, day, hr, minute, sec) = invjdayS(inJulianDay)

		dateComponents = DateComponents()
		dateComponents.year = year
		dateComponents.month = mon
		dateComponents.day = day
		dateComponents.hour = hr
		dateComponents.minute = minute
		dateComponents.second = Int(floor(sec))
		return dateComponents
	}

	//	static func positionSoleilEciAtTime(_ julianDate : Double,
	//										position solar_vector : inout Position) {
	//		solar_vector = positionSoleilEci(julianDate)
	//	}

	func computePosition(_ outPosition :  inout Position,
						 vitesse outVitesse : inout Velocity,
						 latitude outLatitude : inout Double,
						 longitude outLongitude : inout Double,
						 altitude outAltitude : inout Double,
						 visible outVisible : inout Bool,
						 azimuth outAzimuth : inout Double,
						 elevation outElevation : inout Double,
						 distance outDistance : inout Double,
						 vitesseRelative outVitesseRelative : inout Double,
						 rightAscension outRightAscension : inout Double,
						 declination outDeclination : inout Double,
						 atTime inTime : Double)
	{
		var	lObservatorInfo = kVector4Null
		var	lAstronomicInfo = kVector4Null
		var lIsvisible : Bool

		let sgp4Status = self.computePosition(&outPosition,
											  vitesse:&outVitesse,
											  atEpoch:inTime)

		if sgp4Status.isError() {
			outAzimuth = -1.0
			outElevation = -1.0
			outDistance = -1.0
			outVitesseRelative = -1.0
			outVisible = false
			outRightAscension = -1.0
			outDeclination = -1.0
		} else {
			(outLatitude, outLongitude, outAltitude) = self.calculateLatitude(atTime:inTime,
																			  forPosition:outPosition)
		}

		lIsvisible = eciToObservator(gravityConstantType,
									 outPosition,
									 outVitesse,
									 &stationInfoGeodetic,
									 inTime,
									 &lObservatorInfo)
		outAzimuth = lObservatorInfo[0].radToDeg
		outElevation = lObservatorInfo[1].radToDeg
		outDistance = lObservatorInfo[2]
		outVitesseRelative = lObservatorInfo[3]

		outVisible = lIsvisible

		// note : lAstronomicInfo is used in a call to eciToObservator
		lAstronomicInfo = lObservatorInfo
		eciToAstronomic(gravityConstantType,
						outPosition,
						outVitesse,
						&stationInfoGeodetic,
						inTime,
						&lAstronomicInfo)
        outRightAscension = lAstronomicInfo[0].radToHour // sur 24 heures
		outDeclination = lAstronomicInfo[1].radToDeg

		// On enleve la declinaison magnetique
		outAzimuth -= stationInfo[3]
	}

	func visiblePosition(_ outPosition : inout Position,
						 vitesse outVitesse : inout Velocity,
						 latitude outLatitude : inout Double,
						 longitude outLongitude : inout Double,
						 altitude outAltitude : inout Double,
						 azimuth outAzimuth : inout Double,
						 elevation outElevation : inout Double,
						 distance outDistance : inout Double,
						 vitesseRelative outVitesseRelative : inout Double,
						 rightAscension outRightAscension : inout Double,
						 declination outDeclination : inout Double,
						 atTime inTime : Double) -> Bool
	{
		var visible = false
		self.computePosition(&outPosition,
							 vitesse:&outVitesse,
							 latitude:&outLatitude,
							 longitude:&outLongitude,
							 altitude:&outAltitude,
							 visible:&visible,
							 azimuth:&outAzimuth,
							 elevation:&outElevation,
							 distance:&outDistance,
							 vitesseRelative:&outVitesseRelative,
							 rightAscension:&outRightAscension,
							 declination:&outDeclination,
							 atTime:inTime)
		return visible
	}
}
