//
//  SGP4Lib.swift
//  pxSat3D-SK
//
//  Created by Xav perso on 03/08/2020.
//  Copyright © 2020 P'tit Xav. All rights reserved.
//

/*
*  sgp4lib.c
*  pxTleManagement
*
*  Created by Xav perso on 07/02/09.
*  Copyright 2009 P'tit Xav SAS. All rights reserved.
*
*/
import Foundation



func thetaGJD(_ in_jd : Double) -> Double {
	var jd = in_jd
	// Reference:  The 1992 Astronomical Almanac, page B6.
	var UT,TU,GMST : Double
	let secday : Double  = 86400.0        // Seconds per day
	let omega_E : Double = 1.00273790934  // Earth rotations per sidereal day (non-constant)


//	UT   = Modulus(jd + 0.5,1.0)
	UT   = (jd + 0.5) % 1.0
	jd   = jd - UT
	TU   = (jd - 2451545.0)/36525.0
	GMST = 24110.54841 + TU * (8640184.812866 + TU * (0.093104 - TU * 6.2E-6))
//	GMST = Modulus(GMST + secday*omega_E*UT,secday)
	GMST = (GMST + secday*omega_E*UT) % secday
	return Double.twoπ * GMST/secday
}

func thetaGJD(_ in_jd : Double, _ lon : Double) -> Double {
	return (thetaGJD(in_jd) + lon) %  Double.twoπ
}

func Delta_ET(_ year : Double) -> Double
{
	// { Values determined using data from 1950-1991 in the 1990 Astronomical
	// Almanac.  See DELTA_ET.WQ1 for details. }

	return 26.465 + 0.747622*(year - 1950) + 1.886913 * Sin(Double.twoπ*(year - 1975)/33)
}


func eciToGeodetic(_ whichconst : GravityConstantsType, _ pos : Position, _ time : Double) -> Observator {
	var geodetic = kVector4Null

	// Reference:  The 1992 Astronomical Almanac, page K12.
	var lat,lon,alt, theta,r,e2,phi,c : Double
	var dtheta : Double

	let gravityConstants = GravityConstants(type: whichconst)

	let f = 1.0/298.2572235630 // flattening of WGS 84 Ellipsoid

	theta = ArcTan(pos[1],pos[0])
	//NSLog(@"theta : %lf ThetaG_JD : %lf",theta,ThetaG_JD(time))
	dtheta = theta - thetaGJD(time)

	//(test) lon = Modulus(theta - ThetaGJD(time),Double.twoπ)
	//NSLog(@"lon : %lf",lon)

	while (dtheta > π) {
		dtheta -= Double.twoπ
	}
	while (dtheta < -π) {
		dtheta += Double.twoπ
	}
	lon = dtheta
	//NSLog(@"lon : %lf",lon)
	r = √(Sqr(pos[0]) + Sqr(pos[1]))
	e2 = f*(2 - f)
	lat = ArcTan(pos[2],r)
	repeat {
		phi = lat
		c = 1 / √(1 - e2*Sqr(Sin(phi)))
		lat = ArcTan(pos[2] + gravityConstants.radiusearthkm * c * e2 * Sin(phi),r)
		//NSLog(@"lat : %lf",lat)
	} while (Abs(lat - phi) >= 1E-10)

	alt = r/Cos(lat) - gravityConstants.radiusearthkm*c
	geodetic[0] = lat   // radians
	geodetic[1] = lon   // radians
	geodetic[2] = alt   // kilometers
	geodetic[3] = theta // radians
	return geodetic
}

func userToEci(_ whichconst : GravityConstantsType, _ geodetic : inout Observator, _ time : Double) -> (obs_pos : Position, obs_vel : Velocity) {
	var obs_pos = kVector3Null
	var obs_vel = kVector3Null

	let secday   = 86400.0        // Seconds per day
	let omega_E  = 1.00273790934  // Earth rotations per sidereal day (non-constant)

	let gravityConstants = GravityConstants(type: whichconst)

	let f        = 1.0/298.2572235630 // flattening of WGS 84 Ellipsoid

	let mfactor = Double.twoπ * omega_E / secday
	var	lat,lon,alt,theta,c,s,achcp : Double

	lat = geodetic[0]
	lon = geodetic[1]
	alt = geodetic[2]

//	theta = Modulus(thetaGJD(time) + lon,Double.twoπ)
	theta = thetaGJD(time,lon)
	geodetic[3] = theta
	c = 1 / √(1 + f*(f - 2)*Sqr(Sin(lat)))
	s = Sqr(1 - f)*c
	achcp = (gravityConstants.radiusearthkm * c + alt)*Cos(lat)
	obs_pos[0] = achcp * Cos(theta)         //  {kilometers}
	obs_pos[1] = achcp * Sin(theta)
	obs_pos[2] = (gravityConstants.radiusearthkm * s + alt) * Sin(lat)
//	Magnitude(&obs_pos)
	obs_vel[0] = -mfactor*obs_pos[1]       // {kilometers/second}
	obs_vel[1] =  mfactor*obs_pos[0]
	obs_vel[2] =  0
//	Magnitude(&obs_vel)

	return (obs_pos, obs_vel)
}

func eciToObservator(_ whichconst : GravityConstantsType,
					 _ pos : Position,
					 _ vel : Velocity,
					 _ geodetic : inout Observator,
					 _ time : Double,
					 _ obs_set : inout Observator) -> Bool {

	// double alt,lon // pas utilise
	var lat,theta : Double
	var sin_lat,cos_lat : Double
	var sin_theta,cos_theta : Double
	var el,azim : Double
	var top_s,top_e,top_z : Double
//	var range = kVector4Null
//	var rgvel = kVector4Null
	var visible : Bool

	let (obs_pos, obs_vel) = userToEci(whichconst,&geodetic,time)
	let range = pos - obs_pos
	let rgvel = vel - obs_vel

//	for i in 0...2
//	{
//		range[i] = pos[i] - obs_pos[i]
//		rgvel[i] = vel[i] - obs_vel[i]
//	}

//	Magnitude(&range)

	lat = geodetic[0]
	// lon = geodetic[1] // analyse : contenu pas lu
	// alt = geodetic[2] // analyse : contenu pas lu
	theta = geodetic[3]

	sin_lat = Sin(lat)
	cos_lat = Cos(lat)
	sin_theta = Sin(theta)
	cos_theta = Cos(theta)

	top_s = sin_lat*cos_theta*range[0] + sin_lat*sin_theta*range[1] - cos_lat*range[2]
	top_e = -sin_theta*range[0] + cos_theta*range[1]
	top_z = cos_lat*cos_theta*range[0] + cos_lat*sin_theta*range[1] + sin_lat*range[2]
	azim = ArcTan(-top_e/top_s)

	if (top_s > 0.0) {
		azim = azim + π
	}

	if (azim < 0.0) {
		azim = azim + Double.twoπ
	}

	el = ArcSin(top_z/range.mag)

	obs_set[0] = azim                      //{Azimuth (radians)}
	obs_set[1] = el                        //{Elevation (radians)}
	obs_set[2] = range.mag                  //{Range (kilometers)}
//	obs_set[3] = Dot(range,rgvel)/range[3] //{Range Rate (kilometers/second)}
	obs_set[3] = (range⋅rgvel)/range.mag //{Range Rate (kilometers/second)}

	//{ Corrections for atmospheric refraction }
	//{ Reference:  Astronomical Algorithms by Jean Meeus, pp. 101-104 }
	//{ Note:  Correction is meaningless when apparent elevation is below horizon }
	obs_set[1] = obs_set[1] + Radians((1.02/Tan(Radians(Degrees(el)+10.3/(Degrees(el)+5.11))))/60.0)

	if (obs_set[1] >= 0)
	{
		visible = true
	}
	else
	{
		obs_set[1] = el  //{Reset to true elevation}
		visible = false
	}

	return visible
}

func eciToAstronomic(_ whichconst : GravityConstantsType,
					 _ pos : Position, _ vel : Velocity,
					 _ geodetic : inout Observator,
					 _ time : Double,
					 _ obs_set : inout Observator) {
	var phi,theta : Double
	var sin_theta,cos_theta : Double
	var sin_phi,cos_phi : Double
	var az,el : Double
	var Lxh,Lyh,Lzh : Double
	var Sx,Ex,Zx : Double
	var Sy,Ey,Zy : Double
	var Sz,Ez,Zz : Double
	var Lx,Ly,Lz : Double
	var cos_delta : Double
	var sin_alpha,cos_alpha : Double
	var visible : Bool

	visible = eciToObservator(whichconst,pos,vel,&geodetic,time,&obs_set)

	if (visible)
	{
		az = obs_set[0]
		el = obs_set[1]
		phi   = geodetic[0]
//		theta = Modulus(thetaGJD(time) + geodetic[1],Double.twoπ)
		theta = thetaGJD(time,geodetic[1])
		sin_theta = Sin(theta)
		cos_theta = Cos(theta)
		sin_phi = Sin(phi)
		cos_phi = Cos(phi)
		Lxh = -Cos(az)*Cos(el)
		Lyh =  Sin(az)*Cos(el)
		Lzh =  Sin(el)
		Sx = sin_phi*cos_theta
		Ex = -sin_theta
		Zx = cos_theta*cos_phi
		Sy = sin_phi*sin_theta
		Ey = cos_theta
		Zy = sin_theta*cos_phi
		Sz = -cos_phi
		Ez = 0
		Zz = sin_phi
		Lx = Sx*Lxh + Ex*Lyh + Zx*Lzh
		Ly = Sy*Lxh + Ey*Lyh + Zy*Lzh
		Lz = Sz*Lxh + Ez*Lyh + Zz*Lzh
		//printf("*** Lz=%f\n",Lz)
		obs_set[1] = ArcSin(Lz)                 //{Declination (radians)}
		cos_delta = √(1 - Sqr(Lz))
		//printf("*** cos_delta=%f\n",cos_delta)
		sin_alpha = Ly/cos_delta
		cos_alpha = Lx/cos_delta
		obs_set[0] = ArcTan(sin_alpha,cos_alpha) //{Right Ascension (radians)}
//		obs_set[0] = Modulus(obs_set[0],Double.twoπ)
		obs_set[0] = obs_set[0] % Double.twoπ
	}
}

