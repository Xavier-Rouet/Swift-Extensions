//
//  SGP4Sun.swift
//  pxSat3D-SK
//
//  Created by Xav perso on 03/08/2020.
//  Copyright © 2020 P'tit Xav. All rights reserved.
//

/*
*  sgp4sun.c
*  PxSat3D
*
*  Created by Xav perso on 26/10/10.
*  Copyright 2010 P'tit Xav SAS. All rights reserved.
*
* C tradcuction of SOLAR.PAS
* {           Author:  Dr TS Kelso }
* { Original Version:  1990 Jul 29 }
* { Current Revision:  1999 Nov 27 }
* {          Version:  1.30 }
* {        Copyright:  1990-1999, All Rights Reserved }
*/

let AU : Double = 1.49597870E8  // {Astronomical unit - kilometers (IAU 76)}
let secday : Double	= 86400.0        // Seconds per day

func positionSoleilEci(_ time : Double) -> Position
{
	var solar_vector : Position = kVector3Null
	var mjd,year,T,M,L,e,C,O,Lsa,nu,R,eps : Double
	mjd = time - 2415020.0
	year = 1900 + mjd/365.25
	T = (mjd + Delta_ET(year)/secday)/36525.0
//	M = Radians(Modulus(358.47583 + Modulus(35999.04975*T,360.0)
	M = Radians((358.47583 + ((35999.04975*T) % 360.0)
	- (0.000150 + 0.0000033*T)*Sqr(T)) % 360.0)
//	L = Radians(Modulus(279.69668 + Modulus(36000.76892*T,360.0)
	L = Radians((279.69668 + ((36000.76892*T) % 360.0)
	+ 0.0003025*Sqr(T)) % 360.0)
	e = 0.01675104 - (0.0000418 + 0.000000126*T)*T
	C = Radians((1.919460 - (0.004789 + 0.000014*T)*T)*Sin(M)
	+ (0.020094 - 0.000100*T)*Sin(2*M) + 0.000293*Sin(3*M))
//	O = Radians(Modulus(259.18 - 1934.142*T,360.0))
	O = Radians((259.18 - 1934.142*T) % 360.0)
//	Lsa = Modulus(L + C - Radians(0.00569 - 0.00479*Sin(O)),Double.twoπ)
	Lsa = (L + C - Radians(0.00569 - 0.00479*Sin(O))) % Double.twoπ
//	nu = Modulus(M + C,Double.twoπ)
	nu = (M + C) % Double.twoπ
	R = 1.0000002*(1 - Sqr(e))/(1 + e*Cos(nu))
	eps = Radians(23.452294 - (0.0130125 + (0.00000164 - 0.000000503*T)*T)*T
	+ 0.00256*Cos(O))
	R = AU * R
	solar_vector[0] = R*Cos(Lsa)
	solar_vector[1] = R*Sin(Lsa)*Cos(eps)
	solar_vector[2] = R*Sin(Lsa)*Sin(eps)
	return solar_vector
}
