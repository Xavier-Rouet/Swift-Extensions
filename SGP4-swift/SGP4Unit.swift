//
// SGP4Unit.swift
//  pxSat3D-SK
//
//  Created by Xav perso on 03/08/2020.
//  Copyright © 2020 P'tit Xav. All rights reserved.
//

/*     ----------------------------------------------------------------
*
*                                 sgp4unit.h
*
*    this file contains the sgp4 procedures for analytical propagation
*    of a satellite. the code was originally released in the 1980 and 1986
*    spacetrack papers. a detailed discussion of the theory and history
*    may be found in the 2006 aiaa paper by vallado, crawford, hujsak,
*    and kelso.
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
*                           performance mode in elsetrec
*    changes :
*              20 apr 07  david vallado
*                           misc fixes for constants
*              11 aug 06  david vallado
*                           chg lyddane choice back to strn3, constants, misc doc
*              15 dec 05  david vallado
*                           misc fixes
*              26 jul 05  david vallado
*                           fixes for paper
*                           note that each fix is preceded by a
*                           comment with "sgp4fix" and an explanation of
*                           what was changed
*              10 aug 04  david vallado
*                           2nd printing baseline working
*              14 may 01  david vallado
*                           2nd edition baseline
*                     80  norad
*                           original baseline
*       ----------------------------------------------------------------      */
import Foundation

let SGP4Version = "SGP4 Version 2008-09-03"


// -------------------------- structure declarations ----------------------------
enum GravityConstantsType
{
	case wgs72old
	case wgs72
	case wgs84
}

/* -----------------------------------------------------------------------------
*
*  onstants for the propagator. note that mu is identified to
*    facilitiate comparisons with newer models. the common useage is wgs72.
*
*  author        : david vallado                  719-573-2600   21 jul 2006
*
*  inputs        :
*    gravityConstantsType  - which set of constants to use  wgs72old, wgs72, wgs84
*
*  outputs       :
*    tumin       - minutes in one time unit
*    mu          - earth gravitational parameter
*    radiusearthkm - radius of the earth in km
*    xke         - reciprocal of tumin
*    j2, j3, j4  - un-normalized zonal harmonic values
*    j3oj2       - j3 divided by j2
*
*  locals        :
*
*  coupling      :
*    none
*
*  references    :
*    norad spacetrack report #3
*    vallado, crawford, hujsak, kelso  2006
--------------------------------------------------------------------------- */
struct GravityConstants {
	var tumin : Double = 0.0
	var mu : Double = 0.0
	var radiusearthkm : Double = 0.0
	var xke : Double = 0.0
	var j2 : Double = 0.0
	var j3 : Double = 0.0
	var j4 : Double = 0.0
	var j3oj2 : Double = 0.0

	init(type gravityConstantsType : GravityConstantsType) {

		switch gravityConstantsType {
		// -- wgs-72 low precision str#3 constants --
		case .wgs72old:
			mu     = 398600.79964        // in km3 / s2
			radiusearthkm = 6378.135     // km
			xke    = 0.0743669161
			tumin  = 1.0 / xke
			j2     =   0.001082616
			j3     =  -0.00000253881
			j4     =  -0.00000165597
			j3oj2  =  j3 / j2
		// ------------ wgs-72 constants ------------
		case .wgs72:
			mu     = 398600.8            // in km3 / s2
			radiusearthkm = 6378.135     // km
			xke    = 60.0 / √(radiusearthkm * radiusearthkm * radiusearthkm / mu)
			tumin  = 1.0 / xke
			j2     =   0.001082616
			j3     =  -0.00000253881
			j4     =  -0.00000165597
			j3oj2  =  j3 / j2
		case .wgs84:
			// ------------ wgs-84 constants ------------
			mu     = 398600.5            // in km3 / s2
			radiusearthkm = 6378.137     // km
			xke    = 60.0 / √(radiusearthkm * radiusearthkm * radiusearthkm / mu)
			tumin  = 1.0 / xke
			j2     =   0.00108262998905
			j3     =  -0.00000253215306
			j4     =  -0.00000161098761
			j3oj2  =  j3 / j2
		}

	}
}

struct Elsetrec
{
	var satnum : Int
	var epochyr, epochtynumrev : Int
	var error : SGP4Status
	var operationmode : Character
	var inited : Bool
	var method  : Character

	/* Near Earth */
	var isimp : Int
	var aycof  , con41  , cc1    , cc4      , cc5    , d2      , d3   , d4    ,
	delmo  , eta    , argpdot, omgcof   , sinmao , t       , t2cof, t3cof ,
	t4cof  , t5cof  , x1mth2 , x7thm1   , mdot   , nodedot, xlcof , xmcof ,
	nodecf : Double

	/* Deep Space */
	var irez : Int
	var d2201  , d2211  , d3210  , d3222    , d4410  , d4422   , d5220 , d5232 ,
	d5421  , d5433  , dedt   , del1     , del2   , del3    , didt  , dmdt  ,
	dnodt  , domdt  , e3     , ee2      , peo    , pgho    , pho   , pinco ,
	plo    , se2    , se3    , sgh2     , sgh3   , sgh4    , sh2   , sh3   ,
	si2    , si3    , sl2    , sl3      , sl4    , gsto    , xfact , xgh2  ,
	xgh3   , xgh4   , xh2    , xh3      , xi2    , xi3     , xl2   , xl3   ,
	xl4    , xlamo  , zmol   , zmos     , atime  , xli     , xni : Double

	var a      , altp   , alta   , epochdays, jdsatepoch       , nddot , ndot  ,
	bstar  , rcse   , inclo  , nodeo    , ecco             , argpo , mo    ,
	no : Double
}

let kNullElsetrec = Elsetrec(satnum: 0,
							 epochyr: 0,
							 epochtynumrev: 0,
							 error: .success,
							 operationmode: " ",
							 inited: false,
							 method: " ",
							 isimp: 0,
							 aycof: 0.0, con41: 0.0, cc1: 0.0, cc4: 0.0, cc5: 0.0,
							 d2: 0.0, d3: 0.0, d4: 0.0, delmo: 0.0,
							 eta: 0.0, argpdot: 0.0, omgcof: 0.0, sinmao: 0.0,
							 t: 0.0, t2cof: 0.0, t3cof: 0.0, t4cof: 0.0, t5cof: 0.0,
							 x1mth2: 0.0, x7thm1: 0.0,
							 mdot: 0.0, nodedot: 0.0, xlcof: 0.0, xmcof: 0.0, nodecf: 0.0,
							 irez: 0,
							 d2201: 0.0, d2211: 0.0, d3210: 0.0, d3222: 0.0, d4410: 0.0, d4422: 0.0, d5220: 0.0, d5232: 0.0, d5421: 0.0, d5433: 0.0,
							 dedt: 0.0, del1: 0.0, del2: 0.0, del3: 0.0, didt: 0.0, dmdt: 0.0, dnodt: 0.0, domdt: 0.0,
							 e3: 0.0, ee2: 0.0,
							 peo: 0.0, pgho: 0.0, pho: 0.0, pinco: 0.0, plo: 0.0,
							 se2: 0.0, se3: 0.0, sgh2: 0.0, sgh3: 0.0, sgh4: 0.0, sh2: 0.0, sh3: 0.0,
							 si2: 0.0, si3: 0.0, sl2: 0.0, sl3: 0.0, sl4: 0.0,
							 gsto: 0.0, xfact: 0.0, xgh2: 0.0, xgh3: 0.0, xgh4: 0.0,
							 xh2: 0.0, xh3: 0.0, xi2: 0.0, xi3: 0.0, xl2: 0.0, xl3: 0.0, xl4: 0.0, xlamo: 0.0,
							 zmol: 0.0, zmos: 0.0, atime: 0.0, xli: 0.0, xni: 0.0,
							 a: 0.0, altp: 0.0, alta: 0.0,
							 epochdays: 0.0, jdsatepoch: 0.0, nddot: 0.0, ndot: 0.0, bstar: 0.0,
							 rcse: 0.0, inclo: 0.0, nodeo: 0.0, ecco: 0.0, argpo: 0.0, mo: 0.0, no: 0.0)

/*
*	SGP4 errors :
*                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
*                   2 - mean motion less than 0.0
*                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
*                   4 - semi-latus rectum < 0.0
*                   5 - epoch elements are sub-orbital
*                   6 - satellite has decayed
*
*/
enum SGP4Status {
	case success
	case meanElement(Double, Double)
	case meanMotion(Double)
	case perturbationElements(Double)
	case semiLatudRectum(Double)
	case epochElementAreSubOrbital(Double)
	case decayed
	case invalidTle
	case undefined

	func isSuccess() -> Bool {
		switch self {
		case .success:
			return true
		default:
			return false
		}
	}

	func isError() -> Bool {
		return !isSuccess()
	}

	func isDecayed() -> Bool {
		switch self {
		case .decayed:
			return true
		default:
			return false
		}
	}

	var description : String {
		switch self {
		case .success:
			return "success"
		case .meanElement(let ecc, let a):
			return String(format:"mean elements, ecc (%.4f) >= 1.0 or < -0.001 or a (%.3f) < 0.95 er",ecc,a)
		case .meanMotion(let nm):
			return String(format: "mean motion (%.2f) less than 0.0", nm)
		case .perturbationElements(let ecc):
			return String(format:"pert elements, ecc (%.2f) < 0.0  or > 1.0",ecc)
		case .semiLatudRectum(let slr):
			return String(format:"semi-latus rectum (%.3f) < 0.0",slr)
		case .epochElementAreSubOrbital(let rp):
			return String(format:"epoch elements are sub-orbital (rp=%.2f)",rp)
		case .decayed:
			return "object has decayed"
		case .invalidTle:
			return "TLE is not correctly formated"
		case .undefined:
			return "Status is not defined"
		}
	}
}

/*     ----------------------------------------------------------------
*
*                               sgp4unit.cpp
*
*    this file contains the sgp4 procedures for analytical propagation
*    of a satellite. the code was originally released in the 1980 and 1986
*    spacetrack papers. a detailed discussion of the theory and history
*    may be found in the 2006 aiaa paper by vallado, crawford, hujsak,
*    and kelso.
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
*                           fix atime for faster operation in dspace
*                           add operationmode for afspc (a) or improved (i)
*                           performance mode
*    changes :
*              16 jun 08  david vallado
*                           update small eccentricity check
*              16 nov 07  david vallado
*                           misc fixes for better compliance
*              20 apr 07  david vallado
*                           misc fixes for constants
*              11 aug 06  david vallado
*                           chg lyddane choice back to strn3, constants, misc doc
*              15 dec 05  david vallado
*                           misc fixes
*              26 jul 05  david vallado
*                           fixes for paper
*                           note that each fix is preceded by a
*                           comment with "sgp4fix" and an explanation of
*                           what was changed
*              10 aug 04  david vallado
*                           2nd printing baseline working
*              14 may 01  david vallado
*                           2nd edition baseline
*                     80  norad
*                           original baseline
*       ----------------------------------------------------------------      */

let help = "n"

/* -----------------------------------------------------------------------------
*
*                           procedure dpper
*
*  this procedure provides deep space long period periodic contributions
*    to the mean elements.  by design, these periodics are zero at epoch.
*    this used to be dscom which included initialization, but it's really a
*    recurring function.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    e3          -
*    ee2         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    se2 , se3 , sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4 -
*    t           -
*    xh2, xh3, xi2, xi3, xl2, xl3, xl4 -
*    zmol        -
*    zmos        -
*    ep          - eccentricity                           0.0 - 1.0
*    inclo       - inclination - needed for lyddane modification
*    nodep       - right ascension of ascending node
*    argpp       - argument of perigee
*    mp          - mean anomaly
*
*  outputs       :
*    ep          - eccentricity                           0.0 - 1.0
*    inclp       - inclination
*    nodep        - right ascension of ascending node
*    argpp       - argument of perigee
*    mp          - mean anomaly
*
*  locals        :
*    alfdp       -
*    betdp       -
*    cosip  , sinip  , cosop  , sinop  ,
*    dalf        -
*    dbet        -
*    dls         -
*    f2, f3      -
*    pe          -
*    pgh         -
*    ph          -
*    pinc        -
*    pl          -
*    sel   , ses   , sghl  , sghs  , shl   , shs   , sil   , sinzf , sis   ,
*    sll   , sls
*    xls         -
*    xnoh        -
*    zf          -
*    zm          -
*
*  coupling      :
*    none.
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/

func dpper
(
_ e3 : Double,     _ ee2 : Double,    _ peo : Double,     _ pgho : Double,   _ pho : Double,
_ pinco : Double,  _ plo : Double,    _ se2 : Double,     _ se3 : Double,    _ sgh2 : Double,
_ sgh3 : Double,   _ sgh4 : Double,   _ sh2 : Double,     _ sh3 : Double,    _ si2 : Double,
_ si3 : Double,    _ sl2 : Double,    _ sl3 : Double,     _ sl4 : Double,    _ t : Double,
_ xgh2 : Double,   _ xgh3 : Double,   _ xgh4 : Double,    _ xh2 : Double,    _ xh3 : Double,
_ xi2 : Double,    _ xi3 : Double,    _ xl2 : Double,     _ xl3 : Double,    _ xl4 : Double,
_ zmol : Double,   _ zmos : Double,   _ inclo : Double,
_ inited : Bool,
_ out_ep : inout Double, _ out_inclp : inout Double, _ out_nodep : inout Double, _ out_argpp : inout Double, _ out_mp : inout Double,
_ opsmode : Character
	) -> ()
{
	/* --------------------- local variables ------------------------ */
	var alfdp, betdp, cosip, cosop, dalf, dbet, dls,
	f2,    f3,    pe,    pgh,   ph,   pinc, pl ,
	sel,   ses,   sghl,  sghs,  shll, shs,  sil,
	sinip, sinop, sinzf, sis,   sll,  sls,  xls,
	xnoh,  zf,    zm,    zel,   zes,  znl,  zns : Double

	/* ---------------------- constants ----------------------------- */
	zns   = 1.19459e-5
	zes   = 0.01675
	znl   = 1.5835218e-4
	zel   = 0.05490

	/* --------------- calculate time varying periodics ----------- */
	zm    = zmos + zns * t
	// be sure that the initial call has time set to zero
	if inited {
		zm = zmos
	}
	zf    = zm + 2.0 * zes * sin(zm)
	sinzf = sin(zf)
	f2    =  0.5 * sinzf * sinzf - 0.25
	f3    = -0.5 * sinzf * cos(zf)
	ses   = se2 * f2 + se3 * f3
	sis   = si2 * f2 + si3 * f3
	sls   = sl2 * f2 + sl3 * f3 + sl4 * sinzf
	sghs  = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf
	shs   = sh2 * f2 + sh3 * f3
	zm    = zmol + znl * t
	if inited {
		zm = zmol
	}
	zf    = zm + 2.0 * zel * sin(zm)
	sinzf = sin(zf)
	f2    =  0.5 * sinzf * sinzf - 0.25
	f3    = -0.5 * sinzf * cos(zf)
	sel   = ee2 * f2 + e3 * f3
	sil   = xi2 * f2 + xi3 * f3
	sll   = xl2 * f2 + xl3 * f3 + xl4 * sinzf
	sghl  = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf
	shll  = xh2 * f2 + xh3 * f3
	pe    = ses + sel
	pinc  = sis + sil
	pl    = sls + sll
	pgh   = sghs + sghl
	ph    = shs + shll

	if inited{
		pe    = pe - peo
		pinc  = pinc - pinco
		pl    = pl - plo
		pgh   = pgh - pgho
		ph    = ph - pho
		out_inclp = out_inclp + pinc
		out_ep    = out_ep + pe
		sinip = sin(out_inclp)
		cosip = cos(out_inclp)

		/* ----------------- apply periodics directly ------------ */
		//  sgp4fix for lyddane choice
		//  strn3 used original inclination - this is technically feasible
		//  gsfc used perturbed inclination - also technically feasible
		//  probably best to readjust the 0.2 limit value and limit discontinuity
		//  0.2 rad = 11.45916 deg
		//  use next line for original strn3 approach and original inclination
		//  if (inclo >= 0.2)
		//  use next line for gsfc version and perturbed inclination
		if out_inclp >= 0.2 {
			ph     = ph / sinip
			pgh    = pgh - cosip * ph
			out_argpp  = out_argpp + pgh
			out_nodep  = out_nodep + ph
			out_mp     = out_mp + pl
		}
		else
		{
			/* ---- apply periodics with lyddane modification ---- */
			sinop  = sin(out_nodep)
			cosop  = cos(out_nodep)
			alfdp  = sinip * sinop
			betdp  = sinip * cosop
			dalf   =  ph * cosop + pinc * cosip * sinop
			dbet   = -ph * sinop + pinc * cosip * cosop
			alfdp  = alfdp + dalf
			betdp  = betdp + dbet
			out_nodep  = fmod(out_nodep, Double.twoπ)
			//  sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			if ((out_nodep < 0.0) && (opsmode == "a")) {
				out_nodep = out_nodep + Double.twoπ
			}
			xls    = (out_mp) + (out_argpp) + cosip * (out_nodep)
			dls    = pl + pgh - pinc * (out_nodep) * sinip
			xls    = xls + dls
			xnoh   = out_nodep
			out_nodep  = ArcTan(alfdp, betdp)
			//  sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			if (out_nodep < 0.0) && (opsmode == "a") {
				out_nodep = out_nodep + Double.twoπ
			}
			if fabs(xnoh - out_nodep) > π {
				if (out_nodep < xnoh) {
					out_nodep = out_nodep + Double.twoπ
				} else {
					out_nodep = out_nodep - Double.twoπ
				}
			}
			out_mp    = out_mp + pl
			out_argpp = xls - out_mp - cosip * (out_nodep)
		}
	}   // if init == 'n'
}  // end dpper

/*-----------------------------------------------------------------------------
*
*                           procedure dscom
*
*  this procedure provides deep space common items used by both the secular
*    and periodics subroutines.  input is provided as shown. this routine
*    used to be called dpper, but the functions inside weren't well organized.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    epoch       -
*    ep          - eccentricity
*    argpp       - argument of perigee
*    tc          -
*    inclp       - inclination
*    nodep       - right ascension of ascending node
*    np          - mean motion
*
*  outputs       :
*    sinim  , cosim  , sinomm , cosomm , snodm  , cnodm
*    day         -
*    e3          -
*    ee2         -
*    em          - eccentricity
*    emsq        - eccentricity squared
*    gam         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    rtemsq      -
*    se2, se3         -
*    sgh2, sgh3, sgh4        -
*    sh2, sh3, si2, si3, sl2, sl3, sl4         -
*    s1, s2, s3, s4, s5, s6, s7          -
*    ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3         -
*    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
*    xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4         -
*    nm          - mean motion
*    z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33         -
*    zmol        -
*    zmos        -
*
*  locals        :
*    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10         -
*    betasq      -
*    cc          -
*    ctem, stem        -
*    x1, x2, x3, x4, x5, x6, x7, x8          -
*    xnodce      -
*    xnoi        -
*    zcosg  , zsing  , zcosgl , zsingl , zcosh  , zsinh  , zcoshl , zsinhl ,
*    zcosi  , zsini  , zcosil , zsinil ,
*    zx          -
*    zy          -
*
*  coupling      :
*    none.
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/

func dscom
(
	_ epoch : Double,  _ ep : Double,     _ argpp : Double,   _ tc : Double,     _ inclp : Double,
_ nodep : Double,  _ np : Double,
_ out_snodm : inout Double, _ out_cnodm : inout Double, _ out_sinim : inout Double,  _ out_cosim : inout Double, _ out_sinomm : inout Double,
_ out_cosomm : inout Double,_ out_day : inout Double,   _ out_e3 : inout Double,     _ out_ee2 : inout Double,   _ out_em : inout Double,
_ out_emsq : inout Double,  _ out_gam : inout Double,   _ out_peo : inout Double,    _ out_pgho : inout Double,  _ out_pho : inout Double,
_ out_pinco : inout Double, _ out_plo : inout Double,   _ out_rtemsq : inout Double, _ out_se2 : inout Double,   _ out_se3 : inout Double,
_ out_sgh2 : inout Double,  _ out_sgh3 : inout Double,  _ out_sgh4 : inout Double,   _ out_sh2 : inout Double,   _ out_sh3 : inout Double,
_ out_si2 : inout Double,   _ out_si3 : inout Double,   _ out_sl2 : inout Double,    _ out_sl3 : inout Double,   _ out_sl4 : inout Double,
_ out_s1 : inout Double,    _ out_s2 : inout Double,    _ out_s3 : inout Double,     _ out_s4 : inout Double,    _ out_s5 : inout Double,
_ out_s6 : inout Double,    _ out_s7 : inout Double,    _ out_ss1 : inout Double,    _ out_ss2 : inout Double,   _ out_ss3 : inout Double,
_ out_ss4 : inout Double,   _ out_ss5 : inout Double,   _ out_ss6 : inout Double,    _ out_ss7 : inout Double,   _ out_sz1 : inout Double,
_ out_sz2 : inout Double,   _ out_sz3 : inout Double,   _ out_sz11 : inout Double,   _ out_sz12 : inout Double,  _ out_sz13 : inout Double,
_ out_sz21 : inout Double,  _ out_sz22 : inout Double,  _ out_sz23 : inout Double,   _ out_sz31 : inout Double,  _ out_sz32 : inout Double,
_ out_sz33 : inout Double,  _ out_xgh2 : inout Double,  _ out_xgh3 : inout Double,   _ out_xgh4 : inout Double,  _ out_xh2 : inout Double,
_ out_xh3 : inout Double,   _ out_xi2 : inout Double,   _ out_xi3 : inout Double,    _ out_xl2 : inout Double,   _ out_xl3 : inout Double,
_ out_xl4 : inout Double,   _ out_nm : inout Double,    _ out_z1 : inout Double,     _ out_z2 : inout Double,    _ out_z3 : inout Double,
_ out_z11 : inout Double,   _ out_z12 : inout Double,   _ out_z13 : inout Double,    _ out_z21 : inout Double,   _ out_z22 : inout Double,
_ out_z23 : inout Double,   _ out_z31 : inout Double,   _ out_z32 : inout Double,    _ out_z33 : inout Double,   _ out_zmol : inout Double,
_ out_zmos : inout Double
)
{
	/* -------------------------- constants ------------------------- */
	let zes     =  0.01675
	let zel     =  0.05490
	let c1ss    =  2.9864797e-6
	let c1l     =  4.7968065e-7
	let zsinis  =  0.39785416
	let zcosis  =  0.91744867
	let zcosgs  =  0.1945905
	let zsings  = -0.98088458

	/* --------------------- local variables ------------------------ */
	var a1    , a2    , a3    , a4    , a5    , a6    , a7    ,
	a8    , a9    , a10   , betasq, cc    , ctem  , stem  ,
	x1    , x2    , x3    , x4    , x5    , x6    , x7    ,
	x8    , xnodce, xnoi  , zcosg , zcosgl, zcosh , zcoshl,
	zcosi , zcosil, zsing , zsingl, zsinh , zsinhl, zsini ,
	zsinil, zx    , zy : Double

	out_nm     = np
	out_em     = ep
	out_snodm  = sin(nodep)
	out_cnodm  = cos(nodep)
	out_sinomm = sin(argpp)
	out_cosomm = cos(argpp)
	out_sinim  = sin(inclp)
	out_cosim  = cos(inclp)
	out_emsq   = (out_em) * (out_em)
	betasq = 1.0 - out_emsq
	out_rtemsq = √(betasq)

	/* ----------------- initialize lunar solar terms --------------- */
	out_peo    = 0.0
	out_pinco  = 0.0
	out_plo    = 0.0
	out_pgho   = 0.0
	out_pho    = 0.0
	out_day    = epoch + 18261.5 + tc / 1440.0
	xnodce = fmod(4.5236020 - 9.2422029e-4 * (out_day), Double.twoπ)
	stem   = sin(xnodce)
	ctem   = cos(xnodce)
	zcosil = 0.91375164 - 0.03568096 * ctem
	zsinil = √(1.0 - zcosil * zcosil)
	zsinhl = 0.089683511 * stem / zsinil
	zcoshl = √(1.0 - zsinhl * zsinhl)
	out_gam    = 5.8351514 + 0.0019443680 * (out_day)
	zx     = 0.39785416 * stem / zsinil
	zy     = zcoshl * ctem + 0.91744867 * zsinhl * stem
	zx     = ArcTan(zx, zy)
	zx     = out_gam + zx - xnodce
	zcosgl = cos(zx)
	zsingl = sin(zx)

	/* ------------------------- do solar terms --------------------- */
	zcosg = zcosgs
	zsing = zsings
	zcosi = zcosis
	zsini = zsinis
	zcosh = out_cnodm
	zsinh = out_snodm
	cc    = c1ss
	xnoi  = 1.0 / out_nm

	for lsflg in 1...2 {
		a1  =   zcosg * zcosh + zsing * zcosi * zsinh
		a3  =  -zsing * zcosh + zcosg * zcosi * zsinh
		a7  =  -zcosg * zsinh + zsing * zcosi * zcosh
		a8  =   zsing * zsini
		a9  =   zsing * zsinh + zcosg * zcosi * zcosh
		a10 =   zcosg * zsini
		a2  =   out_cosim * a7 + out_sinim * a8
		a4  =   out_cosim * a9 + out_sinim * a10
		a5  =  -out_sinim * a7 + out_cosim * a8
		a6  =  -out_sinim * a9 + out_cosim * a10

		x1  =  a1 * (out_cosomm) + a2 * (out_sinomm)
		x2  =  a3 * (out_cosomm) + a4 * (out_sinomm)
		x3  = -a1 * (out_sinomm) + a2 * (out_cosomm)
		x4  = -a3 * (out_sinomm) + a4 * (out_cosomm)
		x5  =  a5 * (out_sinomm)
		x6  =  a6 * (out_sinomm)
		x7  =  a5 * (out_cosomm)
		x8  =  a6 * (out_cosomm)

		out_z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3
		out_z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4
		out_z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4
		out_z1  =  3.0 *  (a1 * a1 + a2 * a2) + out_z31 * (out_emsq)
		out_z2  =  6.0 *  (a1 * a3 + a2 * a4) + out_z32 * (out_emsq)
		out_z3  =  3.0 *  (a3 * a3 + a4 * a4) + out_z33 * (out_emsq)
		out_z11 = -6.0 * a1 * a5 + out_emsq *  (-24.0 * x1 * x7-6.0 * x3 * x5)
		out_z12 = -6.0 *  (a1 * a6 + a3 * a5) + out_emsq *
			(-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5))
		out_z13 = -6.0 * a3 * a6 + out_emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6)
		out_z21 =  6.0 * a2 * a5 + out_emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7)
		out_z22 =  6.0 *  (a4 * a5 + a2 * a6) + out_emsq *
			(24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8))
		out_z23 =  6.0 * a4 * a6 + out_emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8)
		out_z1  = out_z1 + out_z1 + betasq * (out_z31)
		out_z2  = out_z2 + out_z2 + betasq * (out_z32)
		out_z3  = out_z3 + out_z3 + betasq * (out_z33)
		out_s3  = cc * xnoi
		out_s2  = -0.5 * (out_s3) / (out_rtemsq)
		out_s4  = out_s3 * (out_rtemsq)
		out_s1  = -15.0 * (out_em) * (out_s4)
		out_s5  = x1 * x3 + x2 * x4
		out_s6  = x2 * x3 + x1 * x4
		out_s7  = x2 * x4 - x1 * x3

		/* ----------------------- do lunar terms ------------------- */
		if (lsflg == 1) {
			out_ss1   = out_s1
			out_ss2   = out_s2
			out_ss3   = out_s3
			out_ss4   = out_s4
			out_ss5   = out_s5
			out_ss6   = out_s6
			out_ss7   = out_s7
			out_sz1   = out_z1
			out_sz2   = out_z2
			out_sz3   = out_z3
			out_sz11  = out_z11
			out_sz12  = out_z12
			out_sz13  = out_z13
			out_sz21  = out_z21
			out_sz22  = out_z22
			out_sz23  = out_z23
			out_sz31  = out_z31
			out_sz32  = out_z32
			out_sz33  = out_z33
			zcosg = zcosgl
			zsing = zsingl
			zcosi = zcosil
			zsini = zsinil
			zcosh = zcoshl * (out_cnodm) + zsinhl * (out_snodm)
			zsinh = (out_snodm) * zcoshl - (out_cnodm) * zsinhl
			cc    = c1l
		}
	}

	out_zmol = fmod(4.7199672 + 0.22997150  * (out_day) - out_gam, Double.twoπ)
	out_zmos = fmod(6.2565837 + 0.017201977 * (out_day), Double.twoπ)

	/* ------------------------ do solar terms ---------------------- */
	out_se2  =   2.0 * out_ss1 * out_ss6
	out_se3  =   2.0 * out_ss1 * out_ss7
	out_si2  =   2.0 * out_ss2 * out_sz12
	out_si3  =   2.0 * out_ss2 * (out_sz13 - out_sz11)
	out_sl2  =  -2.0 * out_ss3 * out_sz2
	out_sl3  =  -2.0 * out_ss3 * (out_sz3 - out_sz1)
	out_sl4  =  -2.0 * out_ss3 * (-21.0 - 9.0 * out_emsq) * zes
	out_sgh2 =   2.0 * out_ss4 * out_sz32
	out_sgh3 =   2.0 * out_ss4 * (out_sz33 - out_sz31)
	out_sgh4 = -18.0 * out_ss4 * zes
	out_sh2  =  -2.0 * out_ss2 * out_sz22
	out_sh3  =  -2.0 * out_ss2 * (out_sz23 - out_sz21)

	/* ------------------------ do lunar terms ---------------------- */
	out_ee2  =   2.0 * out_s1 * out_s6
	out_e3   =   2.0 * out_s1 * out_s7
	out_xi2  =   2.0 * out_s2 * out_z12
	out_xi3  =   2.0 * out_s2 * (out_z13 - out_z11)
	out_xl2  =  -2.0 * out_s3 * out_z2
	out_xl3  =  -2.0 * out_s3 * (out_z3 - out_z1)
	out_xl4  =  -2.0 * out_s3 * (-21.0 - 9.0 * out_emsq) * zel
	out_xgh2 =   2.0 * out_s4 * out_z32
	out_xgh3 =   2.0 * out_s4 * (out_z33 - out_z31)
	out_xgh4 = -18.0 * out_s4 * zel
	out_xh2  =  -2.0 * out_s2 * out_z22
	out_xh3  =  -2.0 * out_s2 * (out_z23 - out_z21)

	//#include "debug2.cpp"
}  // end dscom

/*-----------------------------------------------------------------------------
*
*                           procedure dsinit
*
*  this procedure provides deep space contributions to mean motion dot due
*    to geopotential resonance with half day and one day orbits.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    cosim, sinim-
*    emsq        - eccentricity squared
*    argpo       - argument of perigee
*    s1, s2, s3, s4, s5      -
*    ss1, ss2, ss3, ss4, ss5 -
*    sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33 -
*    t           - time
*    tc          -
*    gsto        - greenwich sidereal time                   rad
*    mo          - mean anomaly
*    mdot        - mean anomaly dot (rate)
*    no          - mean motion
*    nodeo       - right ascension of ascending node
*    nodedot     - right ascension of ascending node dot (rate)
*    xpidot      -
*    z1, z3, z11, z13, z21, z23, z31, z33 -
*    eccm        - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    mm          - mean anomaly
*    xn          - mean motion
*    nodem       - right ascension of ascending node
*
*  outputs       :
*    em          - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    nodem       - right ascension of ascending node
*    irez        - flag for resonance           0-none, 1-one day, 2-half day
*    atime       -
*    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433    -
*    dedt        -
*    didt        -
*    dmdt        -
*    dndt        -
*    dnodt       -
*    domdt       -
*    del1, del2, del3        -
*    ses  , sghl , sghs , sgs  , shl  , shs  , sis  , sls
*    theta       -
*    xfact       -
*    xlamo       -
*    xli         -
*    xni
*
*  locals        :
*    ainv2       -
*    aonv        -
*    cosisq      -
*    eoc         -
*    f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543  -
*    g200, g201, g211, g300, g310, g322, g410, g422, g520, g521, g532, g533  -
*    sini2       -
*    temp        -
*    temp1       -
*    theta       -
*    xno2        -
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/

func dsinit(_ whichconst : GravityConstantsType,
_ cosim: Double,  _ emsq: inout Double,   _ argpo: Double,   _ s1: Double,     _ s2: Double,
_ s3: Double,     _ s4: Double,     _ s5: Double,      _ sinim: Double,  _ ss1: Double,
_ ss2: Double,    _ ss3: Double,    _ ss4: Double,     _ ss5: Double,    _ sz1: Double,
_ sz3: Double,    _ sz11: Double,   _ sz13: Double,    _ sz21: Double,   _ sz23: Double,
_ sz31: Double,   _ sz33: Double,   _ t: Double,       _ tc: Double,     _ gsto: Double,
_ mo: Double,     _ mdot: Double,   _ no: Double,      _ nodeo: Double,  _ nodedot: Double,
_ xpidot: Double, _ z1: Double,     _ z3: Double,      _ z11: Double,    _ z13: Double,
_ z21: Double,    _ z23: Double,    _ z31: Double,     _ z33: Double,    _ ecco: Double,
_ eccsq: Double,  _ out_em: inout Double,    _ out_argpm: inout Double,  _ out_inclm: inout Double, _ out_mm: inout Double,
_ out_nm: inout Double,    _ out_nodem: inout Double,
_ out_irez : inout Int,
_ out_atime: inout Double, _ out_d2201: inout Double, _ out_d2211: inout Double,  _ out_d3210: inout Double, _ out_d3222: inout Double,
_ out_d4410: inout Double, _ out_d4422: inout Double, _ out_d5220: inout Double,  _ out_d5232: inout Double, _ out_d5421: inout Double,
_ out_d5433: inout Double, _ out_dedt: inout Double,  _ out_didt: inout Double,   _ out_dmdt: inout Double,  _ out_dndt: inout Double,
_ out_dnodt: inout Double, _ out_domdt: inout Double, _ out_del1: inout Double,   _ out_del2: inout Double,  _ out_del3: inout Double,
	_ out_xfact: inout Double, _ out_xlamo: inout Double, _ out_xli: inout Double,    _ out_xni : inout Double
)
{
	/* --------------------- local variables ------------------------ */

	var ainv2 , cosisq, eoc, f220 , f221  , f311  ,
	f321  , f322  , f330  , f441  , f442  , f522  , f523  ,
	f542  , f543  , g200  , g201  , g211  , g300  , g310  ,
	g322  , g410  , g422  , g520  , g521  , g532  , g533  ,
	ses   , sgs   , sghl  , sghs  , shs   , shll  , sis   ,
	sini2 , sls   , temp  , temp1 , theta , xno2  , q22   ,
	q31   , q33   , root22, root44, root54, rptim , root32,
	root52, x2o3  , znl   , emo   , zns   , emsqo : Double
	var aonv = 0.0


	q22    = 1.7891679e-6
	q31    = 2.1460748e-6
	q33    = 2.2123015e-7
	root22 = 1.7891679e-6
	root44 = 7.3636953e-9
	root54 = 2.1765803e-9
	rptim  = 4.37526908801129966e-3 // this equates to 7.29211514668855e-5 rad/sec
	root32 = 3.7393792e-7
	root52 = 1.1428639e-7
	x2o3   = 2.0 / 3.0
	znl    = 1.5835218e-4
	zns    = 1.19459e-5

	// sgp4fix identify constants and allow alternate values
	let gravityConstants = GravityConstants(type: whichconst)

	/* -------------------- deep space initialization ------------ */
	out_irez = 0
	if (out_nm < 0.0052359877) && (out_nm > 0.0034906585) {
		out_irez = 1
	}
	if (out_nm >= 8.26e-3) && (out_nm <= 9.24e-3) && (out_em >= 0.5) {
		out_irez = 2
	}

	/* ------------------------ do solar terms ------------------- */
	ses  =  ss1 * zns * ss5
	sis  =  ss2 * zns * (sz11 + sz13)
	sls  = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq)
	sghs =  ss4 * zns * (sz31 + sz33 - 6.0)
	shs  = -zns * ss2 * (sz21 + sz23)
	// sgp4fix for 180 deg incl
	if (out_inclm < 5.2359877e-2) || (out_inclm > π - 5.2359877e-2) {
		shs = 0.0
	}
	if sinim != 0.0 {
		shs = shs / sinim
	}
	sgs  = sghs - cosim * shs

	/* ------------------------- do lunar terms ------------------ */
	out_dedt = ses + s1 * znl * s5
	out_didt = sis + s2 * znl * (z11 + z13)
	out_dmdt = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq)
	sghl = s4 * znl * (z31 + z33 - 6.0)
	shll = -znl * s2 * (z21 + z23)
	// sgp4fix for 180 deg incl
	if (out_inclm < 5.2359877e-2) || (out_inclm > π - 5.2359877e-2) {
		shll = 0.0
	}
	out_domdt = sgs + sghl
	out_dnodt = shs
	if sinim != 0.0 {
		out_domdt = out_domdt - cosim / sinim * shll
		out_dnodt = out_dnodt + shll / sinim
	}

	/* ----------- calculate deep space resonance effects -------- */
	out_dndt   = 0.0
	theta  = fmod(gsto + tc * rptim, Double.twoπ)
	out_em     = out_em + out_dedt * t
	out_inclm  = out_inclm + out_didt * t
	out_argpm  = out_argpm + out_domdt * t
	out_nodem  = out_nodem + out_dnodt * t
	out_mm     = out_mm + out_dmdt * t
	//   sgp4fix for negative inclinations
	//   the following if statement should be commented out
	//if (inclm < 0.0)
	//  {
	//    inclm  = -inclm
	//    argpm  = argpm - pi
	//    nodem = nodem + pi
	//  }

	/* -------------- initialize the resonance terms ------------- */
	if out_irez != 0 {
		aonv = pow(out_nm / gravityConstants.xke, x2o3)

		/* ---------- geopotential resonance for 12 hour orbits ------ */
		if out_irez == 2 {
			cosisq = cosim * cosim
			emo    = out_em
			out_em     = ecco
			emsqo  = emsq
			emsq   = eccsq
			eoc    = out_em * emsq
			g201   = -0.306 - (out_em - 0.64) * 0.440

			if out_em <= 0.65 {
				g211 =    3.616  -  13.2470 * out_em +  16.2900 * emsq
				g310 =  -19.302  + 117.3900 * out_em - 228.4190 * emsq +  156.5910 * eoc
				g322 =  -18.9068 + 109.7927 * out_em - 214.6334 * emsq +  146.5816 * eoc
				g410 =  -41.122  + 242.6940 * out_em - 471.0940 * emsq +  313.9530 * eoc
				g422 = -146.407  + 841.8800 * out_em - 1629.014 * emsq + 1083.4350 * eoc
				g520 = -532.114  + 3017.977 * out_em - 5740.032 * emsq + 3708.2760 * eoc
			} else {
				g211 =   -72.099 +   331.819 * out_em -   508.738 * emsq +   266.724 * eoc
				g310 =  -346.844 +  1582.851 * out_em -  2415.925 * emsq +  1246.113 * eoc
				g322 =  -342.585 +  1554.908 * out_em -  2366.899 * emsq +  1215.972 * eoc
				g410 = -1052.797 +  4758.686 * out_em -  7193.992 * emsq +  3651.957 * eoc
				g422 = -3581.690 + 16178.110 * out_em - 24462.770 * emsq + 12422.520 * eoc
				if out_em > 0.715 {
					g520 = -5149.66 + 29936.92 * out_em - 54087.36 * emsq + 31324.56 * eoc
				} else {
					g520 = 1464.74 -  4664.75 * out_em +  3763.64 * emsq
				}
			}
			if out_em < 0.7 {
				g533 = -919.22770 + 4988.6100 * out_em - 9064.7700 * emsq + 5542.21  * eoc
				g521 = -822.71072 + 4568.6173 * out_em - 8491.4146 * emsq + 5337.524 * eoc
				g532 = -853.66600 + 4690.2500 * out_em - 8624.7700 * emsq + 5341.4  * eoc
			} else {
				g533 = -37995.780 + 161616.52 * out_em - 229838.20 * emsq + 109377.94 * eoc
				g521 = -51752.104 + 218913.95 * out_em - 309468.16 * emsq + 146349.42 * eoc
				g532 = -40023.880 + 170470.89 * out_em - 242699.48 * emsq + 115605.82 * eoc
			}

			sini2 =  sinim * sinim
			f220 =  0.75 * (1.0 + 2.0 * cosim+cosisq)
			f221 =  1.5 * sini2
			f321 =  1.875 * sinim  *  (1.0 - 2.0 * cosim - 3.0 * cosisq)
			f322 = -1.875 * sinim  *  (1.0 + 2.0 * cosim - 3.0 * cosisq)
			f441 = 35.0 * sini2 * f220
			f442 = 39.3750 * sini2 * sini2
			f522 =  9.84375 * sinim * (sini2 * (1.0 - 2.0 * cosim - 5.0 * cosisq) +
				0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq) )
			f523 = sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim +
				10.0 * cosisq) + 6.56250012 * (1.0+2.0 * cosim - 3.0 * cosisq))
			f542 = 29.53125 * sinim * (2.0 - 8.0 * cosim+cosisq *
				(-12.0 + 8.0 * cosim + 10.0 * cosisq))
			f543 = 29.53125 * sinim * (-2.0 - 8.0 * cosim+cosisq *
				(12.0 + 8.0 * cosim - 10.0 * cosisq))
			xno2  =  out_nm * out_nm
			ainv2 =  aonv * aonv
			temp1 =  3.0 * xno2 * ainv2
			temp  =  temp1 * root22
			out_d2201 =  temp * f220 * g201
			out_d2211 =  temp * f221 * g211
			temp1 =  temp1 * aonv
			temp  =  temp1 * root32
			out_d3210 =  temp * f321 * g310
			out_d3222 =  temp * f322 * g322
			temp1 =  temp1 * aonv
			temp  =  2.0 * temp1 * root44
			out_d4410 =  temp * f441 * g410
			out_d4422 =  temp * f442 * g422
			temp1 =  temp1 * aonv
			temp  =  temp1 * root52
			out_d5220 =  temp * f522 * g520
			out_d5232 =  temp * f523 * g532
			temp  =  2.0 * temp1 * root54
			out_d5421 =  temp * f542 * g521
			out_d5433 =  temp * f543 * g533
			out_xlamo =  fmod(mo + nodeo + nodeo-theta - theta, Double.twoπ)
			out_xfact =  mdot + out_dmdt + 2.0 * (nodedot + out_dnodt - rptim) - no
			out_em    = emo
			emsq = emsqo
		}

		/* ---------------- synchronous resonance terms -------------- */
		if out_irez == 1 {
			g200  = 1.0 + emsq * (-2.5 + 0.8125 * emsq)
			g310  = 1.0 + 2.0 * emsq
			g300  = 1.0 + emsq * (-6.0 + 6.60937 * emsq)
			f220  = 0.75 * (1.0 + cosim) * (1.0 + cosim)
			f311  = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim)
			f330  = 1.0 + cosim
			f330  = 1.875 * f330 * f330 * f330
			out_del1  = 3.0 * out_nm * out_nm * aonv * aonv
			out_del2  = 2.0 * out_del1 * f220 * g200 * q22
			out_del3  = 3.0 * out_del1 * f330 * g300 * q33 * aonv
			out_del1  = out_del1 * f311 * g310 * q31 * aonv
			out_xlamo = fmod(mo + nodeo + argpo - theta, Double.twoπ)
			out_xfact = mdot + xpidot - rptim + out_dmdt + out_domdt + out_dnodt - no
		}

		/* ------------ for sgp4, initialize the integrator ---------- */
		out_xli   = out_xlamo
		out_xni   = no
		out_atime = 0.0
		out_nm    = no + out_dndt
	}

	//#include "debug3.cpp"
}  // end dsinit

/*-----------------------------------------------------------------------------
*
*                           procedure dspace
*
*  this procedure provides deep space contributions to mean elements for
*    perturbing third body.  these effects have been averaged over one
*    revolution of the sun and moon.  for earth resonance effects, the
*    effects have been averaged over no revolutions of the satellite.
*    (mean motion)
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433 -
*    dedt        -
*    del1, del2, del3  -
*    didt        -
*    dmdt        -
*    dnodt       -
*    domdt       -
*    irez        - flag for resonance           0-none, 1-one day, 2-half day
*    argpo       - argument of perigee
*    argpdot     - argument of perigee dot (rate)
*    t           - time
*    tc          -
*    gsto        - gst
*    xfact       -
*    xlamo       -
*    no          - mean motion
*    atime       -
*    em          - eccentricity
*    ft          -
*    argpm       - argument of perigee
*    inclm       - inclination
*    xli         -
*    mm          - mean anomaly
*    xni         - mean motion
*    nodem       - right ascension of ascending node
*
*  outputs       :
*    atime       -
*    em          - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    xli         -
*    mm          - mean anomaly
*    xni         -
*    nodem       - right ascension of ascending node
*    dndt        -
*    nm          - mean motion
*
*  locals        :
*    delt        -
*    ft          -
*    theta       -
*    x2li        -
*    x2omi       -
*    xl          -
*    xldot       -
*    xnddt       -
*    xndt        -
*    xomi        -
*
*  coupling      :
*    none        -
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/

func dspace(_ irez : Int,
			_ d2201 : Double,  _ d2211 : Double,  _ d3210 : Double,   _ d3222 : Double,  _ d4410 : Double,
			_ d4422 : Double,  _ d5220 : Double,  _ d5232 : Double,   _ d5421 : Double,  _ d5433 : Double,
			_ dedt : Double,   _ del1 : Double,   _ del2 : Double,    _ del3 : Double,   _ didt : Double,
			_ dmdt : Double,   _ dnodt : Double,  _ domdt : Double,   _ argpo : Double,  _ argpdot : Double,
			_ t : Double,      _ tc : Double,     _ gsto : Double,    _ xfact : Double,  _ xlamo : Double,
			_ no : Double,
			_ out_atime : inout Double, _ out_em : inout Double,    _ out_argpm : inout Double,  _ out_inclm : inout Double, _ out_xli : inout Double,
			_ out_mm : inout Double,    _ out_xni : inout Double,   _ out_nodem : inout Double,  _ out_dndt : inout Double,  _ out_nm : inout Double
	) {
	var iretn /* , iret // analyse pas utilise */ : Int
	var delt, ft, theta, x2li, x2omi, xl, xomi, g22, g32,
	g44, g52, g54, fasx2, fasx4, fasx6, rptim , step2, stepn , stepp : Double
	var xldot = 0.0
	var xnddt = 0.0
	var xndt = 0.0

	fasx2 = 0.13130908
	fasx4 = 2.8843198
	fasx6 = 0.37448087
	g22   = 5.7686396
	g32   = 0.95240898
	g44   = 1.8014998
	g52   = 1.0508330
	g54   = 4.4108898
	rptim = 4.37526908801129966e-3 // this equates to 7.29211514668855e-5 rad/sec
	stepp =    720.0
	stepn =   -720.0
	step2 = 259200.0

	/* ----------- calculate deep space resonance effects ----------- */
	out_dndt   = 0.0
	theta  = fmod(gsto + tc * rptim, Double.twoπ)
	out_em     = out_em + dedt * t

	out_inclm  = out_inclm + didt * t
	out_argpm  = out_argpm + domdt * t
	out_nodem  = out_nodem + dnodt * t
	out_mm     = out_mm + dmdt * t

	//   sgp4fix for negative inclinations
	//   the following if statement should be commented out
	//  if (inclm < 0.0)
	// {
	//    inclm = -inclm
	//    argpm = argpm - pi
	//    nodem = nodem + pi
	//  }

	/* - update resonances : numerical (euler-maclaurin) integration - */
	/* ------------------------- epoch restart ----------------------  */
	//   sgp4fix for propagator problems
	//   the following integration works for negative time steps and periods
	//   the specific changes are unknown because the original code was so convoluted

	// sgp4fix take out atime = 0.0 and fix for faster operation
	ft    = 0.0
	if (irez != 0)
	{
		// sgp4fix streamline check
		if ((out_atime == 0.0) || (t * out_atime <= 0.0) || (fabs(t) < fabs(out_atime)) )
		{
			out_atime  = 0.0
			out_xni    = no
			out_xli    = xlamo
		}
		// sgp4fix move check outside loop
		if t > 0.0 {
			delt = stepp
		} else {
			delt = stepn
		}

		iretn = 381 // added for do loop
		// iret  =   0 // added for loop // analyse pas utilise
		while (iretn == 381)
		{
			/* ------------------- dot terms calculated ------------- */
			/* ----------- near - synchronous resonance terms ------- */
			if (irez != 2)
			{
				xndt  = del1 * sin(out_xli - fasx2) + del2 * sin(2.0 * (out_xli - fasx4)) +
					del3 * sin(3.0 * (out_xli - fasx6))
				xldot = out_xni + xfact
				xnddt = del1 * cos(out_xli - fasx2) +
					2.0 * del2 * cos(2.0 * (out_xli - fasx4)) +
					3.0 * del3 * cos(3.0 * (out_xli - fasx6))
				xnddt = xnddt * xldot
			}
			else
			{
				/* --------- near - half-day resonance terms -------- */
				xomi  = argpo + argpdot * out_atime
				x2omi = xomi + xomi
				x2li  = out_xli + out_xli
				xndt  = d2201 * sin(x2omi + out_xli - g22) + d2211 * sin(out_xli - g22) +
					d3210 * sin(xomi + out_xli - g32)  + d3222 * sin(-xomi + out_xli - g32) +
				d4410 * sin(x2omi + x2li - g44) + d4422 * sin(x2li - g44) +
					d5220 * sin(xomi + out_xli - g52)  + d5232 * sin(-xomi + out_xli - g52) +
				d5421 * sin(xomi + x2li - g54) + d5433 * sin(-xomi + x2li - g54)
				xldot = out_xni + xfact
				xnddt = d2201 * cos(x2omi + out_xli - g22) + d2211 * cos(out_xli - g22) +
					d3210 * cos(xomi + out_xli - g32) + d3222 * cos(-xomi + out_xli - g32) +
					d5220 * cos(xomi + out_xli - g52) + d5232 * cos(-xomi + out_xli - g52) +
					2.0 * (d4410 * cos(x2omi + x2li - g44) +
						d4422 * cos(x2li - g44) + d5421 * cos(xomi + x2li - g54) +
						d5433 * cos(-xomi + x2li - g54))
				xnddt = xnddt * xldot
			}

			/* ----------------------- integrator ------------------- */
			// sgp4fix move end checks to end of routine
			if (fabs(t - out_atime) >= stepp)
			{
				// iret  = 0 // analyse pas utilise
				iretn = 381
			}
			else // exit here
			{
				ft    = t - out_atime
				iretn = 0
			}

			if (iretn == 381)
			{
				out_xli   = out_xli + xldot * delt + xndt * step2
				out_xni   = out_xni + xndt * delt + xnddt * step2
				out_atime = out_atime + delt
			}
		}  // while iretn = 381

		out_nm = out_xni + xndt * ft + xnddt * ft * ft * 0.5
		xl = out_xli + xldot * ft + xndt * ft * ft * 0.5
		if (irez != 1)
		{
			out_mm   = xl - 2.0 * out_nodem + 2.0 * theta
			out_dndt = out_nm - no
		}
		else
		{
			out_mm   = xl - out_nodem - out_argpm + theta
			out_dndt = out_nm - no
		}
		out_nm = no + out_dndt
	}

	//#include "debug4.cpp"
}  // end dsspace

/*-----------------------------------------------------------------------------
*
*                           procedure initl
*
*  this procedure initializes the spg4 propagator. all the initialization is
*    consolidated here instead of having multiple loops inside other routines.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    ecco        - eccentricity                           0.0 - 1.0
*    epoch       - epoch time in days from jan 0, 1950. 0 hr
*    inclo       - inclination of satellite
*    no          - mean motion of satellite
*    satn        - satellite number
*
*  outputs       :
*    ainv        - 1.0 / a
*    ao          - semi major axis
*    con41       -
*    con42       - 1.0 - 5.0 cos(i)
*    cosio       - cosine of inclination
*    cosio2      - cosio squared
*    eccsq       - eccentricity squared
*    method      - flag for deep space                    'd', 'n'
*    omeosq      - 1.0 - ecco * ecco
*    posq        - semi-parameter squared
*    rp          - radius of perigee
*    rteosq      - square root of (1.0 - ecco*ecco)
*    sinio       - sine of inclination
*    gsto        - gst at time of observation               rad
*    no          - mean motion of satellite
*
*  locals        :
*    ak          -
*    d1          -
*    del         -
*    adel        -
*    po          -
*
*  coupling      :
*    gstime      - find greenwich sidereal time from the julian date
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/

func initl(_ satn : Int, _ whichconst : GravityConstantsType,
		   _ ecco : Double,   _ epoch : Double,  _ inclo : Double,   _ out_no : inout Double,
		   _ out_method : inout Character,
		   _ out_ainv : inout Double,  _ out_ao : inout Double,    _ out_con41 : inout Double,  _ out_con42 : inout Double, _ out_cosio : inout Double,
		   _ out_cosio2 : inout Double,_ out_eccsq : inout Double, _ out_omeosq : inout Double, _ out_posq : inout Double,
		   _ out_rp : inout Double,    _ out_rteosq : inout Double,_ out_sinio  : inout Double, _ out_gsto : inout Double,
		   _ opsmode : Character) {
	/* --------------------- local variables ------------------------ */
	var ak, d1, del, adel, po, x2o3  : Double

	// sgp4fix use old way of finding gst
	var ds70 : Double
	var ts70, tfrac, c1, thgr70, fk5r, c1p2p /*, thgr*/ /*, thgro*/ : Double

	/* ----------------------- earth constants ---------------------- */
	// sgp4fix identify constants and allow alternate values
	let gravityConstants = GravityConstants(type: whichconst)
	x2o3   = 2.0 / 3.0

	/* ------------- calculate auxillary epoch quantities ---------- */
	out_eccsq  = ecco * ecco
	out_omeosq = 1.0 - out_eccsq
	out_rteosq = √(out_omeosq)
	out_cosio  = cos(inclo)
	out_cosio2 = out_cosio * out_cosio

	/* ------------------ un-kozai the mean motion ----------------- */
	ak    = pow(gravityConstants.xke / out_no, x2o3)
	d1    = 0.75 * gravityConstants.j2 * (3.0 * out_cosio2 - 1.0) / (out_rteosq * out_omeosq)
	del   = d1 / (ak * ak)
	adel  = ak * (1.0 - del * del - del *
	(1.0 / 3.0 + 134.0 * del * del / 81.0))
	del   = d1/(adel * adel)
	out_no    = out_no / (1.0 + del)

	out_ao    = pow(gravityConstants.xke / out_no, x2o3)
	out_sinio = sin(inclo)
	po    = out_ao * out_omeosq
	out_con42 = 1.0 - 5.0 * out_cosio2
	out_con41 = -out_con42-out_cosio2-out_cosio2
	out_ainv  = 1.0 / out_ao
	out_posq  = po * po
	out_rp    = out_ao * (1.0 - ecco)
	out_method = "n"

	// sgp4fix modern approach to finding sidereal time
	if opsmode == "a" {
		// sgp4fix use old way of finding gst
		// count integer number of days from 0 jan 1970
		ts70   = epoch - 7305.0
		ds70   = floor(ts70 + 1.0e-8)
		tfrac  = ts70 - ds70
		// find greenwich location at epoch
		c1     = 1.72027916940703639e-2
		thgr70 = 1.7321343856509374
		fk5r   = 5.07551419432269442e-15
		c1p2p  = c1 + Double.twoπ
		out_gsto  = fmod( thgr70 + c1*ds70 + c1p2p*tfrac + ts70*ts70*fk5r, Double.twoπ)
		if out_gsto < 0.0 {
			out_gsto = out_gsto + Double.twoπ
		}
	} else {
		out_gsto = gstime(epoch + 2433281.5)
	}


	//#include "debug5.cpp"
}  // end initl

/*-----------------------------------------------------------------------------
*
*                             procedure sgp4init
*
*  this procedure initializes variables for sgp4.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    opsmode     - mode of operation afspc or improved 'a', 'i'
*    whichconst  - which set of constants to use  72, 84
*    satn        - satellite number
*    bstar       - sgp4 type drag coefficient              kg/m2er
*    ecco        - eccentricity
*    epoch       - epoch time in days from jan 0, 1950. 0 hr
*    argpo       - argument of perigee (output if ds)
*    inclo       - inclination
*    mo          - mean anomaly (output if ds)
*    no          - mean motion
*    nodeo       - right ascension of ascending node
*
*  outputs       :
*    satrec      - common values for subsequent calls
*    return code - non-zero on error.
*                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
*                   2 - mean motion less than 0.0
*                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
*                   4 - semi-latus rectum < 0.0
*                   5 - epoch elements are sub-orbital
*                   6 - satellite has decayed
*
*  locals        :
*    cnodm  , snodm  , cosim  , sinim  , cosomm , sinomm
*    cc1sq  , cc2    , cc3
*    coef   , coef1
*    cosio4      -
*    day         -
*    dndt        -
*    em          - eccentricity
*    emsq        - eccentricity squared
*    eeta        -
*    etasq       -
*    gam         -
*    argpm       - argument of perigee
*    nodem       -
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    perige      - perigee
*    pinvsq      -
*    psisq       -
*    qzms24      -
*    rtemsq      -
*    s1, s2, s3, s4, s5, s6, s7          -
*    sfour       -
*    ss1, ss2, ss3, ss4, ss5, ss6, ss7         -
*    sz1, sz2, sz3
*    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
*    tc          -
*    temp        -
*    temp1, temp2, temp3       -
*    tsi         -
*    xpidot      -
*    xhdot1      -
*    z1, z2, z3          -
*    z11, z12, z13, z21, z22, z23, z31, z32, z33         -
*
*  coupling      :
*    initl       -
*    dscom       -
*    dpper       -
*    dsinit      -
*    sgp4        -
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/

func sgp4init
(
	_ whichconst : GravityConstantsType, _ opsmode : Character, _ satn : Int,  _ epoch : Double,
	_ xbstar : Double,  _ xecco : Double, _ xargpo : Double,
	_ xinclo : Double,  _ xmo : Double,   _ xno : Double,
	_ xnodeo : Double, _ satrec : Elsetrec
	) ->  Elsetrec {
	var out_satrec = satrec
	/* --------------------- local variables ------------------------ */
	var ao = 0.0, ainv = 0.0,   con42 = 0.0, cosio = 0.0, sinio = 0.0, cosio2 = 0.0, eccsq = 0.0,
	omeosq = 0.0, posq = 0.0,   rp = 0.0,     rteosq = 0.0,
	cnodm  = 0.0, snodm  = 0.0, cosim  = 0.0, sinim  = 0.0, cosomm = 0.0, sinomm = 0.0, cc1sq  = 0.0,
	cc2    = 0.0, cc3    = 0.0, coef   = 0.0, coef1  = 0.0, cosio4 = 0.0, day    = 0.0, dndt   = 0.0,
	em     = 0.0, emsq   = 0.0, eeta   = 0.0, etasq  = 0.0, gam    = 0.0, argpm  = 0.0, nodem  = 0.0,
	inclm  = 0.0, mm     = 0.0, nm     = 0.0, perige = 0.0, pinvsq = 0.0, psisq  = 0.0, qzms24 = 0.0,
	rtemsq = 0.0, s1     = 0.0, s2     = 0.0, s3     = 0.0, s4     = 0.0, s5     = 0.0, s6     = 0.0,
	s7     = 0.0, sfour  = 0.0, ss1    = 0.0, ss2    = 0.0, ss3    = 0.0, ss4    = 0.0, ss5    = 0.0,
	ss6    = 0.0, ss7    = 0.0, sz1    = 0.0, sz2    = 0.0, sz3    = 0.0, sz11   = 0.0, sz12   = 0.0,
	sz13   = 0.0, sz21   = 0.0, sz22   = 0.0, sz23   = 0.0, sz31   = 0.0, sz32   = 0.0, sz33   = 0.0,
	tc     = 0.0, temp   = 0.0, temp1  = 0.0, temp2  = 0.0, temp3  = 0.0, tsi    = 0.0, xpidot = 0.0,
	xhdot1 = 0.0, z1     = 0.0, z2     = 0.0, z3     = 0.0, z11    = 0.0, z12    = 0.0, z13    = 0.0,
	z21    = 0.0, z22    = 0.0, z23    = 0.0, z31    = 0.0, z32    = 0.0, z33 = 0.0,
	qzms2t = 0.0, ss = 0.0, x2o3 = 0.0

	/* ------------------------ initialization --------------------- */
	// sgp4fix divisor for divide by zero check on inclination
	// the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
	// 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
	let temp4    =   1.5e-12

	/* ----------- set all near earth variables to zero ------------ */
	out_satrec.isimp   = 0;   out_satrec.method = "n"; out_satrec.aycof    = 0.0;
	out_satrec.con41   = 0.0; out_satrec.cc1    = 0.0; out_satrec.cc4      = 0.0;
	out_satrec.cc5     = 0.0; out_satrec.d2     = 0.0; out_satrec.d3       = 0.0;
	out_satrec.d4      = 0.0; out_satrec.delmo  = 0.0; out_satrec.eta      = 0.0;
	out_satrec.argpdot = 0.0; out_satrec.omgcof = 0.0; out_satrec.sinmao   = 0.0;
	out_satrec.t       = 0.0; out_satrec.t2cof  = 0.0; out_satrec.t3cof    = 0.0;
	out_satrec.t4cof   = 0.0; out_satrec.t5cof  = 0.0; out_satrec.x1mth2   = 0.0;
	out_satrec.x7thm1  = 0.0; out_satrec.mdot   = 0.0; out_satrec.nodedot  = 0.0;
	out_satrec.xlcof   = 0.0; out_satrec.xmcof  = 0.0; out_satrec.nodecf   = 0.0;

	/* ----------- set all deep space variables to zero ------------ */
	out_satrec.irez  = 0;   out_satrec.d2201 = 0.0; out_satrec.d2211 = 0.0;
	out_satrec.d3210 = 0.0; out_satrec.d3222 = 0.0; out_satrec.d4410 = 0.0;
	out_satrec.d4422 = 0.0; out_satrec.d5220 = 0.0; out_satrec.d5232 = 0.0;
	out_satrec.d5421 = 0.0; out_satrec.d5433 = 0.0; out_satrec.dedt  = 0.0;
	out_satrec.del1  = 0.0; out_satrec.del2  = 0.0; out_satrec.del3  = 0.0;
	out_satrec.didt  = 0.0; out_satrec.dmdt  = 0.0; out_satrec.dnodt = 0.0;
	out_satrec.domdt = 0.0; out_satrec.e3    = 0.0; out_satrec.ee2   = 0.0;
	out_satrec.peo   = 0.0; out_satrec.pgho  = 0.0; out_satrec.pho   = 0.0;
	out_satrec.pinco = 0.0; out_satrec.plo   = 0.0; out_satrec.se2   = 0.0;
	out_satrec.se3   = 0.0; out_satrec.sgh2  = 0.0; out_satrec.sgh3  = 0.0;
	out_satrec.sgh4  = 0.0; out_satrec.sh2   = 0.0; out_satrec.sh3   = 0.0;
	out_satrec.si2   = 0.0; out_satrec.si3   = 0.0; out_satrec.sl2   = 0.0;
	out_satrec.sl3   = 0.0; out_satrec.sl4   = 0.0; out_satrec.gsto  = 0.0;
	out_satrec.xfact = 0.0; out_satrec.xgh2  = 0.0; out_satrec.xgh3  = 0.0;
	out_satrec.xgh4  = 0.0; out_satrec.xh2   = 0.0; out_satrec.xh3   = 0.0;
	out_satrec.xi2   = 0.0; out_satrec.xi3   = 0.0; out_satrec.xl2   = 0.0;
	out_satrec.xl3   = 0.0; out_satrec.xl4   = 0.0; out_satrec.xlamo = 0.0;
	out_satrec.zmol  = 0.0; out_satrec.zmos  = 0.0; out_satrec.atime = 0.0;
	out_satrec.xli   = 0.0; out_satrec.xni   = 0.0;

	// sgp4fix - note the following variables are also passed directly via satrec.
	// it is possible to streamline the sgp4init call by deleting the "x"
	// variables, but the user would need to set the satrec.* values first. we
	// include the additional assignments in case twoline2rv is not used.
	out_satrec.bstar   = xbstar
	out_satrec.ecco    = xecco
	out_satrec.argpo   = xargpo
	out_satrec.inclo   = xinclo
	out_satrec.mo	    = xmo
	out_satrec.no	    = xno
	out_satrec.nodeo   = xnodeo

	// sgp4fix add opsmode
	out_satrec.operationmode = opsmode

	/* ------------------------ earth constants ----------------------- */
	// sgp4fix identify constants and allow alternate values
	let gravityConstants = GravityConstants(type: whichconst)
	ss     = 78.0 / gravityConstants.radiusearthkm + 1.0
	qzms2t = pow(((120.0 - 78.0) / gravityConstants.radiusearthkm), 4)
	x2o3   =  2.0 / 3.0

	out_satrec.inited = false
	out_satrec.t	 = 0.0

	initl(satn, whichconst, out_satrec.ecco, epoch, out_satrec.inclo, &(out_satrec.no), &(out_satrec.method),
		  &ainv, &ao, &(out_satrec.con41), &con42, &cosio, &cosio2, &eccsq, &omeosq,
		  &posq, &rp, &rteosq, &sinio, &(out_satrec.gsto), out_satrec.operationmode
	)
	out_satrec.error = .success

	if rp < 1.0 {
		//         printf("# *** satn%d epoch elts sub-orbital ***\n", satn)
		out_satrec.error = .epochElementAreSubOrbital(rp)
	}

	if (omeosq >= 0.0 ) || ( out_satrec.no >= 0.0) {
		out_satrec.isimp = 0
		if rp < (220.0 / gravityConstants.radiusearthkm + 1.0) {
			out_satrec.isimp = 1
		}
		sfour  = ss
		qzms24 = qzms2t
		perige = (rp - 1.0) * gravityConstants.radiusearthkm

		/* - for perigees below 156 km, s and qoms2t are altered - */
		if perige < 156.0 {
			sfour = perige - 78.0
			if perige < 98.0 {
				sfour = 20.0
			}
			qzms24 = pow(((120.0 - sfour) / gravityConstants.radiusearthkm), 4.0)
			sfour  = sfour / gravityConstants.radiusearthkm + 1.0
		}
		pinvsq = 1.0 / posq

		tsi  = 1.0 / (ao - sfour)
		out_satrec.eta  = ao * out_satrec.ecco * tsi
		etasq = out_satrec.eta * out_satrec.eta
		eeta  = out_satrec.ecco * out_satrec.eta
		psisq = fabs(1.0 - etasq)
		coef  = qzms24 * pow(tsi, 4.0)
		coef1 = coef / pow(psisq, 3.5)
		cc2   = coef1 * out_satrec.no * (ao * (1.0 + 1.5 * etasq + eeta *
		(4.0 + etasq)) + 0.375 * gravityConstants.j2 * tsi / psisq * out_satrec.con41 *
		(8.0 + 3.0 * etasq * (8.0 + etasq)))
		out_satrec.cc1   = out_satrec.bstar * cc2
		cc3   = 0.0
		if out_satrec.ecco > 1.0e-4 {
			cc3 = -2.0 * coef * tsi * gravityConstants.j3oj2 * out_satrec.no * sinio / out_satrec.ecco
		}
		out_satrec.x1mth2 = 1.0 - cosio2
		out_satrec.cc4    = 2.0 * out_satrec.no * coef1 * ao * omeosq *
		(out_satrec.eta * (2.0 + 0.5 * etasq) + out_satrec.ecco *
		(0.5 + 2.0 * etasq) - gravityConstants.j2 * tsi / (ao * psisq) *
		(-3.0 * out_satrec.con41 * (1.0 - 2.0 * eeta + etasq *
		(1.5 - 0.5 * eeta)) + 0.75 * out_satrec.x1mth2 *
		(2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * out_satrec.argpo)))
		out_satrec.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *
		(etasq + eeta) + eeta * etasq)
		cosio4 = cosio2 * cosio2
		temp1  = 1.5 * gravityConstants.j2 * pinvsq * out_satrec.no
		temp2  = 0.5 * temp1 * gravityConstants.j2 * pinvsq
		temp3  = -0.46875 * gravityConstants.j4 * pinvsq * pinvsq * out_satrec.no
		out_satrec.mdot     = out_satrec.no + 0.5 * temp1 * rteosq * out_satrec.con41 + 0.0625 *
		temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4)
		out_satrec.argpdot  = -0.5 * temp1 * con42 + 0.0625 * temp2 *
		(7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
		temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4)
		xhdot1            = -temp1 * cosio
		out_satrec.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
		2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio
		xpidot            =  out_satrec.argpdot + out_satrec.nodedot
		out_satrec.omgcof   = out_satrec.bstar * cc3 * cos(out_satrec.argpo)
		out_satrec.xmcof    = 0.0
		if out_satrec.ecco > 1.0e-4 {
			out_satrec.xmcof = -x2o3 * coef * out_satrec.bstar / eeta
		}
		out_satrec.nodecf = 3.5 * omeosq * xhdot1 * out_satrec.cc1
		out_satrec.t2cof   = 1.5 * out_satrec.cc1
		// sgp4fix for divide by zero with xinco = 180 deg
		if fabs(cosio+1.0) > 1.5e-12 {
			out_satrec.xlcof = -0.25 * gravityConstants.j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio)
		} else {
			out_satrec.xlcof = -0.25 * gravityConstants.j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4
		}
		out_satrec.aycof   = -0.5 * gravityConstants.j3oj2 * sinio
		out_satrec.delmo   = pow((1.0 + out_satrec.eta * cos(out_satrec.mo)), 3)
		out_satrec.sinmao  = sin(out_satrec.mo)
		out_satrec.x7thm1  = 7.0 * cosio2 - 1.0

		/* --------------- deep space initialization ------------- */
		if (Double.twoπ / out_satrec.no) >= 225.0 {
			out_satrec.method = "d"
			out_satrec.isimp  = 1
			tc    =  0.0
			inclm = out_satrec.inclo

			dscom(epoch, out_satrec.ecco, out_satrec.argpo, tc, out_satrec.inclo, out_satrec.nodeo,
				  out_satrec.no, &snodm, &cnodm,  &sinim, &cosim,&sinomm,     &cosomm,
				  &day, &(out_satrec.e3), &(out_satrec.ee2), &em,         &emsq, &gam,
				  &(out_satrec.peo),  &(out_satrec.pgho),   &(out_satrec.pho), &(out_satrec.pinco),
				  &(out_satrec.plo),  &rtemsq,        &(out_satrec.se2), &(out_satrec.se3),
				  &(out_satrec.sgh2), &(out_satrec.sgh3),   &(out_satrec.sgh4),
				  &(out_satrec.sh2),  &(out_satrec.sh3),    &(out_satrec.si2), &(out_satrec.si3),
				  &(out_satrec.sl2),  &(out_satrec.sl3),    &(out_satrec.sl4), &s1, &s2, &s3, &s4, &s5,
				  &s6,   &s7,   &ss1,  &ss2,  &ss3,  &ss4,  &ss5,  &ss6,  &ss7, &sz1, &sz2, &sz3,
				  &sz11, &sz12, &sz13, &sz21, &sz22, &sz23, &sz31, &sz32, &sz33,
				  &out_satrec.xgh2, &out_satrec.xgh3,   &out_satrec.xgh4, &out_satrec.xh2,
				  &out_satrec.xh3,  &out_satrec.xi2,    &out_satrec.xi3,  &out_satrec.xl2,
				  &out_satrec.xl3,  &out_satrec.xl4,    &nm, &z1, &z2, &z3, &z11,
				  &z12, &z13, &z21, &z22, &z23, &z31, &z32, &z33,
				  &out_satrec.zmol, &out_satrec.zmos
			)
			dpper(out_satrec.e3, out_satrec.ee2, out_satrec.peo, out_satrec.pgho,
				  out_satrec.pho, out_satrec.pinco, out_satrec.plo, out_satrec.se2,
				  out_satrec.se3, out_satrec.sgh2, out_satrec.sgh3, out_satrec.sgh4,
				  out_satrec.sh2, out_satrec.sh3, out_satrec.si2,  out_satrec.si3,
				  out_satrec.sl2, out_satrec.sl3, out_satrec.sl4,  out_satrec.t,
				  out_satrec.xgh2,out_satrec.xgh3,out_satrec.xgh4, out_satrec.xh2,
				  out_satrec.xh3, out_satrec.xi2, out_satrec.xi3,  out_satrec.xl2,
				  out_satrec.xl3, out_satrec.xl4, out_satrec.zmol, out_satrec.zmos, inclm, out_satrec.inited,
				  &out_satrec.ecco, &out_satrec.inclo, &out_satrec.nodeo, &out_satrec.argpo, &out_satrec.mo,
				  out_satrec.operationmode
			)

			argpm  = 0.0
			nodem  = 0.0
			mm     = 0.0

			dsinit(whichconst,
				   cosim, &emsq, out_satrec.argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4,
				   ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, out_satrec.t, tc,
				   out_satrec.gsto, out_satrec.mo, out_satrec.mdot, out_satrec.no, out_satrec.nodeo,
				   out_satrec.nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33,
				   out_satrec.ecco, eccsq, &em, &argpm, &inclm, &mm, &nm, &nodem,
				   &out_satrec.irez,  &out_satrec.atime,
				   &out_satrec.d2201, &out_satrec.d2211, &out_satrec.d3210, &out_satrec.d3222 ,
				   &out_satrec.d4410, &out_satrec.d4422, &out_satrec.d5220, &out_satrec.d5232,
				   &out_satrec.d5421, &out_satrec.d5433, &out_satrec.dedt,  &out_satrec.didt,
				   &out_satrec.dmdt,  &dndt,         &out_satrec.dnodt, &out_satrec.domdt ,
				   &out_satrec.del1,  &out_satrec.del2,  &out_satrec.del3,  &out_satrec.xfact,
				   &out_satrec.xlamo, &out_satrec.xli,   &out_satrec.xni
			)
		}

		/* ----------- set variables if not deep space ----------- */
		if (out_satrec.isimp != 1)
		{
			cc1sq          = out_satrec.cc1 * out_satrec.cc1
			out_satrec.d2    = 4.0 * ao * tsi * cc1sq
			temp           = out_satrec.d2 * tsi * out_satrec.cc1 / 3.0
			out_satrec.d3    = (17.0 * ao + sfour) * temp
			out_satrec.d4    = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) *
				out_satrec.cc1
			out_satrec.t3cof = out_satrec.d2 + 2.0 * cc1sq
			out_satrec.t4cof = 0.25 * (3.0 * out_satrec.d3 + out_satrec.cc1 *
				(12.0 * out_satrec.d2 + 10.0 * cc1sq))
			out_satrec.t5cof = 0.2 * (3.0 * out_satrec.d4 +
				12.0 * out_satrec.cc1 * out_satrec.d3 +
				6.0 * out_satrec.d2 * out_satrec.d2 +
				15.0 * cc1sq * (2.0 * out_satrec.d2 + cc1sq))
		}
	} // if omeosq = 0 ...

	/* finally propogate to zero epoch to initialise all others. */
	if out_satrec.error.isSuccess() {
		(out_satrec.error, _, _) = sgp4(whichconst, out_satrec, 0.0)
	}

	out_satrec.inited = true

	return out_satrec
}  // end sgp4init

/*-----------------------------------------------------------------------------
*
*                             procedure sgp4
*
*  this procedure is the sgp4 prediction model from space command. this is an
*    updated and combined version of sgp4 and sdp4, which were originally
*    published separately in spacetrack report #3. this version follows the
*    methodology from the aiaa paper (2006) describing the history and
*    development of the code.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    satrec	 - initialised structure from sgp4init() call.
*    tsince	 - time eince epoch (minutes)
*
*  outputs       :
*    r           - position vector                     km
*    v           - velocity                            km/sec
*  return code - non-zero on error.
*                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
*                   2 - mean motion less than 0.0
*                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
*                   4 - semi-latus rectum < 0.0
*                   5 - epoch elements are sub-orbital
*                   6 - satellite has decayed
*
*  locals        :
*    am          -
*    axnl, aynl        -
*    betal       -
*    cosim   , sinim   , cosomm  , sinomm  , cnod    , snod    , cos2u   ,
*    sin2u   , coseo1  , sineo1  , cosi    , sini    , cosip   , sinip   ,
*    cosisq  , cossu   , sinsu   , cosu    , sinu
*    delm        -
*    delomg      -
*    dndt        -
*    eccm        -
*    emsq        -
*    ecose       -
*    el2         -
*    eo1         -
*    eccp        -
*    esine       -
*    argpm       -
*    argpp       -
*    omgadf      -
*    pl          -
*    r           -
*    rtemsq      -
*    rdotl       -
*    rl          -
*    rvdot       -
*    rvdotl      -
*    su          -
*    t2  , t3   , t4    , tc
*    tem5, temp , temp1 , temp2  , tempa  , tempe  , templ
*    u   , ux   , uy    , uz     , vx     , vy     , vz
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    nodem       - right asc of ascending node
*    xinc        -
*    xincp       -
*    xl          -
*    xlm         -
*    mp          -
*    xmdf        -
*    xmx         -
*    xmy         -
*    nodedf      -
*    xnode       -
*    nodep       -
*    np          -
*
*  coupling      :
*    dpper
*    dpspace
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/

func sgp4 (
	_ whichconst : GravityConstantsType, _ satrec : Elsetrec,  _ tsince : Double) ->
	(status : SGP4Status, pos : Position, vel : Velocity) {
		var out_satrec = satrec
		var r = kVector3Null
		var v = kVector3Null
		var am = 0.0, axnl = 0.0,  aynl  = 0.0, betal  = 0.0,  cosim  = 0.0, cnod   = 0.0,
		cos2u = 0.0, cosi  = 0.0, cosip  = 0.0,  cosisq = 0.0, cossu  = 0.0, cosu = 0.0,
		delm  = 0.0, delomg = 0.0, em  = 0.0  /* , emsq (pas utilise)*/  ,  ecose  = 0.0, el2    = 0.0, eo1  = 0.0,
		ep    = 0.0, esine  = 0.0, argpm = 0.0, argpp  = 0.0,  argpdf = 0.0, pl = 0.0,
		mvt   = 0.0, rdotl  = 0.0, rl    = 0.0, rvdot  = 0.0,  rvdotl = 0.0, sinim  = 0.0,
		sin2u = 0.0, sini  = 0.0, sinip  = 0.0,  sinsu  = 0.0, sinu   = 0.0,
		snod  = 0.0, su     = 0.0, t2    = 0.0, t3     = 0.0,  t4     = 0.0, tem5   = 0.0, temp = 0.0,
		temp1 = 0.0, temp2  = 0.0, tempa = 0.0, tempe  = 0.0,  templ  = 0.0, u      = 0.0, ux   = 0.0,
		uy    = 0.0, uz     = 0.0, vx    = 0.0, vy     = 0.0,  vz     = 0.0, inclm  = 0.0, mm   = 0.0,
		nm    = 0.0, nodem = 0.0, xinc  = 0.0, xincp  = 0.0,  xl     = 0.0, xlm    = 0.0, mp   = 0.0,
		xmdf  = 0.0, xmx    = 0.0, xmy   = 0.0, nodedf = 0.0, xnode  = 0.0, nodep = 0.0, tc   = 0.0, dndt = 0.0,
		x2o3   = 0.0, vkmpersec  = 0.0
		var coseo1 = 0.0,  mrt = 0.0, sineo1 = 0.0
		var ktr : Int

		/* ------------------ set mathematical constants --------------- */
		// sgp4fix divisor for divide by zero check on inclination
		// the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
		// 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
		let temp4 =   1.5e-12
		x2o3  = 2.0 / 3.0
		// sgp4fix identify constants and allow alternate values
		let gravityConstants = GravityConstants(type: whichconst)
		vkmpersec     = gravityConstants.radiusearthkm * gravityConstants.xke/60.0

		/* --------------------- clear sgp4 error flag ----------------- */
		out_satrec.t     = tsince
		out_satrec.error = .success

		/* ------- update for secular gravity and atmospheric drag ----- */
		xmdf    = out_satrec.mo + out_satrec.mdot * out_satrec.t
		argpdf  = out_satrec.argpo + out_satrec.argpdot * out_satrec.t
		nodedf  = out_satrec.nodeo + out_satrec.nodedot * out_satrec.t
		argpm   = argpdf
		mm      = xmdf
		t2      = out_satrec.t * out_satrec.t
		nodem   = nodedf + out_satrec.nodecf * t2
		tempa   = 1.0 - out_satrec.cc1 * out_satrec.t
		tempe   = out_satrec.bstar * out_satrec.cc4 * out_satrec.t
		templ   = out_satrec.t2cof * t2

		if out_satrec.isimp != 1 {
			delomg = out_satrec.omgcof * out_satrec.t
			delm   = out_satrec.xmcof *
				(pow((1.0 + out_satrec.eta * cos(xmdf)), 3) -
					out_satrec.delmo)
			temp   = delomg + delm
			mm     = xmdf + temp
			argpm  = argpdf - temp
			t3     = t2 * out_satrec.t
			t4     = t3 * out_satrec.t
			tempa  = tempa - out_satrec.d2 * t2 - out_satrec.d3 * t3 -
				out_satrec.d4 * t4
			tempe  = tempe + out_satrec.bstar * out_satrec.cc5 * (sin(mm) -
				out_satrec.sinmao)
			templ  = templ + out_satrec.t3cof * t3 + t4 * (out_satrec.t4cof +
				out_satrec.t * out_satrec.t5cof)
		}

		nm    = out_satrec.no
		em    = out_satrec.ecco
		inclm = out_satrec.inclo
		if out_satrec.method == "d" {
			tc = out_satrec.t
			dspace(out_satrec.irez,
				   out_satrec.d2201, out_satrec.d2211, out_satrec.d3210,
				   out_satrec.d3222, out_satrec.d4410, out_satrec.d4422,
				   out_satrec.d5220, out_satrec.d5232, out_satrec.d5421,
				   out_satrec.d5433, out_satrec.dedt,  out_satrec.del1,
				   out_satrec.del2,  out_satrec.del3,  out_satrec.didt,
				   out_satrec.dmdt,  out_satrec.dnodt, out_satrec.domdt,
				   out_satrec.argpo, out_satrec.argpdot, out_satrec.t, tc,
				   out_satrec.gsto, out_satrec.xfact, out_satrec.xlamo,
				   out_satrec.no, &out_satrec.atime,
				   &em, &argpm, &inclm, &out_satrec.xli, &mm, &out_satrec.xni,
				   &nodem, &dndt, &nm
			)
		} // if method = d

		if nm <= 0.0 {
			//         printf("# error nm %f\n", nm)
			out_satrec.error = .meanMotion(nm)
		}
		am = pow((gravityConstants.xke / nm),x2o3) * tempa * tempa
		nm = gravityConstants.xke / pow(am, 1.5)
		em = em - tempe

		// fix tolerance for error recognition
		if (em >= 1.0) || (em < -0.001) || (am < 0.95) {
			//         printf("# error em %f\n", em)
			out_satrec.error = .meanElement(em,am)
		}

		//   sgp4fix change test condition for eccentricity
		if em < 1.0e-6 {
			em  = 1.0e-6
		}
		mm     = mm + out_satrec.no * templ
		xlm    = mm + argpm + nodem
		// emsq   = em * em // analyse : contenu pas lu
		// temp   = 1.0 - emsq // analyse : contenu pas lu

		nodem  = fmod(nodem, Double.twoπ)
		argpm  = fmod(argpm, Double.twoπ)
		xlm    = fmod(xlm, Double.twoπ)
		mm     = fmod(xlm - argpm - nodem, Double.twoπ)

		/* ----------------- compute extra mean quantities ------------- */
		sinim = sin(inclm)
		cosim = cos(inclm)

		/* -------------------- add lunar-solar periodics -------------- */
		ep     = em
		xincp  = inclm
		argpp  = argpm
		nodep  = nodem
		mp     = mm
		sinip  = sinim
		cosip  = cosim
		if out_satrec.method == "d" {
			dpper(out_satrec.e3,   out_satrec.ee2,  out_satrec.peo,
				  out_satrec.pgho, out_satrec.pho,  out_satrec.pinco,
				  out_satrec.plo,  out_satrec.se2,  out_satrec.se3,
				  out_satrec.sgh2, out_satrec.sgh3, out_satrec.sgh4,
				  out_satrec.sh2,  out_satrec.sh3,  out_satrec.si2,
				  out_satrec.si3,  out_satrec.sl2,  out_satrec.sl3,
				  out_satrec.sl4,  out_satrec.t,    out_satrec.xgh2,
				  out_satrec.xgh3, out_satrec.xgh4, out_satrec.xh2,
				  out_satrec.xh3,  out_satrec.xi2,  out_satrec.xi3,
				  out_satrec.xl2,  out_satrec.xl3,  out_satrec.xl4,
				  out_satrec.zmol, out_satrec.zmos, out_satrec.inclo,
				  false, &ep, &xincp, &nodep, &argpp, &mp, out_satrec.operationmode
			)
			if xincp < 0.0 {
				xincp  = -xincp
				nodep = nodep + π
				argpp  = argpp - π
			}
			if (ep < 0.0) || (ep > 1.0) {
				//            printf("# error ep %f\n", ep)
				out_satrec.error = .perturbationElements(ep)
			}
		} // if method = d

		/* -------------------- long period periodics ------------------ */
		if out_satrec.method == "d" {
			sinip =  sin(xincp)
			cosip =  cos(xincp)
			out_satrec.aycof = -0.5 * gravityConstants.j3oj2 * sinip
			// sgp4fix for divide by zero for xincp = 180 deg
			if fabs(cosip+1.0) > 1.5e-12 {
				out_satrec.xlcof = -0.25 * gravityConstants.j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip)
			} else {
				out_satrec.xlcof = -0.25 * gravityConstants.j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4
			}
		}
		axnl = ep * cos(argpp)
		temp = 1.0 / (am * (1.0 - ep * ep))
		aynl = ep * sin(argpp) + temp * out_satrec.aycof
		xl   = mp + argpp + nodep + temp * out_satrec.xlcof * axnl

		/* --------------------- solve kepler's equation --------------- */
		u    = fmod(xl - nodep, Double.twoπ)
		eo1  = u
		tem5 = 9999.9
		ktr = 1
		//   sgp4fix for kepler iteration
		//   the following iteration needs better limits on corrections
		while ( fabs(tem5) >= 1.0e-12) && (ktr <= 10) {
			sineo1 = sin(eo1)
			coseo1 = cos(eo1)
			tem5   = 1.0 - coseo1 * axnl - sineo1 * aynl
			tem5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5
			if fabs(tem5) >= 0.95 {
				tem5 = tem5 > 0.0 ? 0.95 : -0.95
			}
			eo1    = eo1 + tem5
			ktr = ktr + 1
		}

		/* ------------- short period preliminary quantities ----------- */
		ecose = axnl*coseo1 + aynl*sineo1
		esine = axnl*sineo1 - aynl*coseo1
		el2   = axnl*axnl + aynl*aynl
		pl    = am*(1.0-el2)
		if pl < 0.0 {
			//         printf("# error pl %f\n", pl)
			out_satrec.error = .semiLatudRectum(pl)
		} else {
			rl     = am * (1.0 - ecose)
			rdotl  = √(am) * esine/rl
			rvdotl = √(pl) / rl
			betal  = √(1.0 - el2)
			temp   = esine / (1.0 + betal)
			sinu   = am / rl * (sineo1 - aynl - axnl * temp)
			cosu   = am / rl * (coseo1 - axnl + aynl * temp)
			su     = ArcTan(sinu, cosu)
			sin2u  = (cosu + cosu) * sinu
			cos2u  = 1.0 - 2.0 * sinu * sinu
			temp   = 1.0 / pl
			temp1  = 0.5 * gravityConstants.j2 * temp
			temp2  = temp1 * temp

			/* -------------- update for short period periodics ------------ */
			if out_satrec.method == "d" {
				cosisq                 = cosip * cosip
				out_satrec.con41  = 3.0*cosisq - 1.0
				out_satrec.x1mth2 = 1.0 - cosisq
				out_satrec.x7thm1 = 7.0*cosisq - 1.0
			}
			mrt   = rl * (1.0 - 1.5 * temp2 * betal * out_satrec.con41) +
				0.5 * temp1 * out_satrec.x1mth2 * cos2u
			su    = su - 0.25 * temp2 * out_satrec.x7thm1 * sin2u
			xnode = nodep + 1.5 * temp2 * cosip * sin2u
			xinc  = xincp + 1.5 * temp2 * cosip * sinip * cos2u
			mvt   = rdotl - nm * temp1 * out_satrec.x1mth2 * sin2u / gravityConstants.xke
			rvdot = rvdotl + nm * temp1 * (out_satrec.x1mth2 * cos2u +
				1.5 * out_satrec.con41) / gravityConstants.xke

			/* --------------------- orientation vectors ------------------- */
			sinsu =  sin(su)
			cossu =  cos(su)
			snod  =  sin(xnode)
			cnod  =  cos(xnode)
			sini  =  sin(xinc)
			cosi  =  cos(xinc)
			xmx   = -snod * cosi
			xmy   =  cnod * cosi
			ux    =  xmx * sinsu + cnod * cossu
			uy    =  xmy * sinsu + snod * cossu
			uz    =  sini * sinsu
			vx    =  xmx * cossu - cnod * sinsu
			vy    =  xmy * cossu - snod * sinsu
			vz    =  sini * cossu

			/* --------- position and velocity (in km and km/sec) ---------- */
			r[0] = (mrt * ux) * gravityConstants.radiusearthkm
			r[1] = (mrt * uy) * gravityConstants.radiusearthkm
			r[2] = (mrt * uz) * gravityConstants.radiusearthkm
			v[0] = (mvt * ux + rvdot * vx) * vkmpersec
			v[1] = (mvt * uy + rvdot * vy) * vkmpersec
			v[2] = (mvt * uz + rvdot * vz) * vkmpersec

			/*
			printf("Pos %lf %lf %lf\n",r[0],r[1],r[2])
			printf("Vit %lf %lf %lf\n",v[0],v[1],v[2])
			*/
		}  // if pl > 0

		// sgp4fix for decaying satellites
		if (mrt < 1.0)
		{
			//         printf("# decay condition %11.6f \n",mrt)
			out_satrec.error = .decayed
		}


		//#include "debug7.cpp"
		return (out_satrec.error, r, v)
}  // end sgp4


/* -----------------------------------------------------------------------------
*
*                           function gstime
*
*  this function finds the greenwich sidereal time.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    jdut1       - julian date in ut1             days from 4713 bc
*
*  outputs       :
*    gstime      - greenwich sidereal time        0 to 2pi rad
*
*  locals        :
*    temp        - temporary variable for doubles   rad
*    tut1        - julian centuries from the
*                  jan 1, 2000 12 h epoch (ut1)
*
*  coupling      :
*    none
*
*  references    :
*    vallado       2004, 191, eq 3-45
* --------------------------------------------------------------------------- */

func  gstime(_ jdut1 : Double) -> Double {
	var temp, tut1 : Double

	tut1 = (jdut1 - 2451545.0) / 36525.0
	temp = -6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 +
	(876600.0 * 3600 + 8640184.812866) * tut1 + 67310.54841  // sec
    temp = fmod(temp.degToRad / 240.0, Double.twoπ) //360/86400 = 1/240, to deg, to rad

	// ------------------------ check quadrants ---------------------
	if (temp < 0.0) {
		temp += Double.twoπ
	}

	return temp
}  // end gstime





