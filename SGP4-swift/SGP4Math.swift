//
// SGP4Math.swift
//  pxSat3D-SK
//
//  Created by Xav perso on 23/08/2020.
//  Copyright © 2020 P'tit Xav. All rights reserved.
//

import Foundation

let π = Double.pi
//let	halfπ	= 0.5 * π
//let	twoπ	= 2.0 * π
//let	deg2rad	= π / 180.0
//let	rad2deg	= 180.0 / π
//let rad2hour = 12.0 / π
let	xpdotp	= 1440.0 / (2.0 * π)	// 229.1831180523293

let small  = 0.00000001
let undefined = 999999.1
let infinite  = 999999.9

func ArcTan(_ s: Double,_ c : Double) -> Double {
	return atan2((s),(c))
}

func ArcTan(_ s : Double) -> Double {
	return atan(s)
}

func ArcSin(_ arg : Double) -> Double {
	if arg >= 1 {
        return Double.halfπ
	} else if arg <= -1 {
        return -Double.halfπ
	} else {
		return ArcTan(arg / √(1-Sqr(arg)))
	}
}

func ArcCos(_ arg : Double) -> Double
{
    return Double.halfπ - ArcSin(arg)
}

// Vectors
struct SGP4Vector3 {
	var x, y, z : Double

	init(x : Double, y : Double, z : Double) {
		self.x = x
		self.y = y
		self.z = z
	}

	init(_ d : [Double]) {
		x = d[0]
		y = d[1]
		z = d[2]
	}
	subscript(_ index : Int) -> Double {
		set (double) {
			switch index {
			case 0:
				x = double
			case 1:
				y = double
			case 2:
				z = double
			default:
				print("Index out of bounds")
			}
		}

		get {
			switch index {
			case 0:
				return x
			case 1:
				return y
			case 2:
				return z
			default:
				print("Index out of bounds")
				return undefined
			}
		}
	}

	var mag:  Double {
		return √(x*x + y*y + z*z)
	}
	
	static func -(l : SGP4Vector3, r : SGP4Vector3) -> SGP4Vector3 {
		return SGP4Vector3(x: l.x - r.x,
						   y: l.y - r.y,
						   z: l.z - r.z)
	}
}

//typealias Position = [Double]
//typealias Velocity = [Double]
//
//let kVector3Null = [0.0, 0.0, 0.0]

// EVOLUTION :
typealias Position = SGP4Vector3
typealias Velocity = SGP4Vector3

let kVector3Null = SGP4Vector3(x: 0, y: 0, z: 0)

struct SGP4Vector4 {
	var x, y, z, m : Double

	subscript(_ index : Int) -> Double {
		set (double) {
			switch index {
			case 0:
				x = double
			case 1:
				y = double
			case 2:
				z = double
			case 3:
				m = double
			default:
				print("Index out of bounds")
			}
		}

		get {
			switch index {
			case 0:
				return x
			case 1:
				return y
			case 2:
				return z
			case 4:
				return m
			default:
				print("Index out of bounds")
				return undefined
			}
		}
	}

	mutating func mag() {
		m = √(x*x + y*y + z*z)
	}
}

typealias Observator = [Double]
let kVector4Null = [0.0, 0.0, 0.0, 0.0]

// OPTIMISER
func Magnitude(_ v : inout [Double]) {
	v[3] = √(Sqr(v[0]) + Sqr(v[1]) + Sqr(v[2]))
}
func Magnitude(_ v : inout SGP4Vector4) {
	v.mag()
}

/* -----------------------------------------------------------------------------
*
*                           function dot
*
*  this function finds the dot product of two vectors.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec1        - vector number 1
*    vec2        - vector number 2
*
*  outputs       :
*    dot         - result
*
*  locals        :
*    none.
*
*  coupling      :
*    none.
*
* --------------------------------------------------------------------------- */

//func  dot(_ x : [Double], _ y : [Double]) -> Double
//{
//	return (x[0]*y[0] + x[1]*y[1] + x[2]*y[2])
//}  // end dot

infix operator ⋅: MultiplicationPrecedence
func ⋅(_ v1 : [Double],_ v2: [Double]) -> Double {
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
}
func ⋅(_ v1 : SGP4Vector3,_ v2: SGP4Vector3) -> Double {
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
}

//func Modulus(_ value : Double, _ divider : Double) -> Double {
//	return fmod(value, divider)
//}

func %(value : Double, divider : Double) -> Double {
	return fmod(value, divider)
}

prefix operator √
prefix func √(_ a : Double) -> Double {
	return sqrt(a)
}

func Sqr(_ a : Double) -> Double{
	return ((a)*(a))
}

func Cos(_ a : Double) -> Double {
	return cos(a)
}

func Sin(_ a : Double) -> Double {
	return sin(a)
}

func Abs(_ a : Double) -> Double {
	return fabs(a)
}

func Tan(_ a : Double) -> Double {
	return tan(a)
}

func Radians(_ d : Double) -> Double {
	return ((d)*π/180.0)
}

func Degrees(_ r : Double) -> Double {
	return ((r)*180.0/π)
}

func sgn(_ x : Double) -> Double
{
	if (x < 0.0)
	{
		return -1.0
	}
	else
	{
		return 1.0
	}

}  // end sgn

/* -----------------------------------------------------------------------------
*
*                           function mag
*
*  this procedure finds the magnitude of a vector.  the tolerance is set to
*    0.000001, thus the 1.0e-12 for the squared test of underflows.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec         - vector
*
*  outputs       :
*    vec         - answer stored in fourth component
*
*  locals        :
*    none.
*
*  coupling      :
*    none.
* --------------------------------------------------------------------------- */

func mag(_ x : [Double]) -> Double {
	return √(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
}  // end mag

func mag(_ v : SGP4Vector3) -> Double {
	return v.mag
}  // end mag

/* -----------------------------------------------------------------------------
*
*                           procedure cross
*
*  this procedure crosses two vectors.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec1        - vector number 1
*    vec2        - vector number 2
*
*  outputs       :
*    outvec      - vector result of a x b
*
*  locals        :
*    none.
*
*  coupling      :
*    mag           magnitude of a vector
---------------------------------------------------------------------------- */

//func cross(_ vec1 : [Double], _ vec2 : [Double]) -> [Double]
//{
//	var outvec : [Double] = [0.0, 0.0, 0.0]
//	outvec[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1]
//	outvec[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2]
//	outvec[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0]
//	return outvec
//}  // end cross

infix operator ∧: MultiplicationPrecedence

func ∧(_ vec1 : [Double], _ vec2 : [Double]) -> [Double] {
	var outvec : [Double] = [0.0, 0.0, 0.0]
	outvec[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1]
	outvec[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2]
	outvec[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0]
	return outvec
}  // end cross

func ∧(_ vec1 : SGP4Vector3, _ vec2 : SGP4Vector3) -> SGP4Vector3 {
	return SGP4Vector3(x: vec1[1]*vec2[2] - vec1[2]*vec2[1],
				   y: vec1[2]*vec2[0] - vec1[0]*vec2[2],
				   z: vec1[0]*vec2[1] - vec1[1]*vec2[0])
}  // end cross

/* -----------------------------------------------------------------------------
*
*                           procedure angle
*
*  this procedure calculates the angle between two vectors.  the output is
*    set to 999999.1 to indicate an undefined value.  be sure to check for
*    this at the output phase.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec1        - vector number 1
*    vec2        - vector number 2
*
*  outputs       :
*    theta       - angle between the two vectors  -pi to pi
*
*  locals        :
*    temp        - temporary real variable
*
*  coupling      :
*    dot           dot product of two vectors
* --------------------------------------------------------------------------- */
func angle(_ vec1 : SGP4Vector3, _ vec2 : SGP4Vector3) -> Double {
	var magv1, magv2, temp : Double

	magv1 = mag(vec1)
	magv2 = mag(vec2)

	if (magv1*magv2 > small*small)
	{
		//		temp = dot(vec1,vec2) / (magv1*magv2)
		temp = vec1⋅vec2 / (magv1*magv2)
		if (abs( temp ) > 1.0) {
			temp = sgn(temp) * 1.0
		}
		return acos( temp )
	} else {
		return undefined
	}
}  // end angle

func angle(_ vec1 : [Double], _ vec2 : [Double]) -> Double {
	return angle(SGP4Vector3(vec1), SGP4Vector3(vec2))
}


/* -----------------------------------------------------------------------------
*
*                           function asinh
*
*  this function evaluates the inverse hyperbolic sine function.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    xval        - angle value                                  any real
*
*  outputs       :
*    arcsinh     - result                                       any real
*
*  locals        :
*    none.
*
*  coupling      :
*    none.
*
* --------------------------------------------------------------------------- */

func asinh(_ xval : Double) -> Double
{
	return log( xval + √( xval*xval + 1.0 ) )
}  // end asinh

