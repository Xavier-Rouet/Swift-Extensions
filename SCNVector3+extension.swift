//
//  SCNVector3.swift
//
//  Created by Xav perso on 10/12/2018.
//  Copyright © 2018 P'tit Xav. All rights reserved.
//

import Foundation
import SceneKit

extension SCNVector3 {
	func magnitude() -> Float {
		return √((x * x) + (y * y) + (z * z))
	}

	func normalize() -> SCNVector3 {
		let vecMag : Float = self.magnitude()
		if ( vecMag == 0.0 )
		{
			return SCNVector3(1.0, 0.0, 0.0)
		}
		return SCNVector3(x/vecMag, y/vecMag, z/vecMag)
	}

	static func normalizedVector3(_ vector : Position) -> SCNVector3 {
		return SCNVector3Make(v:vector).normalize()
	}

	// Vector dot product
	static func ⋅(_ vector1 : SCNVector3, _ vector2 : SCNVector3) -> Float{
		return vector1.x*vector2.x + vector1.y*vector2.y + vector1.z*vector2.z
	}

//	static func dotProduct(_ vector1 : SCNVector3, _ vector2 : SCNVector3) -> Float{
//		return vector1 * vector2
////		return vector1.x*vector2.x + vector1.y*vector2.y + vector1.z*vector2.z
//	}

	// Vector cross product
	static func ∧(vector1 : SCNVector3, vector2 : SCNVector3) -> SCNVector3 {
		return SCNVector3(x: (vector1.y * vector2.z) - (vector1.z * vector2.y),
						  y: (vector1.z * vector2.x) - (vector1.x * vector2.z),
						  z: (vector1.x * vector2.y) - (vector1.y * vector2.x))
	}

//	static func crossProduct(_ vector1 : SCNVector3, _ vector2 : SCNVector3) -> SCNVector3 {
//		return vector1 ^ vector2
////		return SCNVector3(x: (vector1.y * vector2.z) - (vector1.z * vector2.y),
////						  y: (vector1.z * vector2.x) - (vector1.x * vector2.z),
////						  z: (vector1.x * vector2.y) - (vector1.y * vector2.x))
//	}

	// Sizing of a vector
	static func *(scale : Float, vector: SCNVector3) -> SCNVector3 {
		let vn = vector.normalize()
		return SCNVector3(x: vn.x * scale,
						  y: vn.y * scale,
						  z: vn.z * scale)
	}

}
