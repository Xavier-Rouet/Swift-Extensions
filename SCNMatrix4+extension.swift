//
//  SCNMatrix4+extension.swift
//
//  Created by Xav perso on 08/11/2020.
//  Copyright © 2020 P'tit Xav. All rights reserved.
//

import Foundation
import SceneKit

extension SCNMatrix4 {
	static func *(_ a : SCNMatrix4, _ b : SCNMatrix4) -> SCNMatrix4 {
		return SCNMatrix4Mult(a, b)
	}

	static func *=(_ a : inout SCNMatrix4, _ b : SCNMatrix4) {
		a = a * b
	}

	static func *(_ s : Float, _ m : SCNMatrix4) -> SCNMatrix4 {
		return SCNMatrix4Scale(m, s, s, s)
	}

	static func *(_ m : SCNMatrix4, _ s : Float) -> SCNMatrix4 {
		return SCNMatrix4Scale(m, s, s, s)
	}

	static func *=(_ m : inout SCNMatrix4, _ s : Float) {
		m = s * m
	}

    mutating func rotateX(_ inAngleX : Float) {
        let rotation = SCNMatrix4MakeRotation(inAngleX.degToRad, 1, 0, 0)
        self *= rotation
    }

    mutating func rotateY(_ inAngleY : Float) {
        let rotation = SCNMatrix4MakeRotation(inAngleY.degToRad, 0, 1, 0)
        self *= rotation
    }

    mutating func rotateZ(_ inAngleZ : Float) {
        let rotation = SCNMatrix4MakeRotation(inAngleZ.degToRad, 0, 0, 1)
        self *= rotation
    }

    mutating func rotate(x inX : Float, y inY : Float) {
        let lAngle = √(inX*inX+inY*inY)
        if lAngle != 0.0
        {
            let rotation = SCNMatrix4MakeRotation(lAngle.degToRad, inX, inY, 0)
            self *= rotation
        }
    }

}
