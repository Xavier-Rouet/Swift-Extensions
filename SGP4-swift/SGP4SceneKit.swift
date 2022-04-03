//
//  SGP4SceneKit.swift
//  pxSat3D-SK
//
//  Created by Xav perso on 24/08/2020.
//  Copyright © 2020 P'tit Xav. All rights reserved.
//

import Foundation
import SceneKit

extension SCNVector3 {
	static func SCNVector3Make(v : SGP4Vector3) -> SCNVector3 {
		return SCNVector3(v.x, v.y, v.z)
	}
}
