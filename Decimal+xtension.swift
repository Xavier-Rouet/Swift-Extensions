//
//  Decimal+extension.swift
//
//  Created by Xavier Rouet on 22/12/2021.
//  Copyright © 2021 P'tit Xav. All rights reserved.
//

import Foundation

extension Float {
    public static let π : Float = pi
    public static let halfπ = 0.5 * π
    public static let twoπ  = 2.0 * π
    var degToRad : Float {self * .π / 180.0}
    var radToDeg : Float {self * 180.0 / .π}
    var radToHour : Float {self * 12.0 / .π}
}

prefix operator √
prefix func √(_ a : Float) -> Float {
    return sqrt(a)
}

extension Double {
    public static let π : Double = pi
    public static let halfπ = 0.5 * π
    public static let twoπ  = 2.0 * π
    var degToRad : Double {self * .π / 180.0}
    var radToDeg : Double {self * 180.0 / .π}
    var radToHour : Double {self * 12.0 / .π}
}

prefix func √(_ a : Double) -> Double {
    return sqrt(a)
}

