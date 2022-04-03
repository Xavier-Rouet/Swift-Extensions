//
//  Array+extension.swift
//
//  Created by Xavier Rouet on 19/03/2022.
//

import Foundation

extension Array {
    /// Usage:
    /// ''''
    ///     var array1 = [1,2,3]
    ///     let array2 = [4,5,6]
    ///     let seven = 7
    ///     array1 += array2 // [1,2,3,4,5,6]
    ///     array1 += seven // [1,2,3,4,5,6,7]
    /// ''''
    static func +=(_ array: inout Array, _ element: Element) {
        array.append(element)
    }
    static func +=(_ array: inout Array, _ newElements: Array) {
        array.append(contentsOf: newElements)
    }

    /// Usage:
    /// ''''
    ///     var array1 = [1,2]
    ///     array1 *= 3    // [1,2,1,2,1,2]
    /// ''''
     static func *=(_ array: inout Array, _ nbTimes: UInt) {
        let save = array
        for _ in 0..<nbTimes {
            array.append(contentsOf: save)
        }
    }

    /// Usage:
    /// ''''
    ///     var array1 = [1,2,3,4,5]
    ///     array1--    // [1,2,3,4]
    ///     var a = array1-- // a=4, array1=[1,2,3]
    ///     --array1    // [2,3]
    ///     var a = --array1-- // a=2, array1=[3]
    /// ''''
    static prefix func --(_ array: inout Array) {
        array.removeFirst()
    }

    static prefix func --(_ array: inout Array) -> Element{
        array.removeFirst()
    }

    static postfix func --(_ array: inout Array) {
        array.removeLast()
    }

    static postfix func --(_ array: inout Array) -> Element {
        array.removeLast()
    }

}
