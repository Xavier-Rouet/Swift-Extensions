//
//  String+extension.swift
//
//  Created by Xavier Rouet on 19/03/2022.
//

import Foundation

extension String {
    /// Usage:
    /// ''''
    ///     print("The answer is %d" % 42)
    ///     print("The answer is %@" % "42")
    ///     print("The %@ is %d" % ["answer",42])
    ///     print("The %@ is %@" % ["answer","42"])
    /// ''''
    static func %(_ key: String, _ arg: CVarArg) -> String {
        String(format: NSLocalizedString(key, comment: key), arg)
    }

    static func %(_ key: String, _ arg: [CVarArg]) -> String {
        String(format: NSLocalizedString(key, comment: key), arguments: arg)
    }
}
