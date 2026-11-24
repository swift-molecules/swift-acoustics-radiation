// swift-tools-version: 6.4

import PackageDescription

let package = Package(
    name: "swift-acoustics-radiation",
    platforms: [
        .macOS(.v27),
        .iOS(.v27),
        .tvOS(.v27),
        .watchOS(.v27),
        .visionOS(.v27),
    ],
    products: [
        .library(name: "Acoustics Radiation", targets: ["Acoustics Radiation"])
    ],
    dependencies: [
        .package(url: "https://github.com/swift-molecules/swift-acoustics.git", branch: "main"),
        .package(url: "https://github.com/swift-molecules/swift-galerkin.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-mesh.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-complex.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-numeric.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-frequency.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-pressure.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-power.git", branch: "main"),
        .package(url: "https://github.com/swift-molecules/swift-acoustics-plate.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-bessel.git", branch: "main"),
    ],
    targets: [
        .target(
            name: "Acoustics Radiation",
            dependencies: [
                .product(name: "Acoustics", package: "swift-acoustics"),
                .product(name: "Galerkin", package: "swift-galerkin"),
                .product(name: "Mesh", package: "swift-mesh"),
                .product(name: "Complex", package: "swift-complex"),
                .product(name: "Numeric", package: "swift-numeric"),
                .product(name: "Frequency", package: "swift-frequency"),
                .product(name: "Pressure", package: "swift-pressure"),
                .product(name: "Power", package: "swift-power"),
            ]
        ),
        .testTarget(
            name: "Acoustics Radiation Tests",
            dependencies: [
                .target(name: "Acoustics Radiation"),
                .product(name: "Acoustics Plate", package: "swift-acoustics-plate"),
                .product(name: "Bessel", package: "swift-bessel"),
                .product(name: "Numeric", package: "swift-numeric"),
            ]
        ),
    ],
    swiftLanguageModes: [.v6]
)

for target in package.targets where ![.system, .binary, .plugin, .macro].contains(target.type) {
    target.swiftSettings =
        (target.swiftSettings ?? []) + [
            .strictMemorySafety(),
            .enableUpcomingFeature("ExistentialAny"),
            .enableUpcomingFeature("InternalImportsByDefault"),
            .enableUpcomingFeature("MemberImportVisibility"),
            .enableUpcomingFeature("NonisolatedNonsendingByDefault"),
            .enableUpcomingFeature("InferIsolatedConformances"),
            .enableExperimentalFeature("Lifetimes"),
            .treatAllWarnings(as: .error),
        ]
}
