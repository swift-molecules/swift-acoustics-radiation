import Acoustics_Radiation
import Numeric
import Testing

extension `Acoustics.Radiation Tests` {
    @Suite struct Sphere {
        static let radius = 0.1
        /// 224 triangles.
        static let body = try! Acoustics.Radiation.Body(
            `Acoustics.Radiation Tests`.sphere(radius: radius, bands: 8, around: 16)
        )
        /// The radius of the sphere with the body's volume.
        static let equivalent = `Acoustics.Radiation Tests`.equivalent(body)
    }
}

extension `Acoustics.Radiation Tests`.Sphere {
    // ka = π is the lowest frequency at which the inside of the sphere, held at zero pressure, resonates; there the
    // Kirchhoff–Helmholtz equation alone has no unique solution.
    @Test(arguments: [0.5, 2, Double.pi])
    func `A pulsating sphere radiates with efficiency (ka)² / (1 + (ka)²) and its closed-form pressure`(
        ka: Double
    ) throws {
        let a = Self.equivalent
        let field = try `Acoustics.Radiation Tests`.field(Self.body, ka: ka * Self.radius / a, radius: Self.radius) {
            _ in 1
        }
        let efficiency = ka * ka / (1 + ka * ka)
        #expect(
            abs(field.efficiency - efficiency) < 0.02 * efficiency,
            "found \(field.efficiency), exact \(efficiency)"
        )
        let surface = `Acoustics.Radiation Tests`.pulsating(ka: ka, radius: a, distance: a)
        let mean = `Acoustics.Radiation Tests`.mean(field)
        #expect(`Acoustics.Radiation Tests`.error(mean, surface) < 0.02, "surface \(mean), exact \(surface)")
        // Near and far, in three directions.
        for direction in [Vector<3, Double>(x: 0, y: 0, z: 1), Vector(x: 1, y: 2, z: -1), Vector(x: -1, y: 0, z: 0)] {
            let unit = direction.scaled(by: 1 / direction.length)
            let near = 2 * Self.radius
            let found = field.pressure(at: Point(coordinates: unit.scaled(by: near))).underlying
            let exact = `Acoustics.Radiation Tests`.pulsating(ka: ka, radius: a, distance: near)
            #expect(`Acoustics.Radiation Tests`.error(found, exact) < 0.03, "at 2a: \(found), exact \(exact)")
            let far = 100 * Self.radius
            let distant = field.pressure(toward: direction, at: Space.Physical.metres(far)).underlying
            let farExact = `Acoustics.Radiation Tests`.pulsating(ka: ka, radius: a, distance: far)
            #expect(`Acoustics.Radiation Tests`.error(distant, farExact) < 0.03, "far: \(distant), exact \(farExact)")
        }
    }

    // ka = 4.4934, the first zero of j₁, is where the inside resonates with the oscillating sphere's pattern.
    @Test(arguments: [1, 4.493_409_457_909_064])
    func `An oscillating sphere radiates with efficiency (ka)⁴ / (4 + (ka)⁴) as a dipole`(ka: Double) throws {
        let a = Self.equivalent
        let field = try `Acoustics.Radiation Tests`.field(Self.body, ka: ka * Self.radius / a, radius: Self.radius) {
            point in point.z / Self.radius
        }
        let k4 = ka * ka * ka * ka
        let efficiency = k4 / (4 + k4)
        #expect(
            abs(field.efficiency - efficiency) < 0.025 * efficiency,
            "found \(field.efficiency), exact \(efficiency)"
        )
        // Far away the pressure goes as cos θ from the axis of motion.
        let distance = Space.Physical.metres(100 * Self.radius)
        let axis = `Acoustics.Radiation Tests`.size(
            field.pressure(toward: Vector(x: 0, y: 0, z: 1), at: distance).underlying
        )
        for angle in [Double.pi / 6, Double.pi / 3, Double.pi / 2, 2 * Double.pi / 3] {
            let direction = Vector<3, Double>(x: Double.sin(angle), y: 0, z: Double.cos(angle))
            let found = `Acoustics.Radiation Tests`.size(field.pressure(toward: direction, at: distance).underlying)
            #expect(abs(found / axis - abs(Double.cos(angle))) < 0.02, "at \(angle): \(found / axis)")
        }
    }

    @Test
    func `A finer sphere comes closer to the closed form, at second order`() throws {
        let ka = 2.0
        var errors: [Double] = []
        for (bands, around) in [(6, 12), (8, 16)] {
            let body = try Acoustics.Radiation.Body(
                `Acoustics.Radiation Tests`.sphere(radius: Self.radius, bands: bands, around: around)
            )
            let a = `Acoustics.Radiation Tests`.equivalent(body)
            let field = try `Acoustics.Radiation Tests`.field(body, ka: ka * Self.radius / a, radius: Self.radius) {
                point in point.z / Self.radius
            }
            errors.append(abs(field.efficiency - 16.0 / 20))
        }
        // From 120 faces to 224 the faces' area falls by 1.87, and so does the error at second order.
        #expect(errors[1] < errors[0] / 1.6, "\(errors)")
    }
}
