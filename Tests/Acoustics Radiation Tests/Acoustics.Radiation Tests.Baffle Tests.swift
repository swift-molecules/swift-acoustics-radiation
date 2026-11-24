import Acoustics_Plate
import Acoustics_Radiation
import Bessel
import Numeric
import Testing

extension `Acoustics.Radiation Tests` {
    @Suite struct Baffle {}
}

extension `Acoustics.Radiation Tests`.Baffle {
    /// A flat mesh in the plane z = 0.
    static func flat(_ mesh: Mesh<Point<2, Double>>) -> Mesh<Point<3, Double>> {
        mesh.map { Point<3, Double>(x: $0.x, y: $0.y, z: 0) }.triangulated
    }

    @Test(arguments: [1, Double.pi])
    func `A pulsating hemisphere on a baffle radiates as the whole sphere`(ka: Double) throws {
        let radius = 0.1
        // 112 triangles.
        let body = try Acoustics.Radiation.Body(
            `Acoustics.Radiation Tests`.sphere(radius: radius, bands: 8, around: 16, upper: true),
            mounting: .baffled
        )
        let a = `Acoustics.Radiation Tests`.equivalent(body)
        let field = try `Acoustics.Radiation Tests`.field(body, ka: ka * radius / a, radius: radius) { _ in 1 }
        let efficiency = ka * ka / (1 + ka * ka)
        #expect(
            abs(field.efficiency - efficiency) < 0.02 * efficiency,
            "found \(field.efficiency), exact \(efficiency)"
        )
        let surface = `Acoustics.Radiation Tests`.pulsating(ka: ka, radius: a, distance: a)
        let mean = `Acoustics.Radiation Tests`.mean(field)
        #expect(`Acoustics.Radiation Tests`.error(mean, surface) < 0.02, "surface \(mean), exact \(surface)")
        // Along the baffle as well as above it.
        for direction in [Vector<3, Double>(x: 0, y: 0, z: 1), Vector(x: 1, y: 0, z: 0), Vector(x: 1, y: 1, z: 1)] {
            let far = 100 * radius
            let found = field.pressure(toward: direction, at: Space.Physical.metres(far)).underlying
            let exact = `Acoustics.Radiation Tests`.pulsating(ka: ka, radius: a, distance: far)
            #expect(`Acoustics.Radiation Tests`.error(found, exact) < 0.03, "far: \(found), exact \(exact)")
        }
    }

    @Test(arguments: [0.5, 1, 2])
    func `A baffled disk moving as one radiates as a piston, 1 − J₁(2ka) / ka`(ka: Double) throws {
        // 160 triangles. The disk's outline is a polygon inside its circle: it is compared with the piston of its
        // own area.
        let body = try Acoustics.Radiation.Body(Self.flat(Mesh.disk(radius: 0.1, cells: 4)), mounting: .baffled)
        let radius = (body.area.underlying / Double.pi).squareRoot()
        let field = try `Acoustics.Radiation Tests`.field(body, ka: ka, radius: radius) { _ in 1 }
        let piston = 1 - 2 * Bessel.First<Double>(order: 1).value(at: 2 * ka) / (2 * ka)
        #expect(abs(field.efficiency - piston) < 0.01 * piston, "found \(field.efficiency), piston \(piston)")
    }

    @Test
    func `A simply supported panel's modes radiate as the plate's Rayleigh integral gives`() throws {
        let a = 0.4
        let b = 0.25
        let columns = 12
        let rows = 8
        let mesh = Mesh(columns: Cardinal(UInt(columns)), rows: Cardinal(UInt(rows))) { column, row in
            Point<2, Double>(
                x: a * Double(Int(bitPattern: column)) / Double(columns),
                y: b * Double(Int(bitPattern: row)) / Double(rows)
            )
        }
        let plate = Acoustics.Plate(
            mesh: mesh,
            thickness: Space.Physical.metres(0.003),
            material: Acoustics.Plate.Material(.sitkaSpruce, cut: .quarter),
            edge: .supported
        )
        let modes = try Acoustics.Plate.Modes(plate, count: 4)
        let body = try Acoustics.Radiation.Body(Self.flat(mesh), mounting: .baffled)
        for mode in modes {
            // The plate's shape is its deflection; the velocity is iω times it, and the efficiency does not
            // depend on the scale.
            let field = try Acoustics.Radiation.Field(
                body,
                frequency: mode.frequency,
                velocity: mode.shape.map { Velocity.amount($0) },
                air: `Acoustics.Radiation Tests`.air
            )
            #expect(
                abs(field.efficiency - mode.efficiency) < 0.03 * mode.efficiency,
                "at \(mode.frequency.underlying) Hz: found \(field.efficiency), plate \(mode.efficiency)"
            )
        }
    }
}
