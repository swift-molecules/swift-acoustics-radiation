import Acoustics_Radiation
import Numeric
import Testing

extension `Acoustics.Radiation Tests` {
    @Suite struct Laws {
        static let radius = 0.1
        /// 120 triangles.
        static let mesh = `Acoustics.Radiation Tests`.sphere(radius: radius, bands: 6, around: 12)
        /// An asymmetric motion: the sphere's north side bulging while it rocks.
        static func motion(_ point: Point<3, Double>) -> Double {
            let z = point.z / radius
            let x = point.x / radius
            return 0.4 + z + 0.7 * x * z
        }
    }
}

extension `Acoustics.Radiation Tests`.Laws {
    @Test
    func `The power through a far sphere is the power the surface gives`() throws {
        let body = try Acoustics.Radiation.Body(Self.mesh)
        let field = try `Acoustics.Radiation Tests`.field(body, ka: 2, radius: Self.radius, velocity: Self.motion)
        // ∫ |p|² / 2ρ₀c r² dΩ, by the midpoint rule in θ and φ.
        let r = 50.0
        let polar = 24
        let around = 24
        var power = 0.0
        for i in 0..<polar {
            let theta = Double.pi * (Double(i) + 0.5) / Double(polar)
            for j in 0..<around {
                let phi = 2 * Double.pi * (Double(j) + 0.5) / Double(around)
                let direction = Vector<3, Double>(
                    x: Double.sin(theta) * Double.cos(phi),
                    y: Double.sin(theta) * Double.sin(phi),
                    z: Double.cos(theta)
                )
                let p = `Acoustics.Radiation Tests`.size(
                    field.pressure(toward: direction, at: Space.Physical.metres(r)).underlying
                )
                let solid = Double.sin(theta) * (Double.pi / Double(polar)) * (2 * Double.pi / Double(around))
                power += p * p / (2 * `Acoustics.Radiation Tests`.impedance) * r * r * solid
            }
        }
        #expect(
            abs(power - field.power.underlying) < 0.02 * field.power.underlying,
            "far \(power), surface \(field.power)"
        )
    }

    @Test
    func `Far away the field falls as one over the distance and agrees with the far-field formula`() throws {
        let body = try Acoustics.Radiation.Body(Self.mesh)
        let field = try `Acoustics.Radiation Tests`.field(body, ka: 1.5, radius: Self.radius, velocity: Self.motion)
        let direction = Vector<3, Double>(x: 0.3, y: -0.5, z: 0.8)
        let unit = direction.scaled(by: 1 / direction.length)
        let r = 200 * Self.radius
        let exact = field.pressure(at: Point(coordinates: unit.scaled(by: r))).underlying
        let far = field.pressure(toward: direction, at: Space.Physical.metres(r)).underlying
        #expect(`Acoustics.Radiation Tests`.error(far, exact) < 0.01, "far \(far), exact \(exact)")
    }

    @Test
    func `Pressure is linear in the velocity and the efficiency does not depend on its scale`() throws {
        let body = try Acoustics.Radiation.Body(Self.mesh)
        let one = try `Acoustics.Radiation Tests`.field(body, ka: 1, radius: Self.radius, velocity: Self.motion)
        let three = try `Acoustics.Radiation Tests`.field(body, ka: 1, radius: Self.radius) { 3 * Self.motion($0) }
        for (a, b) in zip(one.pressure, three.pressure) {
            #expect(`Acoustics.Radiation Tests`.error(b.underlying, a.underlying * Complex.Real(3)) < 1e-9)
        }
        #expect(abs(three.power.underlying - 9 * one.power.underlying) < 1e-9 * three.power.underlying)
        #expect(abs(three.efficiency - one.efficiency) < 1e-9)
    }

    @Test
    func `A body s times as large at 1 / s the frequency radiates alike`() throws {
        let small = try Acoustics.Radiation.Body(Self.mesh)
        let large = try Acoustics.Radiation.Body(
            Self.mesh.map { Point<3, Double>(coordinates: $0.coordinates.scaled(by: 2.5)) }
        )
        let a = try `Acoustics.Radiation Tests`.field(small, ka: 1.2, radius: Self.radius) { Self.motion($0) }
        let b = try `Acoustics.Radiation Tests`.field(large, ka: 1.2, radius: 2.5 * Self.radius) { point in
            Self.motion(Point(coordinates: point.coordinates.scaled(by: 1 / 2.5)))
        }
        #expect(abs(a.efficiency - b.efficiency) < 1e-9)
        for (p, q) in zip(a.pressure, b.pressure) {
            #expect(`Acoustics.Radiation Tests`.error(q.underlying, p.underlying) < 1e-9)
        }
    }

    @Test
    func `A surface listed inside out is turned`() throws {
        let body = try Acoustics.Radiation.Body(Self.mesh)
        let turned = try Acoustics.Radiation.Body(Self.mesh.reversed)
        let a = try `Acoustics.Radiation Tests`.field(body, ka: 1, radius: Self.radius, velocity: Self.motion)
        let b = try `Acoustics.Radiation Tests`.field(turned, ka: 1, radius: Self.radius, velocity: Self.motion)
        #expect(abs(a.efficiency - b.efficiency) < 1e-12)
        #expect(abs(body.area.underlying - turned.area.underlying) < 1e-15)
    }

    @Test
    func `Wrong input fails`() throws {
        let body = try Acoustics.Radiation.Body(Self.mesh)
        #expect(throws: Acoustics.Radiation.Error.velocity) {
            try Acoustics.Radiation.Field(body, frequency: Frequency.hertz(100), velocity: [Velocity.amount(1)])
        }
        let still = body.mesh.vertices.map { _ in Velocity.amount(1.0) }
        #expect(throws: Acoustics.Radiation.Error.frequency) {
            try Acoustics.Radiation.Field(body, frequency: Frequency.hertz(0), velocity: still)
        }
        let square = try Mesh(
            vertices: [
                Point<3, Double>(x: 0, y: 0, z: 0), Point(x: 1, y: 0, z: 0), Point(x: 1, y: 1, z: 0),
                Point(x: 0, y: 1, z: 0),
            ],
            faces: [Polygon(0, 1, 2, 3)]
        )
        #expect(throws: Acoustics.Radiation.Error.face(0)) { try Acoustics.Radiation.Body(square) }
        let below = try Mesh(
            vertices: [Point<3, Double>(x: 0, y: 0, z: 0), Point(x: 1, y: 0, z: 0), Point(x: 0, y: 1, z: -0.1)],
            faces: [Polygon(0, 1, 2)]
        )
        #expect(throws: Acoustics.Radiation.Error.vertex(2)) {
            try Acoustics.Radiation.Body(below, mounting: .baffled)
        }
    }
}
