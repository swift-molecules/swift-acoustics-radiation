import Acoustics_Radiation
import Numeric
import Testing

@Suite struct `Acoustics.Radiation Tests` {}

extension `Acoustics.Radiation Tests` {
    /// A sphere of `radius` about the origin, or its upper half, cut into `bands` bands of latitude from the north
    /// pole to the south pole (or to the equator) and `around` meridians, as triangles.
    static func sphere(radius: Double, bands: Int, around: Int, upper: Bool = false) -> Mesh<Point<3, Double>> {
        let last = upper ? bands / 2 : bands
        func point(_ band: Int, _ meridian: Int) -> Point<3, Double> {
            let polar = Double.pi * Double(band) / Double(bands)
            let longitude = 2 * Double.pi * Double(meridian % around) / Double(around)
            if band == 0 { return Point(x: 0, y: 0, z: radius) }
            if band == bands { return Point(x: 0, y: 0, z: -radius) }
            let z = band * 2 == bands ? 0 : radius * Double.cos(polar)
            return Point(
                x: radius * Double.sin(polar) * Double.cos(longitude),
                y: radius * Double.sin(polar) * Double.sin(longitude),
                z: z
            )
        }
        var polygons: [Polygon<Point<3, Double>>] = []
        for band in 0..<last {
            for meridian in 0..<around {
                let a = point(band, meridian)
                let b = point(band + 1, meridian)
                let c = point(band + 1, meridian + 1)
                let d = point(band, meridian + 1)
                // Counterclockwise seen from outside: north to south, then east.
                if band != 0 { polygons.append(Polygon(a, b, d)) }
                if band != bands - 1 { polygons.append(Polygon(b, c, d)) }
            }
        }
        return Mesh(polygons).welded
    }

    /// The polar angle of a point.
    static func polar(_ point: Point<3, Double>) -> Double {
        let r = (point.x * point.x + point.y * point.y + point.z * point.z).squareRoot()
        return Double.acos(Swift.max(-1, Swift.min(1, point.z / r)))
    }

    /// The field of `body` moving with normal velocity `velocity` (in metres per second at each vertex) at the
    /// frequency where k a = `ka`.
    static func field(
        _ body: Acoustics.Radiation.Body,
        ka: Double,
        radius: Double,
        velocity: (Point<3, Double>) -> Double
    ) throws -> Acoustics.Radiation.Field {
        try Acoustics.Radiation.Field(
            body,
            frequency: Frequency.hertz(ka / radius * air.speed.underlying / (2 * Double.pi)),
            velocity: body.mesh.vertices.map { Velocity.amount(velocity($0)) },
            air: air
        )
    }

    static let air = Acoustics.Air.room

    /// ρ₀ c.
    static var impedance: Double { air.density.underlying * air.speed.underlying }

    /// |z|.
    static func size(_ z: Complex.Number<Double>) -> Double {
        (z.realValue * z.realValue + z.imaginaryValue * z.imaginaryValue).squareRoot()
    }

    /// The relative difference |found − exact| / |exact|.
    static func error(_ found: Complex.Number<Double>, _ exact: Complex.Number<Double>) -> Double {
        size(found - exact) / size(exact)
    }

    /// The radius of the sphere with the volume the surface encloses (with its mirror image, on a baffle).
    ///
    /// A sphere's faces are flat, so its surface lies inside the sphere it was drawn on; it radiates as the sphere
    /// of equal volume, which is what the closed forms are taken at.
    static func equivalent(_ body: Acoustics.Radiation.Body) -> Double {
        var volume = 0.0
        for face in body.mesh.faces {
            let c = face.vertices.map { body.mesh.vertices[Int(bitPattern: $0)].coordinates }
            volume += c[0].dot(c[1].cross(c[2])) / 6
        }
        if body.mounting == .baffled { volume *= 2 }
        return Double.exp(Double.log(3 * volume / (4 * Double.pi)) / 3)
    }

    /// The area-weighted mean of the pressure on the faces.
    static func mean(_ field: Acoustics.Radiation.Field) -> Complex.Number<Double> {
        var sum = Complex.Number<Double>.zero
        var area = 0.0
        for (face, pressure) in zip(field.body.mesh.faces, field.pressure) {
            let corners = face.vertices.map { field.body.mesh.vertices[Int(bitPattern: $0)].coordinates }
            let weight = (corners[1] - corners[0]).cross(corners[2] - corners[0]).length / 2
            sum += pressure.underlying * Complex.Real(weight)
            area += weight
        }
        return sum / Complex.Real(area)
    }

    /// The pressure of a sphere of radius a pulsating with velocity v, at distance r ≥ a from its centre:
    /// ρ₀ c v (ika / (ika − 1)) (a / r) e^{ik(r − a)}.
    static func pulsating(
        ka: Double,
        radius a: Double,
        distance r: Double,
        velocity v: Double = 1
    )
        -> Complex.Number<Double>
    {
        let k = ka / a
        let phase = Complex.Number(Double.cos(k * (r - a)), Double.sin(k * (r - a)))
        return Complex.Number(0, ka) / Complex.Number(-1, ka) * phase * Complex.Real(impedance * v * a / r)
    }
}
