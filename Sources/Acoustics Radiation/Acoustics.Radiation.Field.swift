public import Complex
public import Frequency
internal import Galerkin
internal import Numeric
public import Power
public import Pressure

extension Acoustics.Radiation {
    /// The sound a body sends into the air when its surface moves at one frequency: the pressure on every face,
    /// the power it radiates and how well it radiates, and the pressure anywhere in the air.
    ///
    /// Amplitudes are complex, for time e^{−iωt}: a value A stands for Re(A e^{−iωt}).
    public struct Field {
        /// The body.
        public let body: Body
        /// The air around it.
        public let air: Acoustics.Air
        /// The frequency of the motion.
        public let frequency: Frequency.Value<Double>
        /// The velocity amplitude at each vertex of the body's mesh, along the surface's normal into the air, and
        /// linear over each face.
        public let velocity: [Velocity.Value<Double>]
        /// The pressure amplitude on each face, in the mesh's order.
        public let pressure: [Pressure.Value<Complex.Number<Double>>]
        /// The power radiated into the air, averaged over a period.
        public let power: Power.Value<Double>
        /// The radiation efficiency: `power` over ρ₀ c times the body's area times the mean square normal velocity,
        /// ½ ∫ |v|² dS / S.
        public let efficiency: Double

        /// The field of `body` moving at `frequency` with normal velocity amplitude `velocity` at each vertex of its
        /// mesh, into `air`.
        ///
        /// It fails with `.velocity` unless there is one velocity per vertex, with `.frequency` unless the frequency
        /// is positive and finite, and with `.singular` if the boundary equations cannot be solved.
        public init(
            _ body: Body,
            frequency: Frequency.Value<Double>,
            velocity: [Velocity.Value<Double>],
            air: Acoustics.Air = .room
        ) throws(Acoustics.Radiation.Error) {
            guard velocity.count == body.mesh.vertices.count else { throw .velocity }
            let f = frequency.underlying
            guard f > 0, f.isFinite else { throw .frequency }
            let solution = try Solution(body, wavenumber: 2 * Double.pi * f / air.speed.underlying, velocity: velocity)
            let rho = air.density.underlying
            let omega = 2 * Double.pi * f
            self.body = body
            self.air = air
            self.frequency = frequency
            self.velocity = velocity
            // p = iωρ₀ times the solution for unit ∂p/∂n per (m/s), which is what `Solution` solves for.
            let scale = Complex.Number<Double>(0, omega * rho)
            let pressure = solution.pressure.map { $0 * scale }
            self.pressure = pressure.map { Pressure.pascals($0) }
            // W = ½ Re Σ pⱼ ∫ⱼ v dS, and ∫ |v|² dS with v linear over each face.
            var flux = 0.0
            var square = 0.0
            for (index, face) in body.faces.enumerated() {
                let v = face.vertices.map { velocity[$0].underlying }
                let sum = v[0] + v[1] + v[2]
                flux += pressure[index].realValue * face.area * sum / 3
                square += face.area / 12 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + sum * sum)
            }
            let watts = flux / 2
            self.power = Power.watts(watts)
            self.efficiency = square > 0 ? 2 * watts / (rho * air.speed.underlying * square) : 0
        }
    }
}

extension Acoustics.Radiation.Field {
    /// The wavenumber ω / c.
    var wavenumber: Double { 2 * Double.pi * frequency.underlying / air.speed.underlying }

    /// The faces as sources of the integrals: corners, normal and velocity, with their images where baffled.
    var sources: [(face: Acoustics.Radiation.Kernel.Source, pressure: Complex.Number<Double>)] {
        var sources: [(face: Acoustics.Radiation.Kernel.Source, pressure: Complex.Number<Double>)] = []
        for (index, face) in body.faces.enumerated() {
            let source = Acoustics.Radiation.Kernel.Source(
                corners: face.corners,
                normal: face.normal,
                velocity: face.vertices.map { velocity[$0].underlying }
            )
            sources.append((source, pressure[index].underlying))
            if body.mounting == .baffled { sources.append((source.image, pressure[index].underlying)) }
        }
        return sources
    }

    /// The pressure amplitude at `point` in the air, off the surface, by the Kirchhoff–Helmholtz integral
    /// p(x) = ∫ (p ∂G/∂n_y − G ∂p/∂n) dS, with ∂p/∂n = iωρ₀ v.
    public func pressure(at point: Point<3, Double>) -> Pressure.Value<Complex.Number<Double>> {
        let kernel = Acoustics.Radiation.Kernel(k: wavenumber)
        let x = point.coordinates
        let iOmegaRho = Complex.Number<Double>(0, 2 * Double.pi * frequency.underlying * air.density.underlying)
        var sum = Complex.Number<Double>.zero
        for (source, p) in sources {
            let integrals = kernel.integrals(from: x, over: source)
            sum += p * integrals.double - iOmegaRho * integrals.single
        }
        return Pressure.pascals(sum)
    }

    /// The pressure amplitude far away, at `distance` along `direction` from the origin, where it falls as one
    /// over the distance: e^{ikr} / 4πr ∫ (−ik (r̂ · n) p − iωρ₀ v) e^{−ik r̂ · y} dS. The direction need not be of
    /// unit length; it must not be zero.
    public func pressure(
        toward direction: Vector<3, Double>,
        at distance: Length<Space.Physical, Double>
    ) -> Pressure.Value<Complex.Number<Double>> {
        let length = direction.length
        precondition(length > 0 && length.isFinite, "a direction")
        let unit = direction.scaled(by: 1 / length)
        let k = wavenumber
        let kernel = Acoustics.Radiation.Kernel(k: k)
        let omegaRho = 2 * Double.pi * frequency.underlying * air.density.underlying
        var sum = Complex.Number<Double>.zero
        for (source, p) in sources {
            let tilt = Complex.Number<Double>(0, -k * unit.dot(source.normal)) * p
            for (reference, weight) in zip(
                Acoustics.Radiation.Kernel.rule.points,
                Acoustics.Radiation.Kernel.rule.weights
            ) {
                let v = source.velocity(reference.x, reference.y)
                let integrand = tilt - Complex.Number<Double>(0, omegaRho * v)
                let y = source.point(reference.x, reference.y)
                sum += integrand * kernel.phase(-unit.dot(y)) * Complex.Real(weight * source.determinant)
            }
        }
        let r = distance.underlying
        return Pressure.pascals(sum * kernel.phase(r) * Complex.Real(1 / (4 * Double.pi * r)))
    }
}

extension Acoustics.Radiation.Field: Sendable {}
