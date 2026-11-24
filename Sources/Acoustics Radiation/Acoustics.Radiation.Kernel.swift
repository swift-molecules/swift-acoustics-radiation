public import Complex
internal import Galerkin
internal import Numeric

extension Acoustics.Radiation {
    /// The free-space Green's function of the Helmholtz equation, G = e^{ikR} / 4πR, for time e^{−iωt}, and its
    /// integrals over one face seen from one point.
    ///
    /// With r = x − y from a point y of the face to the point x that sees it, and R = |r|, its derivative along
    /// the face's normal is ∂G/∂n_y = e^{ikR} (1 − ikR) (r · n_y) / 4πR³.
    struct Kernel {
        /// The wavenumber, in radians per metre.
        let k: Double
    }
}

extension Acoustics.Radiation.Kernel {
    /// What a face adds to the equation of one point: ∫ G v dS and ∫ ∂G/∂n_y dS over it, with v the normal
    /// velocity, linear over the face.
    struct Integrals {
        var single: Complex.Number<Double> = .zero
        var double: Complex.Number<Double> = .zero

        static func + (lhs: Self, rhs: Self) -> Self {
            Self(single: lhs.single + rhs.single, double: lhs.double + rhs.double)
        }
    }

    /// A face as the integrals see it: its corners, its normal, and the velocity at each corner.
    struct Source {
        let corners: [Vector<3, Double>]
        let normal: Vector<3, Double>
        let velocity: [Double]
        /// The edges from the first corner to the second and to the third.
        let ab: Vector<3, Double>
        let ac: Vector<3, Double>
        /// Twice the area: the map's determinant from the reference triangle.
        let determinant: Double
        let centroid: Vector<3, Double>
        /// The longest edge.
        let size: Double

        init(corners: [Vector<3, Double>], normal: Vector<3, Double>, velocity: [Double]) {
            self.corners = corners
            self.normal = normal
            self.velocity = velocity
            self.ab = corners[1] - corners[0]
            self.ac = corners[2] - corners[0]
            self.determinant = ab.cross(ac).length
            self.centroid = (corners[0] + corners[1] + corners[2]).scaled(by: 1.0 / 3)
            self.size = Swift.max(ab.length, ac.length, (corners[2] - corners[1]).length)
        }

        /// Its mirror image in the plane z = 0.
        var image: Self {
            Self(
                corners: corners.map { Vector(x: $0.x, y: $0.y, z: -$0.z) },
                normal: Vector(x: normal.x, y: normal.y, z: -normal.z),
                velocity: velocity
            )
        }

        /// The point at (ξ, η) on the reference triangle.
        func point(_ xi: Double, _ eta: Double) -> Vector<3, Double> {
            corners[0] + ab.scaled(by: xi) + ac.scaled(by: eta)
        }

        /// The velocity at (ξ, η), linear between the corners.
        func velocity(_ xi: Double, _ eta: Double) -> Double {
            (1 - xi - eta) * velocity[0] + xi * velocity[1] + eta * velocity[2]
        }
    }

    /// The rule on each face, or on each piece of one, away from the point that sees it.
    static let rule = Galerkin.Reference.Triangle.quadrature(degree: 4, Double.self)
    /// The rule on a face seen from its own centroid, collapsed there.
    static let centred = Galerkin.Reference.Triangle.quadrature(degree: 8, around: Point(x: 1.0 / 3, y: 1.0 / 3))

    /// e^{ikR}.
    func phase(_ distance: Double) -> Complex.Number<Double> {
        Complex.Number(Double.cos(k * distance), Double.sin(k * distance))
    }

    /// The integrals over `source`, seen from x off it, by its rule on pieces small next to their distance from x:
    /// a piece is cut in four while x is nearer to its centroid than twice its longest edge.
    func integrals(from x: Vector<3, Double>, over source: Source) -> Integrals {
        var sums = Integrals()
        /// The rule on the piece with these corners (ξ, η) on the reference triangle, as a share of the whole.
        func rule(_ p: [Point<2, Double>], share: Double) {
            for (reference, weight) in zip(Self.rule.points, Self.rule.weights) {
                let xi = p[0].x + reference.x * (p[1].x - p[0].x) + reference.y * (p[2].x - p[0].x)
                let eta = p[0].y + reference.x * (p[1].y - p[0].y) + reference.y * (p[2].y - p[0].y)
                add(
                    from: x,
                    at: source.point(xi, eta),
                    velocity: source.velocity(xi, eta),
                    weight: weight * share * source.determinant,
                    normal: source.normal,
                    to: &sums
                )
            }
        }
        func piece(_ p: [Point<2, Double>], share: Double, depth: Int) {
            let centroid = source.point((p[0].x + p[1].x + p[2].x) / 3, (p[0].y + p[1].y + p[2].y) / 3)
            // Each cut halves the edges.
            let size = source.size / Double(1 << depth)
            guard depth < 8 && (x - centroid).length < 2 * size else { return rule(p, share: share) }
            func middle(_ i: Int, _ j: Int) -> Point<2, Double> {
                Point(x: (p[i].x + p[j].x) / 2, y: (p[i].y + p[j].y) / 2)
            }
            let m01 = middle(0, 1)
            let m12 = middle(1, 2)
            let m20 = middle(2, 0)
            piece([p[0], m01, m20], share: share / 4, depth: depth + 1)
            piece([m01, p[1], m12], share: share / 4, depth: depth + 1)
            piece([m20, m12, p[2]], share: share / 4, depth: depth + 1)
            piece([m12, m20, m01], share: share / 4, depth: depth + 1)
        }
        if (x - source.centroid).length < 2 * source.size {
            piece(Galerkin.Reference.Triangle.corners(Double.self), share: 1, depth: 0)
        } else {
            for (reference, weight) in zip(Self.rule.points, Self.rule.weights) {
                add(
                    from: x,
                    at: source.point(reference.x, reference.y),
                    velocity: source.velocity(reference.x, reference.y),
                    weight: weight * source.determinant,
                    normal: source.normal,
                    to: &sums
                )
            }
        }
        return sums
    }

    /// G v and ∂G/∂n_y at y, seen from x, times `weight`, added to `sums`.
    func add(
        from x: Vector<3, Double>,
        at y: Vector<3, Double>,
        velocity v: Double,
        weight: Double,
        normal: Vector<3, Double>,
        to sums: inout Integrals
    ) {
        let r = x - y
        let distance = r.length
        let g = phase(distance) * Complex.Real(weight / (4 * Double.pi * distance))
        sums.single += g * Complex.Real(v)
        let slope = Complex.Number(1, -k * distance) * Complex.Real(r.dot(normal) / (distance * distance))
        sums.double += g * slope
    }

    /// The integrals over `source` seen from its own centroid.
    ///
    /// The face is flat, so r · n vanishes on it and ∫ ∂G/∂n_y is zero. Of ∫ G v, the part v̄ / 4πR, with v̄ the
    /// velocity at the centroid, is integrated in closed form over the triangle, and the bounded rest by the rule
    /// collapsed at the centroid.
    func integrals(centredOn source: Source) -> Integrals {
        let x = source.centroid
        let mean = (source.velocity[0] + source.velocity[1] + source.velocity[2]) / 3
        let fourPi = 4 * Double.pi
        var single = Complex.Number<Double>(mean * Self.inverse(over: source.corners, from: x) / fourPi)
        for (reference, weight) in zip(Self.centred.points, Self.centred.weights) {
            let distance = (x - source.point(reference.x, reference.y)).length
            let w = weight * source.determinant
            let v = source.velocity(reference.x, reference.y)
            let e = phase(distance)
            // G v − v̄ / 4πR = (e^{ikR} − 1) v̄ / 4πR + e^{ikR} (v − v̄) / 4πR, bounded.
            single += (e - Complex.Number<Double>.one) * Complex.Real(w * mean / (fourPi * distance))
            single += e * Complex.Real(w * (v - mean) / (fourPi * distance))
        }
        return Integrals(single: single, double: .zero)
    }

    /// ∫ 1/R dA over a flat polygon from a point inside it, in closed form: each edge at distance h from the point,
    /// with its ends at signed positions s₁ < s₂ along it from the foot of the perpendicular, adds
    /// h (asinh(s₂/h) − asinh(s₁/h)).
    static func inverse(over corners: [Vector<3, Double>], from x: Vector<3, Double>) -> Double {
        var sum = 0.0
        for edge in 0..<corners.count {
            let a = corners[edge]
            let b = corners[(edge + 1) % corners.count]
            let tangent = (b - a).scaled(by: 1 / (b - a).length)
            let s1 = (a - x).dot(tangent)
            let s2 = (b - x).dot(tangent)
            let h = ((a - x) - tangent.scaled(by: s1)).length
            sum += h * (asinh(s2 / h) - asinh(s1 / h))
        }
        return sum
    }

    static func asinh(_ value: Double) -> Double {
        let magnitude = Swift.abs(value)
        let result = Double.log(magnitude + (magnitude * magnitude + 1).squareRoot())
        return value < 0 ? -result : result
    }
}
