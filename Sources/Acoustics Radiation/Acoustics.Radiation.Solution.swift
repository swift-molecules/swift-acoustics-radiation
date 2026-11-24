public import Complex

extension Acoustics.Radiation {
    /// The pressure on every face of a body whose surface has normal velocity v, per unit of iωρ₀: the solution
    /// of the boundary integral equation with ∂p/∂n = v.
    ///
    /// Each face carries a constant pressure, the velocity is linear over it, and the Kirchhoff–Helmholtz equation
    /// is collocated at each face's centroid, where the surface is smooth:
    /// ½ p(x) = ∫ p ∂G/∂n_y dS − ∫ G v dS.
    ///
    /// Alone it fails at the frequencies where the body's inside, held at zero pressure on its surface, would
    /// resonate. So the same integral is also asked to vanish at points inside the body, as it must (Schenck's
    /// CHIEF, 1968), and the equations are solved together in least squares. A point inside fails to pin a
    /// frequency only if it lies on a nodal surface of that resonance, so there are several, spread through the
    /// inside.
    ///
    /// In a baffle, the image of every face in the plane z = 0 is added to it, which makes the plane rigid. A face
    /// that lies in the plane is its own image: its ∂G/∂n_y cancels its image's, and a point on it sees the whole
    /// half-space, so its free term is p(x), not ½ p(x). For a flat board that is Rayleigh's integral.
    struct Solution {
        let pressure: [Complex.Number<Double>]

        init(_ body: Body, wavenumber k: Double, velocity: [Velocity.Value<Double>]) throws(Acoustics.Radiation.Error) {
            let kernel = Kernel(k: k)
            let n = body.faces.count
            let sources = body.faces.map { face in
                Kernel.Source(
                    corners: face.corners,
                    normal: face.normal,
                    velocity: face.vertices.map { velocity[$0].underlying }
                )
            }
            let baffled = body.mounting == .baffled
            let images = baffled ? sources.map(\.image) : []
            /// What every face adds to the equation of `x`, which is the centroid of face `own` if not nil.
            func row(_ x: Vector<3, Double>, own: Int?) -> [Kernel.Integrals] {
                (0..<n).map { j in
                    let sums =
                        j == own ? kernel.integrals(centredOn: sources[j]) : kernel.integrals(from: x, over: sources[j])
                    guard baffled else { return sums }
                    guard !body.faces[j].planar else {
                        return Kernel.Integrals(single: sums.single + sums.single, double: .zero)
                    }
                    return sums + kernel.integrals(from: x, over: images[j])
                }
            }
            var matrix = Dense(rows: n + body.interior.count, columns: n)
            var right = [Complex.Number<Double>](repeating: .zero, count: n + body.interior.count)
            for i in 0..<n {
                let sums = row(body.faces[i].centroid, own: i)
                for j in 0..<n {
                    matrix[i, j] = -sums[j].double
                    right[i] -= sums[j].single
                }
                matrix[i, i] += Complex.Number(baffled && body.faces[i].planar ? 1 : 0.5)
            }
            for (index, point) in body.interior.enumerated() {
                let sums = row(point, own: nil)
                for j in 0..<n {
                    matrix[n + index, j] = sums[j].double
                    right[n + index] += sums[j].single
                }
            }
            guard let pressure = matrix.solve(right) else { throw .singular }
            self.pressure = pressure
        }
    }
}
