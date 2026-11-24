public import Mesh

extension Acoustics.Radiation {
    /// The surface of a body that sends sound into the air, as flat triangles, and where that air is.
    public struct Body {
        /// The surface: its vertices in metres, its faces triangles that run counterclockwise seen from the air.
        public let mesh: Mesh<Point<3, Double>>
        /// Whether the air is all around the body or above a rigid plane.
        public let mounting: Mounting
        /// Each face's corners and the quantities the integrals need.
        let faces: [Face]
        /// Points inside the body, where the radiated field's integral must vanish.
        let interior: [Vector<3, Double>]

        /// The body with surface `mesh`, its vertices in metres and its faces triangles, in `mounting`.
        ///
        /// A face that is not a triangle, or has no area, fails with `.face(index)`. A surface listed inside out (its
        /// faces counterclockwise seen from within) is turned. In `.baffled`, a vertex below the plane z = 0 fails
        /// with `.vertex(index)`.
        public init(_ mesh: Mesh<Point<3, Double>>, mounting: Mounting = .free) throws(Acoustics.Radiation.Error) {
            var faces: [Face] = []
            faces.reserveCapacity(mesh.faces.count)
            for (index, polygon) in mesh.faces.enumerated() {
                let vertices = polygon.vertices.map { Int(bitPattern: $0) }
                guard vertices.count == 3 else { throw .face(index) }
                let face = Face(
                    vertices: vertices,
                    corners: vertices.map { mesh.vertices[$0].coordinates }
                )
                guard face.area > 0 else { throw .face(index) }
                faces.append(face)
            }
            if mounting == .baffled {
                for (index, vertex) in mesh.vertices.enumerated() where vertex.z < 0 {
                    throw .vertex(index)
                }
            }
            // The volume the faces enclose, with the baffle's plane closing it, and their area projected on it.
            var volume = 0.0
            var projection = 0.0
            for face in faces {
                volume += face.centroid.dot(face.normal) * face.area / 3
                projection += face.normal.z * face.area
            }
            let size = faces.reduce(0.0) { $0 + $1.area }
            let inward =
                Swift.abs(volume) > 1e-9 * size * size.squareRoot()
                ? volume < 0 : mounting == .baffled && projection < 0
            self.mesh = inward ? mesh.reversed : mesh
            self.mounting = mounting
            self.faces = inward ? faces.map(\.reversed) : faces
            self.interior = Self.interior(self.faces, baffled: mounting == .baffled)
        }
    }
}

extension Acoustics.Radiation.Body {
    /// Where the air is.
    public enum Mounting: Hashable, Sendable {
        /// All around the body, without bound.
        case free
        /// Above an infinite rigid plane, z = 0, the body on or above it: the air fills z > 0.
        case baffled
    }
}

extension Acoustics.Radiation.Body {
    /// The total area of the surface's faces.
    public var area: Area<Space.Physical>.Value<Double> {
        Area<Space.Physical>.amount(faces.reduce(0.0) { $0 + $1.area })
    }
}

extension Acoustics.Radiation.Body {
    /// Points inside the body, spread through it: from evenly spaced faces, a point part of the way along the
    /// inward normal to where it next meets the surface (or the baffle). The share of the way is spread over
    /// 0.2 to 0.8 by the golden ratio, so no two points sit alike. A face whose inward ray leaves the body, as
    /// every face of a flat board does, gives no point.
    static func interior(_ faces: [Face], baffled: Bool) -> [Vector<3, Double>] {
        let count = Swift.min(faces.count, Swift.max(8, faces.count / 8))
        var points: [Vector<3, Double>] = []
        for step in 0..<count {
            let face = faces[step * faces.count / count]
            let direction = Vector<3, Double>.zero - face.normal
            var reach = Double.infinity
            for other in faces {
                if let distance = other.hit(from: face.centroid, along: direction), distance > 1e-9 * face.diameter {
                    reach = Swift.min(reach, distance)
                }
            }
            if baffled && direction.z < 0 { reach = Swift.min(reach, face.centroid.z / -direction.z) }
            guard reach.isFinite, reach > 0 else { continue }
            let golden = Double(step) * 0.618_033_988_749_894_8
            let share = 0.2 + 0.6 * (golden - golden.rounded(.down))
            points.append(face.centroid + direction.scaled(by: share * reach))
        }
        return points
    }
}

extension Acoustics.Radiation.Body {
    /// A flat triangle of the surface.
    struct Face {
        /// The mesh's vertices at its corners.
        let vertices: [Int]
        /// Its corners in metres, counterclockwise seen from the air.
        let corners: [Vector<3, Double>]
        /// The unit normal, into the air.
        let normal: Vector<3, Double>
        let area: Double
        let centroid: Vector<3, Double>
        /// Its longest edge.
        let diameter: Double
        /// Whether it lies in the baffle's plane.
        let planar: Bool

        init(vertices: [Int], corners: [Vector<3, Double>]) {
            self.vertices = vertices
            self.corners = corners
            let cross = (corners[1] - corners[0]).cross(corners[2] - corners[0])
            let length = cross.length
            self.area = length / 2
            self.normal = length > 0 ? cross.scaled(by: 1 / length) : .zero
            self.centroid = (corners[0] + corners[1] + corners[2]).scaled(by: 1.0 / 3)
            self.diameter = Swift.max(
                (corners[1] - corners[0]).length,
                (corners[2] - corners[1]).length,
                (corners[0] - corners[2]).length
            )
            self.planar = corners.allSatisfy { $0.z == 0 }
        }

        /// How far from `origin` along the unit `direction` the ray meets this face, if it does (Möller and
        /// Trumbore, 1997).
        func hit(from origin: Vector<3, Double>, along direction: Vector<3, Double>) -> Double? {
            let ab = corners[1] - corners[0]
            let ac = corners[2] - corners[0]
            let p = direction.cross(ac)
            let determinant = ab.dot(p)
            guard Swift.abs(determinant) > 1e-14 * ab.length * ac.length else { return nil }
            let offset = origin - corners[0]
            let u = offset.dot(p) / determinant
            guard u >= 0, u <= 1 else { return nil }
            let q = offset.cross(ab)
            let v = direction.dot(q) / determinant
            guard v >= 0, u + v <= 1 else { return nil }
            let distance = ac.dot(q) / determinant
            return distance > 0 ? distance : nil
        }

        var reversed: Self {
            Self(vertices: [vertices[0], vertices[2], vertices[1]], corners: [corners[0], corners[2], corners[1]])
        }
    }
}

extension Acoustics.Radiation.Body: Sendable {}
extension Acoustics.Radiation.Body.Face: Sendable {}
