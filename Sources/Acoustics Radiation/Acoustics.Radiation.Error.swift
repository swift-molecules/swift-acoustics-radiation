extension Acoustics.Radiation {
    /// Why a body or its radiation could not be made.
    public enum Error: Swift.Error, Hashable, Sendable {
        /// The face at this index is not a triangle with area.
        case face(Int)
        /// The vertex at this index lies below a baffle.
        case vertex(Int)
        /// The velocities are not one for each vertex of the body.
        case velocity
        /// The frequency is not positive and finite.
        case frequency
        /// The boundary equations have no unique solution, as when faces overlap.
        case singular
    }
}
