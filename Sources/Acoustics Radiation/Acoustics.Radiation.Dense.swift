public import Complex

extension Acoustics.Radiation {
    /// A dense complex matrix with at least as many rows as columns.
    ///
    /// The boundary element method couples every face to every other, so its matrix is full.
    struct Dense {
        let rows: Int
        let columns: Int
        /// Row by row.
        var entries: [Complex.Number<Double>]

        init(rows: Int, columns: Int) {
            precondition(rows >= columns, "at least as many equations as unknowns")
            self.rows = rows
            self.columns = columns
            self.entries = [Complex.Number<Double>](repeating: .zero, count: rows * columns)
        }

        subscript(row: Int, column: Int) -> Complex.Number<Double> {
            get { entries[row * columns + column] }
            set { entries[row * columns + column] = newValue }
        }
    }
}

extension Acoustics.Radiation.Dense {
    /// The x that minimizes |A x − b|, from the normal equations Aᴴ A x = Aᴴ b when there are more rows than
    /// columns, or nil when the columns are dependent.
    func solve(_ b: [Complex.Number<Double>]) -> [Complex.Number<Double>]? {
        precondition(b.count == rows, "one right-hand side per row")
        guard rows > columns else { return Self.eliminate(entries, b, order: columns) }
        let n = columns
        var normal = [Complex.Number<Double>](repeating: .zero, count: n * n)
        var right = [Complex.Number<Double>](repeating: .zero, count: n)
        for row in 0..<rows {
            for i in 0..<n {
                let conjugate = entries[row * n + i].conjugate
                right[i] += conjugate * b[row]
                for j in i..<n { normal[i * n + j] += conjugate * entries[row * n + j] }
            }
        }
        for i in 0..<n {
            for j in 0..<i { normal[i * n + j] = normal[j * n + i].conjugate }
        }
        return Self.eliminate(normal, right, order: n)
    }

    /// The x with A x = b for a square A, by Gaussian elimination with partial pivoting, or nil when a pivot
    /// vanishes.
    static func eliminate(
        _ matrix: [Complex.Number<Double>],
        _ b: [Complex.Number<Double>],
        order n: Int
    ) -> [Complex.Number<Double>]? {
        var a = matrix
        var x = b
        for column in 0..<n {
            var pivot = column
            var largest = 0.0
            for row in column..<n {
                let size = a[row * n + column].magnitude.squared.value
                if size > largest {
                    largest = size
                    pivot = row
                }
            }
            guard largest > 0, largest.isFinite else { return nil }
            if pivot != column {
                for j in 0..<n { a.swapAt(pivot * n + j, column * n + j) }
                x.swapAt(pivot, column)
            }
            let inverse = Complex.Number<Double>.one / a[column * n + column]
            for row in (column + 1)..<n {
                let factor = a[row * n + column] * inverse
                if factor == .zero { continue }
                for j in (column + 1)..<n { a[row * n + j] -= factor * a[column * n + j] }
                x[row] -= factor * x[column]
            }
        }
        for row in stride(from: n - 1, through: 0, by: -1) {
            var sum = x[row]
            for j in (row + 1)..<n { sum -= a[row * n + j] * x[j] }
            x[row] = sum / a[row * n + row]
        }
        return x
    }
}
