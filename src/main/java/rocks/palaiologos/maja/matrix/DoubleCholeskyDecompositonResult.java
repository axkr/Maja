package rocks.palaiologos.maja.matrix;

/**
 * The result of Cholesky decomposition of a matrix.
 *
 * @param l   The lower triangular matrix.
 * @param spd Whether the matrix is symmetric positive definite.
 * @author Palaiologos
 */

import java.util.Objects;

public final class DoubleCholeskyDecompositonResult {
    private final DoubleMatrix l;
    private final boolean spd;

    DoubleCholeskyDecompositonResult(DoubleMatrix l, boolean spd) {
        this.l = l;
        this.spd = spd;
    }

    public DoubleMatrix l() {
        return l;
    }

    public boolean spd() {
        return spd;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == this) return true;
        if (obj == null || obj.getClass() != this.getClass()) return false;
        DoubleCholeskyDecompositonResult that = (DoubleCholeskyDecompositonResult) obj;
        return Objects.equals(this.l, that.l) &&
                this.spd == that.spd;
    }

    @Override
    public int hashCode() {
        return Objects.hash(l, spd);
    }

    @Override
    public String toString() {
        return "DoubleCholeskyDecompositonResult[" +
                "l=" + l + ", " +
                "spd=" + spd + ']';
    }

    /**
     * Solve A * X = B.
     *
     * @return X s.t. L * L' * X = B
     * @throws IllegalArgumentException If matrix row dimensions do not agree or the matrix is not symmetric positive definite.
     */
    public DoubleMatrix solve(DoubleMatrix B) {
        if (B.height() != l.height()) {
            throw new IllegalArgumentException("Matrix row dimensions must be the same.");
        }

        if (!spd) {
            throw new IllegalArgumentException("Matrix is not symmetric positive definite.");
        }

        int n = l.height();
        int nx = B.width();
        DoubleMatrix X = B.copy();

        for (int k = 0; k < n; k++) {
            for (int j = 0; j < nx; j++) {
                for (int i = 0; i < k; i++)
                    X.set(k, j, X.get(k, j) - X.get(i, j) * l().get(k, i));
                X.set(k, j, X.get(k, j) / l().get(k, k));
            }
        }

        for (int k = n - 1; k >= 0; k--) {
            for (int j = 0; j < nx; j++) {
                for (int i = k + 1; i < n; i++)
                    X.set(k, j, X.get(k, j) - X.get(i, j) * l().get(i, k));
                X.set(k, j, X.get(k, j) / l().get(k, k));
            }
        }

        return X;
    }
}
