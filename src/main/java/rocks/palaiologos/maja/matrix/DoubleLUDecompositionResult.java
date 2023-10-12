package rocks.palaiologos.maja.matrix;

import java.util.Objects;

public final class DoubleLUDecompositionResult {
    private final DoubleMatrix lower;
    private final DoubleMatrix upper;
    private final boolean singular;

    /**
     * The result of LU decomposition.
     *
     * @param lower    The lower triangular matrix.
     * @param upper    The upper triangular matrix.
     * @param singular Whether the matrix is singular.
     *                 If the matrix is singular, the upper and lower triangular matrices will be null.
     * @author Palaiologos
     */
    DoubleLUDecompositionResult(DoubleMatrix lower, DoubleMatrix upper, boolean singular) {
        this.lower = lower;
        this.upper = upper;
        this.singular = singular;
    }

    public DoubleMatrix lower() {
        return lower;
    }

    public DoubleMatrix upper() {
        return upper;
    }

    public boolean singular() {
        return singular;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == this) return true;
        if (obj == null || obj.getClass() != this.getClass()) return false;
        DoubleLUDecompositionResult that = (DoubleLUDecompositionResult) obj;
        return Objects.equals(this.lower, that.lower) &&
                Objects.equals(this.upper, that.upper) &&
                this.singular == that.singular;
    }

    @Override
    public int hashCode() {
        return Objects.hash(lower, upper, singular);
    }

    @Override
    public String toString() {
        return "DoubleLUDecompositionResult[" +
                "lower=" + lower + ", " +
                "upper=" + upper + ", " +
                "singular=" + singular + ']';
    }

}
