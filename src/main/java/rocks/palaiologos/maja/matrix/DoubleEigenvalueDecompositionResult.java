package rocks.palaiologos.maja.matrix;

import org.hipparchus.complex.Complex;

public class DoubleEigenvalueDecompositionResult {
  private final DoubleMatrix V;
  private final DoubleMatrix D;
  private final Complex[] e;

  public DoubleEigenvalueDecompositionResult(DoubleMatrix V, DoubleMatrix D, Complex[] e) {
    this.V = V;
    this.D = D;
    this.e = e;
  }

  public DoubleMatrix v() {
    return V;
  }

  public DoubleMatrix d() {
    return D;
  }

  public Complex[] e() {
    return e;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o)
      return true;
    if (o == null || getClass() != o.getClass())
      return false;
    DoubleEigenvalueDecompositionResult that = (DoubleEigenvalueDecompositionResult) o;
    return java.util.Objects.equals(V, that.V) && java.util.Objects.equals(D, that.D)
        && java.util.Arrays.equals(e, that.e);
  }

  @Override
  public int hashCode() {
    int result = java.util.Objects.hash(V, D);
    result = 31 * result + java.util.Arrays.hashCode(e);
    return result;
  }

  @Override
  public String toString() {
      return "DoubleEigenvalueDecompositionResult{" +
              "V=" + V +
              ", D=" + D +
              ", e=" + java.util.Arrays.toString(e) +
              '}';
  }
}
