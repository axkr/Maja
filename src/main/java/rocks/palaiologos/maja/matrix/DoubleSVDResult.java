package rocks.palaiologos.maja.matrix;

/**
 * The result of SV decomposition.
 *
 * @param rank The effective numerical rank.
 * @param norm L_2 norm of the matrix.
 * @param conditionNumber Ratio of largest to smallest singular value.
 * @param inverseConditionNumber Ratio of smallest to largest singular value.
 * @param singularValues The singular values.
 * @param U The left singular vectors.
 * @param V The right singular vectors.
 * @author Palaiologos
 */

public class DoubleSVDResult {
  private final int rank;
  private final double norm;
  private final double conditionNumber;
  private final double inverseConditionNumber;
  private final double[] singularValues;
  private final DoubleMatrix U;
  private final DoubleMatrix V;

  public DoubleSVDResult(int rank, double norm, double conditionNumber,
      double inverseConditionNumber, double[] singularValues, DoubleMatrix U, DoubleMatrix V) {
    this.rank = rank;
    this.norm = norm;
    this.conditionNumber = conditionNumber;
    this.inverseConditionNumber = inverseConditionNumber;
    this.singularValues = singularValues;
    this.U = U;
    this.V = V;
  }

  public int rank() {
    return rank;
  }

  public double norm() {
    return norm;
  }

  public double conditionNumber() {
    return conditionNumber;
  }

  public double inverseConditionNumber() {
    return inverseConditionNumber;
  }

  public double[] singularValues() {
    return singularValues;
  }

  public DoubleMatrix u() {
    return U;
  }

  public DoubleMatrix v() {
    return V;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o)
      return true;
    if (o == null || getClass() != o.getClass())
      return false;
    DoubleSVDResult that = (DoubleSVDResult) o;
    return rank == that.rank && Double.compare(that.norm, norm) == 0
        && Double.compare(that.conditionNumber, conditionNumber) == 0
        && Double.compare(that.inverseConditionNumber, inverseConditionNumber) == 0
        && java.util.Arrays.equals(singularValues, that.singularValues)
        && java.util.Objects.equals(U, that.U) && java.util.Objects.equals(V, that.V);
  }

  @Override
  public int hashCode() {
    int result = java.util.Objects.hash(rank, norm, conditionNumber, inverseConditionNumber, U, V);
    result = 31 * result + java.util.Arrays.hashCode(singularValues);
    return result;
  }

  @Override
  public String toString() {
    return "DoubleSVDResult{" + "rank=" + rank + ", norm=" + norm + ", conditionNumber="
        + conditionNumber + ", inverseConditionNumber=" + inverseConditionNumber
        + ", singularValues=" + java.util.Arrays.toString(singularValues) + ", U=" + U + ", V=" + V
        + '}';
  }
}
