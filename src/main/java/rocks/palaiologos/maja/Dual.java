package rocks.palaiologos.maja;

/**
 * A record representing a dual number a + bε.
 *
 * @param a real part
 * @param b dual part
 */

public class Dual {
  public static final Dual NaN = new Dual(Double.NaN);
  public static final Dual INFINITY = new Dual(Double.POSITIVE_INFINITY);
  public static final Dual ZERO = new Dual(0);
  public static final Dual ONE = new Dual(1);
  public static final Dual EPS = new Dual(0, 1);

  private final double a;
  private final double b;

  public Dual(double a, double b) {
    this.a = a;
    this.b = b;
  }

  public double a() {
    return a;
  }

  public double b() {
    return b;
  }

  @Override
  public int hashCode() {
    int result = Double.hashCode(a);
    result = 31 * result + Double.hashCode(b);
    return result;
  }


  /**
   * Construct a complex value with the given real part and zero dual parts.
   */
  public Dual(double re) {
    this(re, 0);
  }

  /**
   * Construct a complex value with real and dual parts.
   */
  public Dual() {
    this(0, 0);
  }

  public static boolean isNaN(Dual c) {
    return Double.isNaN(c.a) || Double.isNaN(c.b);
  }

  public static boolean isInfinite(Dual c) {
    return Double.isInfinite(c.a) || Double.isInfinite(c.b);
  }

  @Override
  public boolean equals(Object o) {
    if (this == o)
      return true;
    if (o == null || getClass() != o.getClass())
      return false;
    Dual dual = (Dual) o;
    return Double.compare(dual.a, a) == 0 && Double.compare(dual.b, b) == 0;
  }

  @Override
  public String toString() {
    if (b < 0) {
      return a + " - " + -b + "ε";
    } else {
      return a + " + " + b + "ε";
    }
  }
}
