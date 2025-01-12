package rocks.palaiologos.maja;

/**
 * A two element tuple class.
 *
 * @param first
 * @param second
 * @param <X>
 * @param <Y>
 */

public class Pair<X, Y> {
  private final X first;
  private final Y second;

  public Pair(X first, Y second) {
    this.first = first;
    this.second = second;
  }

  public X first() {
    return first;
  }

  public Y second() {
    return second;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o)
      return true;
    if (o == null || getClass() != o.getClass())
      return false;

    Pair<?, ?> pair = (Pair<?, ?>) o;

    if (!first.equals(pair.first))
      return false;
    return second.equals(pair.second);
  }

  @Override
  public int hashCode() {
    int result = first.hashCode();
    result = 31 * result + second.hashCode();
    return result;
  }

  @Override
  public String toString() {
    return "Pair{" + "first=" + first + ", second=" + second + '}';
  }
}