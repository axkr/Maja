package rocks.palaiologos.maja.structure;

/**
 * Basic algebraic structure with a closed binary operation. Formally speaking for all a, b in T, dot(a, b) is present in T.
 * No other requirements are imposed.
 */
public interface MultiplicativeGroupoid<T> {
    T dot(T a, T b);

    static <T> MultiplicativeGroupoid<T> of(AdditiveGroupoid<T> groupoid) {
        return groupoid::plus;
    }
}
