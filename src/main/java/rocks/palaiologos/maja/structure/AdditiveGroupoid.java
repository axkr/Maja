package rocks.palaiologos.maja.structure;

/**
 * Basic algebraic structure with a closed binary operation. Formally speaking for all a, b in T, plus(a, b) is present in T.
 * No other requirements are imposed.
 */
public interface AdditiveGroupoid<T> {
    T plus(T a, T b);

    static <T> AdditiveGroupoid<T> of(MultiplicativeGroupoid<T> groupoid) {
        return groupoid::dot;
    }
}
