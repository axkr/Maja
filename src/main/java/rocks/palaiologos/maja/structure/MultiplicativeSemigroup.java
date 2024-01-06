package rocks.palaiologos.maja.structure;

/**
 * Basic algebraic structure with an associative binary operation.
 * Formally speaking, for all a, b, c in T, dot(dot(a, b), c) = dot(a, dot(b, c)) must hold.
 * Further, the operation is closed.
 */
public interface MultiplicativeSemigroup<T> extends MultiplicativeGroupoid<T> {
    T dot(T a, T b);

    static <T> MultiplicativeSemigroup<T> of(AdditiveSemigroup<T> semigroup) {
        return semigroup::plus;
    }
}
