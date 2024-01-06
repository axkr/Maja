package rocks.palaiologos.maja.structure;

/**
 * Basic algebraic structure with an associative binary operation.
 * Formally speaking, for all a, b, c in T, plus(plus(a, b), c) = plus(a, plus(b, c)) must hold.
 * Further, the operation is closed.
 */
public interface AdditiveSemigroup<T> extends AdditiveGroupoid<T> {
    T plus(T a, T b);

    static <T> AdditiveSemigroup<T> of(MultiplicativeSemigroup<T> semigroup) {
        return semigroup::dot;
    }
}
