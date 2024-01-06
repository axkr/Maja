package rocks.palaiologos.maja.structure;

/**
 * An algebraic structure with an associative binary operation (implied by semigroup properties) and an identity element.
 * Formally speaking for all a in T, dot(a, one()) = dot(one(), a) = a must hold.
 */
public interface MultiplicativeMonoid<T> extends MultiplicativeSemigroup<T> {
    T one();

    static <T> MultiplicativeMonoid<T> of(AdditiveMonoid<T> semigroup) {
        return new MultiplicativeMonoid<T>() {
            @Override
            public T one() {
                return semigroup.zero();
            }

            @Override
            public T dot(T a, T b) {
                return semigroup.plus(a, b);
            }
        };
    }
}
