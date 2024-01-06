package rocks.palaiologos.maja.structure;

/**
 * An algebraic structure with an associative binary operation (implied by semigroup properties) and an identity element.
 * Formally speaking for all a in T, add(a, zero()) = add(zero(), a) = a must hold.
 */
public interface AdditiveMonoid<T> extends AdditiveSemigroup<T> {
    T zero();

    static <T> AdditiveMonoid<T> of(MultiplicativeMonoid<T> monoid) {
        return new AdditiveMonoid<T>() {
            @Override
            public T zero() {
                return monoid.one();
            }

            @Override
            public T plus(T a, T b) {
                return monoid.dot(a, b);
            }
        };
    }
}
