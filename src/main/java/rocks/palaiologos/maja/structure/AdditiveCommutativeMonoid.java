package rocks.palaiologos.maja.structure;

/**
 * An algebraic structure with an associative and commutative binary operation and an identity element.
 * Formally speaking for all a in T, add(a, zero()) = add(zero(), a) = a must hold.
 */
public interface AdditiveCommutativeMonoid<T> extends AdditiveSemigroup<T> {
    T zero();

    static <T> AdditiveCommutativeMonoid<T> of(MultiplicativeCommutativeMonoid<T> monoid) {
        return new AdditiveCommutativeMonoid<T>() {
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
