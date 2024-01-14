package rocks.palaiologos.maja.structure;

/**
 * An algebraic structure with an associative and commutative binary operation and an identity element.
 * Formally speaking for all a in T, dot(a, one()) = dot(one(), a) = a must hold.
 */
public interface MultiplicativeCommutativeMonoid<T> extends MultiplicativeSemigroup<T> {
    T one();

    static <T> MultiplicativeCommutativeMonoid<T> of(AdditiveCommutativeMonoid<T> semigroup) {
        return new MultiplicativeCommutativeMonoid<T>() {
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
