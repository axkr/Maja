package rocks.palaiologos.maja.structure;

/**
 * A ring is a set R equipped with two binary operations + and ·, where (R, +) is a commutative monoid and (R, ·) is a monoid.
 * The multiplication distributes over addition (left and right distributivity).
 */
public interface SemiRing<T> extends AdditiveCommutativeMonoid<T>, MultiplicativeMonoid<T> {
    static <T> SemiRing<T> of(Ring<T> ring) {
        return new SemiRing<T>() {
            @Override
            public T zero() {
                return ring.zero();
            }

            @Override
            public T plus(T a, T b) {
                return ring.plus(a, b);
            }

            @Override
            public T one() {
                return ring.one();
            }

            @Override
            public T dot(T a, T b) {
                return ring.dot(a, b);
            }
        };
    }

    static <T> SemiRing<T> of(AdditiveCommutativeMonoid<T> additiveCommutativeMonoid, MultiplicativeMonoid<T> multiplicativeMonoid) {
        return new SemiRing<T>() {
            @Override
            public T zero() {
                return additiveCommutativeMonoid.zero();
            }

            @Override
            public T plus(T a, T b) {
                return additiveCommutativeMonoid.plus(a, b);
            }

            @Override
            public T one() {
                return multiplicativeMonoid.one();
            }

            @Override
            public T dot(T a, T b) {
                return multiplicativeMonoid.dot(a, b);
            }
        };
    }
}
