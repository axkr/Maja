package rocks.palaiologos.maja.structure;

/**
 * A ring is a set R equipped with two binary operations + and ·, where (R, +) is an abelian group and (R, ·) is a monoid.
 * The multiplication distributes over addition (left and right distributivity).
 */
public interface Ring<T> extends AdditiveAbelianGroup<T>, MultiplicativeMonoid<T> {
    static <T> Ring<T> of(DivisionRing<T> divisionRing) {
        return new Ring<T>() {
            @Override
            public T addInv(T a) {
                return divisionRing.addInv(a);
            }

            @Override
            public T zero() {
                return divisionRing.zero();
            }

            @Override
            public T plus(T a, T b) {
                return divisionRing.plus(a, b);
            }

            @Override
            public T one() {
                return divisionRing.one();
            }

            @Override
            public T dot(T a, T b) {
                return divisionRing.dot(a, b);
            }
        };
    }

    static <T> Ring<T> of(CommutativeRing<T> commutativeRing) {
        return new Ring<T>() {
            @Override
            public T addInv(T a) {
                return commutativeRing.addInv(a);
            }

            @Override
            public T zero() {
                return commutativeRing.zero();
            }

            @Override
            public T plus(T a, T b) {
                return commutativeRing.plus(a, b);
            }

            @Override
            public T one() {
                return commutativeRing.one();
            }

            @Override
            public T dot(T a, T b) {
                return commutativeRing.dot(a, b);
            }
        };
    }

    static <T> Ring<T> of(AdditiveAbelianGroup<T> additiveAbelianGroup, MultiplicativeMonoid<T> multiplicativeMonoid) {
        return new Ring<T>() {
            @Override
            public T addInv(T a) {
                return additiveAbelianGroup.addInv(a);
            }

            @Override
            public T zero() {
                return additiveAbelianGroup.zero();
            }

            @Override
            public T plus(T a, T b) {
                return additiveAbelianGroup.plus(a, b);
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
