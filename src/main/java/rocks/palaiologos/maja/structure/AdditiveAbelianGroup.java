package rocks.palaiologos.maja.structure;

/**
 * Group with a commutative binary operation.
 */
public interface AdditiveAbelianGroup<T> extends AdditiveGroup<T> {
    static <T> AdditiveAbelianGroup<T> of(MultiplicativeAbelianGroup<T> group) {
        return new AdditiveAbelianGroup<T>() {
            @Override
            public T zero() {
                return group.one();
            }

            @Override
            public T plus(T a, T b) {
                return group.dot(a, b);
            }

            @Override
            public T addInv(T a) {
                return group.mulInv(a);
            }
        };
    }
}
