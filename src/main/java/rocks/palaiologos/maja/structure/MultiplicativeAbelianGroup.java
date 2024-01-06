package rocks.palaiologos.maja.structure;

/**
 * Group with a commutative binary operation.
 */
public interface MultiplicativeAbelianGroup<T> extends MultiplicativeGroup<T> {
    static <T> MultiplicativeAbelianGroup<T> of(AdditiveAbelianGroup<T> group) {
        return new MultiplicativeAbelianGroup<>() {
            @Override
            public T dot(T a, T b) {
                return group.plus(a, b);
            }

            @Override
            public T one() {
                return group.zero();
            }

            @Override
            public T mulInv(T a) {
                return group.addInv(a);
            }
        };
    }
}
