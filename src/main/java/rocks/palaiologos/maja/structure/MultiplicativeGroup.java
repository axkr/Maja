package rocks.palaiologos.maja.structure;

/**
 * An extension to the concept of monoid. Every element of T must have an inverse element
 * specified by the mulInv() method which is also in T. The inverse element is required to be unique.
 */
public interface MultiplicativeGroup<T> extends MultiplicativeMonoid<T> {
    T mulInv(T a);

    static <T> MultiplicativeGroup<T> of(AdditiveGroup<T> group) {
        return new MultiplicativeGroup<T>() {
            @Override
            public T one() {
                return group.zero();
            }

            @Override
            public T dot(T a, T b) {
                return group.plus(a, b);
            }

            @Override
            public T mulInv(T a) {
                return group.addInv(a);
            }
        };
    }
}
