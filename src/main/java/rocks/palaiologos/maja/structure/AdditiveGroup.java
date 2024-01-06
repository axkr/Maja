package rocks.palaiologos.maja.structure;

/**
 * An extension to the concept of monoid. Every element of T must have an inverse element
 * specified by the addInv() method which is also in T. The inverse element is required to be unique.
 */
public interface AdditiveGroup<T> extends AdditiveMonoid<T> {
    T addInv(T a);

    static <T> AdditiveGroup<T> of(MultiplicativeGroup<T> group) {
        return new AdditiveGroup<T>() {
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
