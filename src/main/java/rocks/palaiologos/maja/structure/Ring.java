package rocks.palaiologos.maja.structure;

/**
 * A ring is a set R equipped with two binary operations + and ·, where (R, +) is an abelian group and (R, ·) is a monoid.
 * The multiplication distributes over addition (left and right distributivity).
 */
public interface Ring<T> extends AdditiveAbelianGroup<T>, MultiplicativeMonoid<T> {
}
