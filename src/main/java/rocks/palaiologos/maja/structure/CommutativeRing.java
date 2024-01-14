package rocks.palaiologos.maja.structure;

/**
 * A division ring is a set R equipped with two binary operations + and ·, where (R, +) and (R, ·) are abelian groups.
 * The multiplication distributes over addition (left and right distributivity).
 */
public interface CommutativeRing<T> extends AdditiveAbelianGroup<T>, MultiplicativeAbelianGroup<T> {
}
