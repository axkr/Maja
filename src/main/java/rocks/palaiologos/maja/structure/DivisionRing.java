package rocks.palaiologos.maja.structure;

/**
 * A division ring is a set R equipped with two binary operations + and ·, where (R, +) is an abelian group and (R, ·) is a group.
 * The multiplication distributes over addition (left and right distributivity).
 */
public interface DivisionRing<T> extends AdditiveAbelianGroup<T>, MultiplicativeGroup<T> {
}
