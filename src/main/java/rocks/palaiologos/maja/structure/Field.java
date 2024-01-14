package rocks.palaiologos.maja.structure;

/**
 * A a commutative division ring (i.e. a commutative ring which contains a multiplicative inverse for every nonzero element).
 */
public interface Field<T> extends AdditiveAbelianGroup<T>, MultiplicativeAbelianGroup<T> {
}
