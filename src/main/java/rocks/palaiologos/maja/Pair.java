package rocks.palaiologos.maja;

/**
 * A two element tuple class.
 *
 * @param first
 * @param second
 * @param <X>
 * @param <Y>
 */

import com.github.bsideup.jabel.Desugar;

@Desugar
public record Pair<X, Y>(X first, Y second) {
}
