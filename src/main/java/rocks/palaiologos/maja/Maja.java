package rocks.palaiologos.maja;

import java.math.BigDecimal;
import java.util.Arrays;
import java.util.Random;
import java.util.function.BiFunction;
import java.util.function.Function;
import org.hipparchus.complex.Complex;
import rocks.palaiologos.maja.structure.AdditiveGroup;
import rocks.palaiologos.maja.structure.AdditiveGroupoid;
import rocks.palaiologos.maja.structure.MultiplicativeGroup;
import rocks.palaiologos.maja.structure.MultiplicativeGroupoid;

/**
 * A slick numerics-oriented Mathematical library for Java.
 *
 * @author Palaiologos
 */
public class Maja {
    /**
     * SI prefix for 10^1.
     */
    public static final double DECA = 10.0e1;
    /**
     * SI prefix for 10^2.
     */
    public static final double HECTO = 10.0e2;
    /**
     * SI prefix for 10^3.
     */
    public static final double KILO = 10.0e3;
    /**
     * SI prefix for 10^6.
     */
    public static final double MEGA = 10.0e6;
    /**
     * SI prefix for 10^9.
     */
    public static final double GIGA = 10.0e9;
    /**
     * SI prefix for 10^12.
     */
    public static final double TERA = 10.0e12;
    /**
     * SI prefix for 10^15.
     */
    public static final double PETA = 10.0e15;
    /**
     * SI prefix for 10^18.
     */
    public static final double EXA = 10.0e18;
    /**
     * SI prefix for 10^21.
     */
    public static final double ZETTA = 10.0e21;
    /**
     * SI prefix for 10^24.
     */
    public static final double YOTTA = 10.0e24;
    /**
     * SI prefix for 10^-1.
     */
    public static final double DECI = 10.0e-1;
    /**
     * SI prefix for 10^-2.
     */
    public static final double CENTI = 10.0e-2;
    /**
     * SI prefix for 10^-3.
     */
    public static final double MILLI = 10.0e-3;
    /**
     * SI prefix for 10^-6.
     */
    public static final double MICRO = 10.0e-6;
    /**
     * SI prefix for 10^-9.
     */
    public static final double NANO = 10.0e-9;
    /**
     * SI prefix for 10^-12.
     */
    public static final double PICO = 10.0e-12;
    /**
     * SI prefix for 10^-15.
     */
    public static final double FEMTO = 10.0e-15;
    /**
     * SI prefix for 10^-18.
     */
    public static final double ATTO = 10.0e-18;
    /**
     * SI prefix for 10^-21.
     */
    public static final double ZEPTO = 10.0e-21;
    /**
     * SI prefix for 10^-24.
     */
    public static final double YOCTO = 10.0e-24;

    /**
     * The value of ln(2).
     */
    public static final double LN2 = 0.6931471805599453094172321214581765680755001343602552541206800094;
    /**
     * The value of ln(10).
     */
    public static final double LN10 = 2.3025850929940456840179914546843642076011014886287729760333279009;
    /**
     * The value of log2(e).
     */
    public static final double LOG2E = 1.4426950408889634073599246810018921374266459541529859341354494069;
    /**
     * The value of pi.
     */
    public static final double PI = 3.1415926535897932384626433832795028841971693993751058209749445923;
    /**
     * The value of pi/2.
     */
    public static final double PI_2 = 1.5707963267948966192313216916397514420985846996875529104874722961;
    /**
     * The value of pi/4.
     */
    public static final double PI_4 = 0.7853981633974483096156608458198757210492923498437764552437361480;
    /**
     * The value of 2*PI.
     */
    public static final double TWO_PI = 6.2831853071795864769252867665590057683943387987502116419498891846;
    /**
     * The value of 1/pi.
     */
    public static final double ONE_OVER_PI = 0.3183098861837906715377675267450287240689192914809128974953346881;
    /**
     * The value of e.
     */
    public static final double E = 2.7182818284590452353602874713526624977572470936999595749669676277;
    /**
     * The value of 1/e.
     */
    public static final double ONE_OVER_E = 0.3678794411714423215955237701614608674458111310317678345078368016;
    /**
     * The epsilon value for double precision numbers.
     */
    public static final double EPSILON = Math.ulp(1.0d);
    /**
     * The Euler-Mascheroni constant.
     */
    public static final double EULER_GAMMA = 0.5772156649015328606065120900824024310421593359399235988057672348;
    /**
     * The golden ratio.
     */
    public static final double GOLDEN_RATIO = 1.6180339887498948482045868343656381177203091798057628621354486227;
    /**
     * Apery's constant - zeta(3).
     */
    public static final double APERY_CONSTANT = 1.2020569031595942853997381615114499907649862923404988817922715553;
    /**
     * Glaisher-Kinkelin constant.
     */
    public static final double GLAISHER_CONSTANT = 1.2824271291006226368753425688697917277676889273250011920637400217;
    /**
     * Catalan's constant.
     */
    public static final double CATALAN_CONSTANT = 0.9159655941772190150546035149323841107741493742816721342664981196;
    /**
     * Golomb-Dickman constant.
     */
    public static final double GOLOMB_DICKMAN_CONSTANT = 0.6243299885435508709929363831008372441796426201805292869735519024;
    /**
     * Mills' constant.
     */
    public static final double MILLS_CONSTANT = 1.3063778838630806904686144926026057129167845851567136443680537599;
    /**
     * Feigenbaum constant.
     */
    public static final double FEIGENBAUM_CONSTANT = 4.6692016091029906718532038204662016172581855774757686327456513430;
    /**
     * Khinchin's constant.
     */
    public static final double KHINCHIN_CONSTANT = 2.6854520010653064453097148354817956938203822939944629530511523455;
    /**
     * The imaginary unit.
     */
    public static final Complex I = new Complex(0.0, 1.0);
    /**
     * The random number generator used by this class.
     */
    private static final Random random = new Random();

    private Maja() {
    }

    /**
     * Adds two double precision numbers together.
     *
     * @param x
     * @param y
     * @return x + y
     */
    public static double add(double x, double y) {
        return x + y;
    }

    /**
     * Adds two values through an additive groupoid.
     */
    public static <T> T add(AdditiveGroupoid<T> groupoid, T x, T y) {
        return groupoid.plus(x, y);
    }

    /**
     * Subtracts two double precision numbers.
     *
     * @param x
     * @param y
     * @return x - y
     */
    public static double sub(double x, double y) {
        return x - y;
    }

    /**
     * Subtracts two values through an additive group.
     */
    public static <T> T add(AdditiveGroup<T> group, T x, T y) {
        return group.plus(x, group.addInv(y));
    }

    /**
     * Multiplies two double precision numbers.
     *
     * @param x
     * @param y
     * @return x * y
     */
    public static double mul(double x, double y) {
        return x * y;
    }

    /**
     * Multiplies two values through a multiplicative groupoid.
     */
    public static <T> T mul(MultiplicativeGroupoid<T> groupoid, T x, T y) {
        return groupoid.dot(x, y);
    }

    /**
     * Divides two double precision numbers.
     *
     * @param x
     * @param y
     * @return x / y
     */
    public static double div(double x, double y) {
        return x / y;
    }

    /**
     * Divides two values through a multiplicative group.
     */
    public static <T> T mul(MultiplicativeGroup<T> group, T x, T y) {
        return group.dot(x, group.mulInv(y));
    }

    /**
     * Returns the modulus of two double precision numbers.
     *
     * @param x
     * @param y
     * @return x mod y
     */
    public static double mod(double x, double y) {
        final double r = x % y;
        if (r != 0 && Math.signum(x) != Math.signum(y))
            return r + y;
        return r;
    }

    /**
     * Returns the modulus of two numbers.
     *
     * @param x
     * @param y
     * @return x mod y
     */
    public static long mod(long x, long y) {
        final long r = x % y;
        if (r != 0 && Math.signum(x) != Math.signum(y))
            return r + y;
        return r;
    }

    /**
     * Returns the remainder that results from dividing two double precision numbers.
     *
     * @param x
     * @param y
     * @return x % y
     */
    public static double rem(double x, double y) {
        return x % y;
    }

    /**
     * Returns the absolute value of a double precision number.
     *
     * @param x
     * @return |x|
     * @see java.lang.Math#abs(double)
     */
    public static double abs(double x) {
        return Math.abs(x);
    }

    /**
     * Returns the absolute value of a number.
     *
     * @param x
     * @return |x|
     */
    public static long abs(long x) {
        return Math.abs(x);
    }

    /**
     * Returns the value of the inverse cosine of a double precision number.
     *
     * @param x
     * @return cos^-1(x)
     * @see java.lang.Math#acos(double)
     */
    public static double acos(double x) {
        return Math.acos(x);
    }

    /**
     * Returns the value of the inverse sine of a double precision number.
     *
     * @param x
     * @return sin^-1(x)
     * @see java.lang.Math#asin(double)
     */
    public static double asin(double x) {
        return Math.asin(x);
    }

    /**
     * Returns the value of the inverse tangent of a double precision number.
     *
     * @param x
     * @return tan^-1(x)
     * @see java.lang.Math#atan(double)
     */
    public static double atan(double x) {
        return Math.atan(x);
    }

    /**
     * Returns the value of the inverse tangent of two double precision numbers.
     *
     * @param y
     * @param x
     * @return atan2(y, x)
     * @see java.lang.Math#atan2(double, double)
     */
    public static double atan2(double y, double x) {
        return Math.atan2(y, x);
    }

    /**
     * Returns the value of the cube root of a double precision number.
     *
     * @param x
     * @return x^(1/3)
     * @see java.lang.Math#cbrt(double)
     */
    public static double cbrt(double x) {
        return Math.cbrt(x);
    }

    /**
     * Returns the value of the smallest integer greater than or equal to a double precision number.
     *
     * @param x
     * @return ceil(x)
     * @see java.lang.Math#ceil(double)
     */
    public static double ceil(double x) {
        return Math.ceil(x);
    }

    /**
     * Returns the value of the smallest integer greater than or equal to a double precision number.
     *
     * @param x
     * @return ceil(x)
     * @see java.lang.Math#ceil(double)
     */
    public static double cos(double x) {
        return Math.cos(x);
    }

    /**
     * Returns the hyperbolic cosine of a double precision number.
     *
     * @param x
     * @return cosh(x)
     * @see java.lang.Math#cosh(double)
     */
    public static double cosh(double x) {
        return Math.cosh(x);
    }

    /**
     * Returns the value of the exponential function of a double precision number.
     *
     * @param x
     * @return e^x
     * @see java.lang.Math#exp(double)
     */
    public static double exp(double x) {
        return Math.exp(x);
    }

    /**
     * Returns the value of the exponential function of a double precision number minus one.
     *
     * @param x
     * @return e^x - 1
     * @see java.lang.Math#expm1(double)
     */
    public static double expm1(double x) {
        return Math.expm1(x);
    }

    /**
     * Returns the value of the largest integer less than or equal to a double precision number.
     *
     * @param x
     * @return floor(x)
     * @see java.lang.Math#floor(double)
     */
    public static double floor(double x) {
        return Math.floor(x);
    }

    /**
     * Returns the length of the hypotenuse of a right triangle with sides of length x and y.
     *
     * @param x
     * @param y
     * @return sqrt(x ^ 2 + y ^ 2)
     * @see java.lang.Math#hypot(double, double)
     */
    public static double hypot(double x, double y) {
        return Math.hypot(x, y);
    }

    /**
     * Returns the value of the natural logarithm of a double precision number.
     *
     * @param x
     * @return ln(x)
     * @see java.lang.Math#log(double)
     */
    public static double log(double x) {
        return Math.log(x);
    }

    /**
     * Returns the value of the base 10 logarithm of a double precision number.
     *
     * @param x
     * @return log10(x)
     * @see java.lang.Math#log10(double)
     */
    public static double log10(double x) {
        return Math.log10(x);
    }

    /**
     * Returns the value of the natural logarithm of a double precision number plus one.
     *
     * @param x
     * @return ln(x + 1)
     * @see java.lang.Math#log1p(double)
     */
    public static double log1p(double x) {
        return Math.log1p(x);
    }

    /**
     * Returns the value of the base 2 logarithm of a double precision number.
     *
     * @param x
     * @return log2(x)
     */
    public static double log2(double x) {
        return Math.log(x) / LN2;
    }

    /**
     * Returns the larger of two numbers (maximum).
     *
     * @param x
     * @param y
     * @return max(x, y)
     * @see java.lang.Math#max(double, double)
     */
    public static double max(double x, double y) {
        return Math.max(x, y);
    }

    /**
     * Returns the smaller of two numbers (minimum).
     *
     * @param x
     * @param y
     * @return min(x, y)
     * @see java.lang.Math#min(double, double)
     */
    public static double min(double x, double y) {
        return Math.min(x, y);
    }

    /**
     * Returns the value of the first argument raised to the power of the second argument.
     *
     * @param x
     * @param y
     * @return x^y
     * @see java.lang.Math#pow(double, double)
     */
    public static double pow(double x, double y) {
        return Math.pow(x, y);
    }

    /**
     * Returns the sign of a double precision number.
     *
     * @param x
     * @return -1 if x &lt; 0, 0 if x == 0, 1 if x &gt; 0
     * @see java.lang.Math#signum(double)
     */
    public static double signum(double x) {
        return Math.signum(x);
    }

    /**
     * Returns the sign of a number.
     *
     * @param x
     * @return -1 if x &lt; 0, 1 if x &gt; 0
     */
    public static long signum(long x) {
        if (x < 0) return -1;
        if (x > 0) return 1;
        return 0;
    }

    /**
     * Returns the sign of a complex number.
     *
     * @param x
     * @return A complex number where the real and imaginary parts correspond to the sign of the real and imaginary parts of x.
     */
    public static Complex signum(Complex x) {
        return new Complex(signum(x.getReal()), signum(x.getImaginary()));
    }

    /**
     * Returns the value of the sine of a double precision number.
     *
     * @param x
     * @return sin(x)
     * @see java.lang.Math#sin(double)
     */
    public static double sin(double x) {
        return Math.sin(x);
    }

    /**
     * Rounds both the real and imaginary parts of a complex number.
     *
     * @param x
     * @return ceil(x)
     */
    public static Complex ceil(Complex x) {
        return new Complex(Math.ceil(x.getReal()), Math.ceil(x.getImaginary()));
    }

    /**
     * Rounds both the real and imaginary parts of a complex number.
     *
     * @param x
     * @return floor(x)
     */
    public static Complex floor(Complex x) {
        return new Complex(Math.floor(x.getReal()), Math.floor(x.getImaginary()));
    }

    /**
     * Rounds both the real and imaginary parts of a complex number.
     *
     * @param x
     * @return round(x)
     */
    public static Complex round(Complex x) {
        return new Complex(Math.round(x.getReal()), Math.round(x.getImaginary()));
    }

    /**
     * Returns the sinc function of a double precision number, defined
     * as sin(x) / x, except for x = 0, where sinc(x) = 1.
     *
     * @param x
     * @return sinc(x)
     */
    public static double sinc(double x) {
        if (x == 0)
            return 1;
        return Math.sin(x) / x;
    }

    /**
     * Returns the sinc function of a complex number, defined
     * as sin(x) / x, except for x = 0, where sinc(x) = 1.
     *
     * @param x
     * @return sinc(x)
     */
    public static Complex sinc(Complex x) {
        if (eq(x, Complex.ZERO))
            return Complex.ONE;
        return div(sin(x), x);
    }

    /**
     * Returns the hyperbolic sine of a double precision number.
     *
     * @param x
     * @return sinh(x)
     * @see java.lang.Math#sinh(double)
     */
    public static double sinh(double x) {
        return Math.sinh(x);
    }

    /**
     * Returns the square root of a double precision number.
     *
     * @param x
     * @return x^(1/2)
     * @see java.lang.Math#sqrt(double)
     */
    public static double sqrt(double x) {
        return Math.sqrt(x);
    }

    /**
     * Returns the tangent of a double precision number.
     *
     * @param x
     * @return tan(x)
     * @see java.lang.Math#tan(double)
     */
    public static double tan(double x) {
        return Math.tan(x);
    }

    /**
     * Returns the hyperbolic tangent of a double precision number.
     *
     * @param x
     * @return tanh(x)
     * @see java.lang.Math#tanh(double)
     */
    public static double tanh(double x) {
        return Math.tanh(x);
    }

    /**
     * Converts the value to degrees from radians.
     *
     * @param x
     * @return x * 180 / PI
     * @see java.lang.Math#toDegrees(double)
     */
    public static double toDegrees(double x) {
        return Math.toDegrees(x);
    }

    /**
     * Converts the value to radians from degrees.
     *
     * @param x
     * @return x * PI / 180
     * @see java.lang.Math#toRadians(double)
     */
    public static double toRadians(double x) {
        return Math.toRadians(x);
    }

    /**
     * Returns the value of the <i>unit in last place</i> (ULP) of a double precision number.
     *
     * @param x
     * @return ULP(x)
     * @see java.lang.Math#ulp(double)
     */
    public static double ulp(double x) {
        return Math.ulp(x);
    }

    /**
     * Returns the value of the fused multiply-add operation.
     *
     * @param a
     * @param b
     * @param c
     * @return a * b + c
     */
    public static double fma(double a, double b, double c) {
        if (Double.isNaN(a) || Double.isNaN(b) || Double.isNaN(c)) {
            return Double.NaN;
        } else { // All inputs non-NaN
            boolean infiniteA = Double.isInfinite(a);
            boolean infiniteB = Double.isInfinite(b);
            boolean infiniteC = Double.isInfinite(c);
            double result;

            if (infiniteA || infiniteB || infiniteC) {
                if (infiniteA && b == 0.0 ||
                        infiniteB && a == 0.0 ) {
                    return Double.NaN;
                }
                double product = a * b;
                if (Double.isInfinite(product) && !infiniteA && !infiniteB) {
                    // Intermediate overflow; might cause a
                    // spurious NaN if added to infinite c.
                    assert Double.isInfinite(c);
                    return c;
                } else {
                    result = product + c;
                    assert !Double.isFinite(result);
                    return result;
                }
            } else { // All inputs finite
                BigDecimal product = (new BigDecimal(a)).multiply(new BigDecimal(b));
                if (c == 0.0) { // Positive or negative zero
                    // If the product is an exact zero, use a
                    // floating-point expression to compute the sign
                    // of the zero final result. The product is an
                    // exact zero if and only if at least one of a and
                    // b is zero.
                    if (a == 0.0 || b == 0.0) {
                        return a * b + c;
                    } else {
                        // The sign of a zero addend doesn't matter if
                        // the product is nonzero. The sign of a zero
                        // addend is not factored in the result if the
                        // exact product is nonzero but underflows to
                        // zero; see IEEE-754 2008 section 6.3 "The
                        // sign bit".
                        return product.doubleValue();
                    }
                } else {
                    return product.add(new BigDecimal(c)).doubleValue();
                }
            }
        }
    }

    /**
     * Returns the value of the complex fused multiply-add operation.
     *
     * @param x
     * @param y
     * @param z
     * @return x * y + z
     */
    public static Complex fma(Complex x, Complex y, Complex z) {
        return add(mul(x, y), z);
    }

    /**
     * Returns the value of the double precision number adjacent to the first argument in the direction of the second argument.
     *
     * @param x
     * @param y
     * @return nextAfter(x, y)
     * @see java.lang.Math#nextAfter(double, double)
     */
    public static double nextAfter(double x, double y) {
        return Math.nextAfter(x, y);
    }

    /**
     * Returns the value of the double precision number adjacent to the first argument in the direction of positive infinity.
     *
     * @param x
     * @return nextUp(x)
     * @see java.lang.Math#nextUp(double)
     */
    public static double nextUp(double x) {
        return Math.nextUp(x);
    }

    /**
     * Returns the value of the double precision number adjacent to the first argument in the direction of negative infinity.
     *
     * @param x
     * @return nextDown(x)
     * @see java.lang.Math#nextDown(double)
     */
    public static double nextDown(double x) {
        return Math.nextDown(x);
    }

    /**
     * Returns the correctly rounded value of x*2^n.
     *
     * @param x
     * @param n
     * @return scalb(x, n)
     * @see java.lang.Math#scalb(double, int)
     */
    public static double scalb(double x, int n) {
        return Math.scalb(x, n);
    }

    /**
     * Copy the sign of the second argument to the first argument.
     *
     * @param x
     * @param y
     * @return copySign(x, y)
     * @see java.lang.Math#copySign(double, double)
     */
    public static double copySign(double x, double y) {
        return Math.copySign(x, y);
    }

    /**
     * Copy the sign of the second argument to the first argument.
     *
     * @param x
     * @param y
     * @return copySign(x, y)
     */
    public static long copySign(long x, long y) {
        return signum(y) * abs(x);
    }

    /**
     * Copy the sign of the second argument to the first argument.
     *
     * @param x
     * @param y
     * @return copySign(x, y)
     */
    public static Complex copySign(Complex x, Complex y) {
        return mul(signum(y), absparts(x));
    }

    /**
     * Compute the absolute values of parts of a complex number.
     *
     * @param x
     * @return |re(x)| + i|im(x)|
     */
    public static Complex absparts(Complex x) {
        return new Complex(abs(x.getReal()), abs(x.getImaginary()));
    }

    /**
     * Returns the exponent of a double precision number.
     *
     * @param x
     * @return exponent(x)
     * @see java.lang.Math#getExponent(double)
     */
    public static int getExponent(double x) {
        return Math.getExponent(x);
    }

    /**
     * Rounds a double precision number to the nearest integer.
     *
     * @param x
     * @return round(x)
     * @see java.lang.Math#round(double)
     */
    public static long round(double x) {
        return Math.round(x);
    }

    /**
     * Returns a random double precision number in range [0, 1).
     *
     * @return random()
     */
    public static double random() {
        return Math.random();
    }

    /**
     * Returns a random double precision number in range [min, max).
     *
     * @param min
     * @param max
     * @return random() * (max - min) + min
     */
    public static double random(double min, double max) {
        return min + (max - min) * Math.random();
    }

    /**
     * Returns a random long integer in range [min, max).
     *
     * @param min
     * @param max
     * @return random() * (max - min) + min
     */
    public static long random(long min, long max) {
        long r = random.nextLong();
        if (min < max) {
            final long n = max - min;
            final long m = n - 1;
            if ((n & m) == 0L) {
                r = (r & m) + min;
            } else if (n > 0L) {
                for (long u = r >>> 1;
                     u + m - (r = u % n) < 0L;
                     u = random.nextLong() >>> 1)
                    ;
                r += min;
            }
            else {
                while (r < min || r >= max)
                    r = random.nextLong();
            }
        }
        return r;
    }

    /**
     * Returns a random double precision number in range [0, max).
     *
     * @param max
     * @return random() * max
     */
    public static double random(double max) {
        return max * Math.random();
    }

    /**
     * Returns a random long integer in range [0, max).
     *
     * @param max
     * @return random() * max
     */
    public static long random(long max) {
        final long m = max - 1;
        long r = random.nextLong();
        if ((max & m) == 0L) {
            r &= m;
        } else {
            for (long u = r >>> 1;
                 u + m - (r = u % max) < 0L;
                 u = random.nextLong() >>> 1)
                ;
        }
        return r;
    }

    /**
     * Return true if x &lt; y.
     *
     * @param x
     * @param y
     * @return x &lt; y
     */
    public static boolean lt(double x, double y) {
        return x < y;
    }

    /**
     * Return true if x &lt;= y.
     *
     * @param x
     * @param y
     * @return x &lt;= y
     */
    public static boolean le(double x, double y) {
        return x <= y;
    }

    /**
     * Return true if x &gt; y.
     *
     * @param x
     * @param y
     * @return x &gt; y
     */
    public static boolean gt(double x, double y) {
        return x > y;
    }

    /**
     * Return true if x &gt;= y.
     *
     * @param x
     * @param y
     * @return x &gt;= y
     */
    public static boolean ge(double x, double y) {
        return x >= y;
    }

    /**
     * Return true if x == y.
     *
     * @param x
     * @param y
     * @return x == y
     */
    public static boolean eq(double x, double y) {
        return x == y;
    }

    /**
     * Return true if x != y.
     *
     * @param x
     * @param y
     * @return x != y
     */
    public static boolean ne(double x, double y) {
        return x != y;
    }

    /**
     * Compare two double precision numbers.
     *
     * @param x
     * @param y
     * @return 1 if x &gt; y, 0 if x == y, -1 if x &lt; y
     */
    public static int compare(double x, double y) {
        return Double.compare(x, y);
    }

    /**
     * Compare two integers.
     *
     * @param x
     * @param y
     * @return 1 if x &gt; y, 0 if x == y, -1 if x &lt; y
     */
    public static int compare(long x, long y) {
        return Long.compare(x, y);
    }

    /**
     * Return true if x is approximately equal to y.
     *
     * @param x
     * @param y
     * @param tol
     * @return abs(x - y) &lt;= tol
     */
    public static boolean eq(double x, double y, double tol) {
        return Math.abs(x - y) <= tol;
    }

    /**
     * Determine whether a number is a perfect square.
     *
     * @param x
     * @return true if x is a perfect square, false otherwise
     */
    public static boolean isPerfectSquare(long x) {
        return IsSquare.isSquare(x);
    }

    /**
     * Linearly map a value from one range to another. Input range must not be empty.
     *
     * @param inRangeStart  Input range start
     * @param inRangeEnd    Input range end
     * @param outRangeStart Output range start
     * @param outRangeEnd   Output range end
     * @param value         Value to map
     * @return Mapped value. Values outside the input range are not clamped to output range
     */
    public static double linearMap(double inRangeStart, double inRangeEnd, double outRangeStart, double outRangeEnd, double value) {
        return outRangeStart + (value - inRangeStart) * (outRangeEnd - outRangeStart) / (inRangeEnd - inRangeStart);
    }

    /**
     * Linearly normalise value from a range. Range must not be empty.
     *
     * @param rangeStart Range start normalized to 0
     * @param rangeEnd   Range end normalized to 1
     * @param value      Value to normalize
     * @return Normalized value. Values outside the range are not clamped to 0 and 1
     */
    public static double linearNorm(double rangeStart, double rangeEnd, double value) {
        return (value - rangeStart) / (rangeEnd - rangeStart);
    }

    /**
     * Linearly interpolate between two values.
     *
     * @param fromValue
     * @param toValue
     * @param progress
     * @return fromValue + (toValue - fromValue) * progress
     */
    public static double linearInterpolate(double fromValue, double toValue, double progress) {
        return fromValue + (toValue - fromValue) * progress;
    }

    /**
     * Clamp a value in the range [min, max].
     *
     * @param value
     * @param min
     * @param max
     * @return min if value &lt; min, max if value &gt; max, value otherwise
     */
    public static double clamp(double value, double min, double max) {
        if (value < min) return min;
        return Math.min(value, max);
    }

    /**
     * Clamp a value in the range [min, max].
     *
     * @param value
     * @param min
     * @param max
     * @return min if value &lt; min, max if value &gt; max, value otherwise
     */
    public static long clamp(long value, long min, long max) {
        if (value < min) return min;
        return Math.min(value, max);
    }

    /**
     * Determine if a value is a power of two.
     *
     * @param value
     * @return true if value is a power of two, false otherwise
     */
    public static boolean isPowerOfTwo(long value) {
        return value != 0 && (value & value - 1) == 0;
    }

    /**
     * Returns the next power of two.
     * If the number is already a power of two, this function acts as an identity function.
     *
     * @param value
     * @return the next power of two
     */
    public static long nextPowerOfTwo(long value) {
        if (value == 0) return 1;
        value--;
        value |= value >> 1;
        value |= value >> 2;
        value |= value >> 4;
        value |= value >> 8;
        value |= value >> 16;
        return value + 1;
    }

    /**
     * Return a random sign (-1 or 1).
     *
     * @return -1 or 1
     */
    public static long randomSign() {
        return random.nextBoolean() ? 1 : -1;
    }

    /**
     * Return a cached value of sin(x) using Raven's method. The cache is 65KB large.
     *
     * @param x angle in radians, between -2*pi and 2*pi inclusive.
     * @return sin(x)
     */
    public static float fastSin(float x) {
        return FastTrigonometry.sin(x);
    }

    /**
     * Return a cached value of cos(x) using Raven's method. The cache is 65KB large.
     *
     * @param x angle in radians, between -2*pi and 2*pi inclusive.
     * @return cos(x)
     */
    public static float fastCos(float x) {
        return FastTrigonometry.cos(x);
    }

    /**
     * Return the integer cube root of a number.
     * If x &lt; 0, -cbrt(-x) is returned.
     *
     * @param x
     * @return floor(cbrt ( x))
     */
    public static int icbrt(int x) {
        long s, y = 0, b, y2 = 0;

        if (x < 0) return -icbrt(-x);

        for (s = 30; s >= 0; s = s - 3) {
            y2 = 4 * y2;
            y = 2 * y;
            b = 3 * (y2 + y) + 1 << s;
            if (x >= b) {
                x = (int) (x - b);
                y2 = y2 + 2 * y + 1;
                y = y + 1;
            }
        }

        return (int) y;
    }

    /**
     * Return the long integer cube root of a number.
     * If x &lt; 0, -cbrt(-x) is returned.
     *
     * @param x
     * @return floor(cbrt ( x))
     */
    public static long icbrt(long x) {
        long s, y = 0, b, y2 = 0;

        if (x < 0) return -icbrt(-x);

        for (s = 60; s >= 0; s = s - 3) {
            y2 = 4 * y2;
            y = 2 * y;
            b = 3 * (y2 + y) + 1 << s;
            if (x >= b) {
                x = x - b;
                y2 = y2 + 2 * y + 1;
                y = y + 1;
            }
        }

        return y;
    }

    /**
     * Return the short integer cube root of a number.
     * If a &lt; 0, -cbrt(-a) is returned.
     *
     * @param a
     * @return floor(cbrt ( a))
     */
    public static short icbrt(short a) {
        long s, y = 0, b, y2 = 0, x = a;

        if (x < 0) return (short) -icbrt(-x);

        for (s = 15; s >= 0; s = s - 3) {
            y2 = 4 * y2;
            y = 2 * y;
            b = 3 * (y2 + y) + 1 << s;
            if (x >= b) {
                x = x - b;
                y2 = y2 + 2 * y + 1;
                y = y + 1;
            }
        }

        return (short) y;
    }

    /**
     * Compute the integer square root of a number.
     * If x &lt; 0, -isqrt(-x) is returned.
     *
     * @param x
     * @return floor(sqrt ( x))
     */
    public static int isqrt(int x) {
        long m = 0x40000000, y = 0, b, t;

        if (x < 0) return -isqrt(-x);

        while (m != 0) {
            b = y | m;
            y = y >> 1;
            t = (int) (x | ~(x - b)) >> 31;
            x = (int) (x - (b & t));
            y = y | m & t;
            m = m >> 2;
        }

        return (int) y;
    }

    /**
     * Compute the integer square root of a number.
     * If x &lt; 0, -isqrt(-x) is returned.
     *
     * @param x
     * @return floor(sqrt ( x))
     */
    public static int isqrt(long x) {
        long m = 0x4000000000000000L, y = 0, b, t;

        if (x < 0) return -isqrt(-x);

        while (m != 0) {
            b = y | m;
            y = y >> 1;
            t = (int) (x | ~(x - b)) >> 31;
            x = (int) (x - (b & t));
            y = y | m & t;
            m = m >> 2;
        }

        return (int) y;
    }

    /**
     * Compute the integer square root of a number.
     * If x &lt; 0, -isqrt(-x) is returned.
     *
     * @param x
     * @return floor(sqrt ( x))
     */
    public static int isqrt(short x) {
        long m = 0x4000, y = 0, b, t;

        if (x < 0) return -isqrt(-x);

        while (m != 0) {
            b = y | m;
            y = y >> 1;
            t = (int) (x | ~(x - b)) >> 31;
            x = (short) (x - (b & t));
            y = y | m & t;
            m = m >> 2;
        }

        return (int) y;
    }

    /**
     * Compute the value of the integer logarithm in base 10 of a number
     *
     * @param x
     * @return floor(log10 ( x))
     */
    public static int ilog10(int x) {
        int y;
        final int[] table = {0, 9, 99, 999, 9999, 99999, 999999, 9999999, 99999999, 999999999, 0xFFFFFFFF};
        y = 19 * (31 - Integer.numberOfLeadingZeros(x)) >> 6;
        y = y + (table[y + 1] - x >>> 31);
        return y;
    }

    /**
     * Compute the value of integer x^n.
     *
     * @param x
     * @param n
     * @return x^n.
     */
    public static int ipow(int x, int n) {
        int p = x, y = 1;
        while (true) {
            if ((n & 1) != 0)
                y = p * y;
            n = n >> 1;
            if (n == 0)
                return y;
            p = p * p;
        }
    }

    /**
     * Compute the value of integer x^n.
     *
     * @param x
     * @param n
     * @return x^n.
     */
    public static long ipow(long x, long n) {
        long p = x, y = 1;
        while (true) {
            if ((n & 1) != 0)
                y = p * y;
            n = n >> 1;
            if (n == 0)
                return y;
            p = p * p;
        }
    }

    /**
     * Compute the value of integer x^n.
     *
     * @param x
     * @param n
     * @return x^n.
     */
    public static short ipow(short x, short n) {
        int p = x, y = 1;
        while (true) {
            if ((n & 1) != 0)
                y = p * y;
            n >>= 1;
            if (n == 0)
                return (short) y;
            p = p * p;
        }
    }

    /**
     * Break floating-point number down into exponent and mantissa
     *
     * @param value
     * @return A pair of the exponent and the mantissa.
     */
    public static Pair<Integer, Double> frexp(double value) {
        if (value == 0.0) {
            return new Pair<>(0, 0.0);
        }

        if (Double.isNaN(value)) {
            return new Pair<>(-1, Double.NaN);
        }

        if (Double.isInfinite(value)) {
            return new Pair<>(-1, value);
        }

        double mantissa = value;
        int exponent = 0;
        int sign = 1;

        if (mantissa < 0f) {
            sign = -1;
            mantissa = -mantissa;
        }
        while (mantissa < 0.5f) {
            mantissa *= 2.0f;
            exponent -= 1;
        }
        while (mantissa >= 1.0f) {
            mantissa *= 0.5f;
            exponent++;
        }
        mantissa *= sign;
        return new Pair<>(exponent, mantissa);
    }

    /**
     * Compute the value of x^z where x is a double precision
     * floating point number and z is an integer.
     *
     * @param x
     * @param z
     * @return x^z.
     * @throws ArithmeticException in case of overflow.
     */
    public static double pow(double x, int z) {
        int n, e, sign, asign, lx;
        double w, y, s;
        if (x == 0.0) {
            if (z == 0)
                return 1.0;
            else if (z < 0)
                return Double.POSITIVE_INFINITY;
            else
                return 0.0;
        }
        if (z == 0)
            return 1.0;
        if (x < 0.0) {
            asign = -1;
            x = -x;
        } else
            asign = 0;
        if (z < 0) {
            sign = -1;
            n = -z;
        } else {
            sign = 1;
            n = z;
        }

        Pair<Integer, Double> ep = frexp(x);
        lx = ep.first();
        s = ep.second();

        e = (lx - 1) * n;
        if (e == 0 || e > 64 || e < -64) {
            s = (s - 7.0710678118654752e-1) / (s + 7.0710678118654752e-1);
            s = (2.9142135623730950 * s - 0.5 + lx) * z * 1.4426950408889634073599;
        } else {
            s = 1.4426950408889634073599 * e;
        }

        if (s > 7.09782712893383996843E2) {
            throw new ArithmeticException("pow: overflow");
        }
        if (s < -7.09782712893383996843E2)
            return 0.0;
        if ((n & 1) != 0)
            y = x;
        else {
            y = 1.0;
            asign = 0;
        }
        w = x;
        n >>= 1;
        while (n != 0) {
            w = w * w;
            if ((n & 1) != 0)
                y *= w;
            n >>= 1;
        }
        if (asign != 0)
            y = -y;
        if (sign < 0)
            y = 1.0 / y;
        return y;
    }

    /**
     * Compute the value of the Airy Ai function at the specified point.
     *
     * @param x
     * @return Ai(x)
     */
    public static double airyAi(double x) {
        return Airy.airy(x)[0];
    }

    /**
     * Compute the value of the Airy Ai function's first derivative at the specified point.
     *
     * @param x
     * @return Ai'(x)
     */
    public static double airyAip(double x) {
        return Airy.airy(x)[1];
    }

    /**
     * Compute the value of the Airy Bi function at the specified point.
     *
     * @param x
     * @return Bi(x)
     */
    public static double airyBi(double x) {
        return Airy.airy(x)[2];
    }

    /**
     * Compute the value of the Airy Bi function's first derivative at the specified point.
     *
     * @param x
     * @return Bi'(x)
     */
    public static double airyBip(double x) {
        return Airy.airy(x)[3];
    }

    /**
     * Compute the value of the Airy Ai, Ai', Bi and Bi' functions at the specified point.
     *
     * @param x
     * @return a double array of length 4 containing Ai, Ai', Bi and Bi' in that order.
     */
    public static double[] airy(double x) {
        return Airy.airy(x);
    }

    /**
     * Compute the gamma function of x.
     *
     * @param x
     * @return gamma(x)
     */
    public static double gamma(double x) {
        return Gamma.gamma(x);
    }

    /**
     * Compute the logarithm of the gamma function of x.
     *
     * @param x
     * @return log(gamma ( x))
     */
    public static double loggamma(double x) {
        return Gamma.loggamma(x);
    }

    /**
     * Compute the digamma function of x.
     *
     * @param x
     * @return digamma(x)
     */
    public static double digamma(double x) {
        return Gamma.digamma(x);
    }

    /**
     * Compute the trigamma function of x.
     *
     * @param x
     * @return trigamma(x)
     */
    public static double trigamma(double x) {
        return Gamma.trigamma(x);
    }

    /**
     * Compute the value of the upper incomplete gamma function.
     *
     * @param a
     * @param x
     * @return gamma_u(a, x)
     */
    public static double uiGamma(double a, double x) {
        return Gamma.upperIncomplete(a, x);
    }

    /**
     * Compute the value of the lower incomplete gamma function.
     *
     * @param a
     * @param x
     * @return gamma_l(a, x)
     */
    public static double liGamma(double a, double x) {
        return Gamma.lowerIncomplete(a, x);
    }

    /**
     * Compute the Pochhammer symbol (x)_n.
     *
     * @param x
     * @param n
     * @return (x)_n
     */
    public static double pochhammer(double x, double n) {
        return Gamma.gamma(x + n) / Gamma.gamma(x);
    }

    /**
     * Compute the Pochhammer symbol (x)_n.
     *
     * @param x
     * @param n
     * @return (x)_n
     */
    public static double pochhammer(double x, int n) {
        return Gamma.gamma(x + n) / Gamma.gamma(x);
    }

    /**
     * Compute the Pochhammer symbol (x)_n.
     *
     * @param x
     * @param n
     * @return (x)_n
     */
    public static Complex pochhammer(Complex x, int n) {
        return div(Gamma.gamma(add(x, n)), Gamma.gamma(x));
    }

    /**
     * Compute the Pochhammer symbol (x)_n.
     *
     * @param x
     * @param n
     * @return (x)_n
     */
    public static Complex pochhammer(Complex x, double n) {
        return div(Gamma.gamma(add(x, n)), Gamma.gamma(x));
    }

    /**
     * Compute the Pochhammer symbol (x)_n.
     *
     * @param x
     * @param n
     * @return (x)_n
     */
    public static Complex pochhammer(Complex x, Complex n) {
        return div(Gamma.gamma(add(x, n)), Gamma.gamma(x));
    }

    /**
     * Compute the value of the exponential integral at x.
     *
     * @param x
     * @return Ei(x)
     */
    public static double Ei(double x) {
        return Ei.expint(x);
    }

    /**
     * Compute the value of the Riemann zeta function at x.
     *
     * @param x
     * @return zeta(x)
     */
    public static double zeta(double x) {
        return Zeta.riemann_zeta(x);
    }

    /**
     * Compute the value of the complex Riemann zeta function at z.
     *
     * @param z
     * @return zeta(z)
     */
    public static Complex zeta(Complex z) {
        return Zeta.riemann_zeta(z);
    }

    /**
     * Compute the value of the Hurwitz zeta function at x.
     *
     * @param x
     * @param a
     * @return zeta(x, a)
     */
    public static double hurwitzZeta(double x, double a) {
        return Zeta.hurwitz_zeta(x, a);
    }

    /**
     * Compute the value of the n-th polygamma function at x.
     *
     * @param n
     * @param x
     * @return polygamma(n, x)
     */
    public static double polygamma(double n, double x) {
        // Polygamma[n, x] = (-1)^(n+1) * Gamma[n + 1] * HurwitzZeta[n + 1, x]
        return Math.pow(-1, n + 1) * gamma(n + 1) * Zeta.hurwitz_zeta(n + 1, x);
    }

    /**
     * Compute the value of the order-N polygamma function at complex z.
     *
     * @param n
     * @param z
     * @return polygamma(n, z)
     */
    public static Complex polygamma(Complex n, Complex z) {
        Complex np1 = add(n, 1);
        return mul(pow(-1, np1), mul(gamma(np1), Zeta.hurwitz_zeta(np1, z)));
    }

    /**
     * Compute the beta function of two values.
     *
     * @param x
     * @param y
     * @return beta(x, y)
     */
    public static double beta(double x, double y) {
        return Gamma.gamma(x) * Gamma.gamma(y) / Gamma.gamma(x + y);
    }

    /**
     * Compute the logarihtm of the beta function of two values.
     * Uses the identity that log(beta(x, y)) = log(gamma(x)) + log(gamma(y)) - log(gamma(x + y)).
     *
     * @param x
     * @param y
     * @return log(beta ( x, y))
     */
    public static double logbeta(double x, double y) {
        return Gamma.loggamma(x) + Gamma.loggamma(y) - Gamma.loggamma(x + y);
    }

    /**
     * Return the factorial of n as a double-precision. n must be positive.
     * Faster than using the gamma function.
     *
     * @param n
     * @return n!
     * @throws ArithmeticException if n is negative
     */
    public static double factorial(long n) {
        return Gamma.factorial(n);
    }

    /**
     * Compute the dilogarithm (the value of the Spence function at 1-x) of x.
     *
     * @param n
     * @return dilog(x)
     */
    public static double dilog(double n) {
        return Spence.dilog(n);
    }

    /**
     * Compute the Spence function of x.
     *
     * @param n
     * @return Spence(x)
     */
    public static double spence(double n) {
        return Spence.spence(n);
    }

    /**
     * Compute the polylogarithm of x.
     *
     * @param n
     * @param x
     * @return Li_n(x)
     */
    public static double polylog(int n, double x) {
        return Spence.polylog(n, x);
    }

    /**
     * Compute the value of the secant (1 / cos(x)) function at x.
     *
     * @param x
     * @return sec(x)
     */
    public static double sec(double x) {
        return 1.0 / Math.cos(x);
    }

    /**
     * Compute the value of the cosecant (1 / sin(x)) function at x.
     *
     * @param x
     * @return csc(x)
     */
    public static double csc(double x) {
        return 1.0 / Math.sin(x);
    }

    /**
     * Compute the value of the cotangent (1 / tan(x)) function at x.
     *
     * @param x
     * @return cot(x)
     */
    public static double cot(double x) {
        return 1.0 / Math.tan(x);
    }

    /**
     * Compute the value of the hyperbolic secant (1 / cosh(x)) function at x.
     *
     * @param x
     * @return sech(x)
     */
    public static double sech(double x) {
        return 1.0 / Math.cosh(x);
    }

    /**
     * Compute the value of the hyperbolic cosecant (1 / sinh(x)) function at x.
     *
     * @param x
     * @return csch(x)
     */
    public static double csch(double x) {
        return 1.0 / Math.sinh(x);
    }

    /**
     * Compute the value of the hyperbolic cotangent (1 / tanh(x)) function at x.
     *
     * @param x
     * @return coth(x)
     */
    public static double coth(double x) {
        return 1.0 / Math.tanh(x);
    }

    /**
     * Compute the value of the inverse secant (1 / acos(x)) function at x.
     *
     * @param x
     * @return asec(x)
     */
    public static double asec(double x) {
        return Math.acos(1.0 / x);
    }

    /**
     * Compute the value of the inverse cosecant (1 / asin(x)) function at x.
     *
     * @param x
     * @return acsc(x)
     */
    public static double acsc(double x) {
        return Math.asin(1.0 / x);
    }

    /**
     * Compute the value of the inverse cotangent (1 / atan(x)) function at x.
     *
     * @param x
     * @return acot(x)
     */
    public static double acot(double x) {
        return Math.atan(1.0 / x);
    }

    /**
     * Compute the value of the inverse hyperbolic sine function at x.
     *
     * @param a
     * @return asinh(a)
     */
    public static double asinh(double a) {
        final double sign;

        if (Double.doubleToRawLongBits(a) < 0) {
            a = Math.abs(a);
            sign = -1.0d;
        } else {
            sign = 1.0d;
        }

        return sign * Math.log(Math.sqrt(a * a + 1.0d) + a);
    }

    private static double safeLog(double x) {
        if (x == 0.0D) {
            return 0.0D;
        } else {
            return Math.log(x);
        }
    }

    /**
     * Compute the value of the inverse hyperbolic cosine function at x.
     *
     * @param x
     * @return acosh(x)
     */
    public static double acosh(double x) {
        double ans;

        if (Double.isNaN(x) || x < 1) {
            ans = Double.NaN;
        } else if (x < 94906265.62) {
            ans = safeLog(x + Math.sqrt(x * x - 1.0D));
        } else {
            ans = 0.69314718055994530941723212145818D + safeLog(x);
        }

        return ans;
    }

    /**
     * Compute the value of the inverse hyperbolic tangent function at x.
     *
     * @param a
     * @return atanh(a)
     */
    public static double atanh(double a) {
        final double mult;

        if (Double.doubleToRawLongBits(a) < 0) {
            a = Math.abs(a);
            mult = -0.5d;
        } else {
            mult = 0.5d;
        }
        return mult * Math.log((1.0d + a) / (1.0d - a));
    }

    /**
     * Compute the value of the inverse hyperbolic secant (1 / acosh(x)) function at x.
     *
     * @param x
     * @return asech(x)
     */
    public static double asech(double x) {
        return acosh(1.0 / x);
    }

    /**
     * Compute the value of the inverse hyperbolic cosecant (1 / asinh(x)) function at x.
     *
     * @param x
     * @return acsch(x)
     */
    public static double acsch(double x) {
        return asinh(1.0 / x);
    }

    /**
     * Compute the value of the inverse hyperbolic cotangent (1 / atanh(x)) function at x.
     *
     * @param x
     * @return acoth(x)
     */
    public static double acoth(double x) {
        return atanh(1.0 / x);
    }

    /**
     * Compute Lambert W_0 (x).
     *
     * @param x
     * @return W_0(x)
     */
    public static double lambertW0(double x) {
        return Lambert.lambert0(x);
    }

    /**
     * Compute Lambert W_(-1) (x).
     *
     * @param x
     * @return W_(- 1) (x)
     */
    public static double lambertWm1(double x) {
        return Lambert.lambertn1(x);
    }

    /**
     * Compute the value of the Lerch transcendent function at z, s, a.
     *
     * @param z
     * @param s
     * @param a
     * @return Lerch(z, s, a)
     * @throws ArithmeticException if the computation fails unexpectedly due to exceeding the amount of allowed iterations.
     */
    public static double lerchPhi(double z, double s, double a) {
        return Zeta.lerch_phi(z, s, a);
    }

    /**
     * Compute the value of the Dawson function (D+) at x.
     *
     * @param x
     * @return D+(x)
     */
    public static double dawsonPlus(double x) {
        return Erf.dawson(x);
    }

    /**
     * Compute the value of the Dawson function (D-) at x.
     *
     * @param x
     * @return D-(x)
     */
    public static double dawsonMinus(double x) {
        return Erf.dawsonm(x);
    }

    /**
     * Compute the value of the error function at x.
     *
     * @param x
     * @return erf(x)
     */
    public static double erf(double x) {
        return Erf.erf(x);
    }

    /**
     * Compute the value of the complementary error function at x.
     *
     * @param x
     * @return erfc(x)
     */
    public static double erfc(double x) {
        return Erf.erfc(x);
    }

    /**
     * Compute the value of the imaginary error function at x.
     *
     * @param x
     * @return erfi(x)
     */
    public static double erfi(double x) {
        return Erf.erfi(x);
    }

    /**
     * Compute the value of the inverse of the logistic sigmoid "squash" function at x.
     *
     * @param x
     * @return stretch(x)
     */
    public static double stretch(double x) {
        return Math.log(x / (1.0 - x));
    }

    /**
     * Compute the value of the logistic sigmoid "squash" function at x.
     *
     * @param x
     * @return squash(x)
     */
    public static double squash(double x) {
        return 1.0 / (1.0 + Math.exp(-x));
    }

    /**
     * Compute the value of the sine integral function at x.
     *
     * @param x
     * @return Si(x)
     */
    public static double Si(double x) {
        return TrigonometricIntegral.Si(x);
    }

    /**
     * Compute the value of the cosine integral function at x.
     *
     * @param x
     * @return Ci(x)
     */
    public static double Ci(double x) {
        return TrigonometricIntegral.Ci(x);
    }

    /**
     * Compute the value of the hyperbolic sine integral function at x.
     *
     * @param x
     * @return
     */
    public static double Shi(double x) {
        return TrigonometricIntegral.Shi(x);
    }

    /**
     * Compute the value of the hyperbolic cosine integral function at x.
     *
     * @param x
     * @return Chi(x)
     */
    public static double Chi(double x) {
        return TrigonometricIntegral.Chi(x);
    }

    /**
     * Compute the value of the sine integral function at x.
     *
     * @param x
     * @return si(x)
     */
    public static double si(double x) {
        return TrigonometricIntegral.si(x);
    }

    /**
     * Compute the value of the cosine integral function at x.
     *
     * @param x
     * @return Cin(x)
     */
    public static double Cin(double x) {
        return TrigonometricIntegral.Cin(x);
    }

    /**
     * Compute both the hyperbolic sine and cosine integral at x.
     * Return the results in a two-element double precision array.
     *
     * @param x
     * @return {Shi(x), Chi(x)}
     */
    public static double[] ShiChi(double x) {
        return TrigonometricIntegral.ShiChi(x);
    }

    /**
     * Compute the Fresnel integral C(x).
     *
     * @param x
     * @return C(x)
     */
    public static double fresnelC(double x) {
        return Fresnel.fresnelC(x);
    }

    /**
     * Compute the Fresnel integral S(x).
     *
     * @param x
     * @return S(x)
     */
    public static double fresnelS(double x) {
        return Fresnel.fresnelS(x);
    }

    /**
     * Compute the bessel Y0 function at x.
     *
     * @param x
     * @return Y0(x)
     */
    public static double besselY0(double x) {
        return Bessel.y0(x);
    }

    /**
     * Compute the bessel Y1 function at x.
     *
     * @param x
     * @return Y1(x)
     */
    public static double besselY1(double x) {
        return Bessel.y1(x);
    }

    /**
     * Compute the bessel Yn function at x.
     *
     * @param n
     * @param x
     * @return Yn(x)
     */
    public static double besselYn(int n, double x) {
        return Bessel.yn(n, x);
    }

    /**
     * Compute the bessel J0 function at x.
     *
     * @param x
     * @return J0(x)
     */
    public static double besselJ0(double x) {
        return Bessel.bessel0(x);
    }

    /**
     * Compute the bessel J1 function at x.
     *
     * @param x
     * @return J1(x)
     */
    public static double besselJ1(double x) {
        return Bessel.bessel1(x);
    }

    /**
     * Compute the bessel Jn function at x.
     *
     * @param n
     * @param x
     * @return Jn(x)
     */
    public static double besselJn(int n, double x) {
        return Bessel.bessel(n, x);
    }

    /**
     * Compute the bessel I0 function at x.
     *
     * @param x
     * @return I0(x)
     */
    public static double besselI0(double x) {
        return Bessel.i0(x);
    }

    /**
     * Compute the bessel I1 function at x.
     *
     * @param x
     * @return I1(x)
     */
    public static double besselI1(double x) {
        return Bessel.i1(x);
    }

    /**
     * Compute the bessel K0 function at x.
     *
     * @param x
     * @return K0(x)
     */
    public static double besselK0(double x) {
        return Bessel.k0(x);
    }

    /**
     * Compute the bessel K1 function at x.
     *
     * @param x
     * @return K1(x)
     */
    public static double besselK1(double x) {
        return Bessel.k1(x);
    }

    /**
     * Compute the bessel Kn function at x.
     *
     * @param n
     * @param x
     * @return Kn(x)
     */
    public static double besselKn(int n, double x) {
        return Bessel.kn(n, x);
    }

    /**
     * Compute the greatest common divisor of two integers.
     *
     * @param a
     * @param b
     * @return gcd(a, b)
     */
    public static long gcd(long a, long b) {
        long result;
        if (a == 0) {
            result = b;
        } else if (b == 0) {
            result = a;
        } else {
            long r = a % b;
            while (r != 0) {
                a = b;
                b = r;
                r = a % b;
            }
            result = b;
        }
        return result;
    }

    /**
     * Domain extension of the greatest common division function onto the real line.
     *
     * @param a
     * @param b
     * @return gcd(a, b)
     */
    public static double gcd(double a, double b) {
        double result;
        if (a == 0) {
            result = b;
        } else if (b == 0) {
            result = a;
        } else {
            double r = a % b;
            while (r != 0) {
                a = b;
                b = r;
                r = a % b;
            }
            result = b;
        }
        return result;
    }

    /**
     * Domain extension of the greatest common division function onto the complex plane.
     *
     * @param a
     * @param b
     * @return gcd(a, b)
     */
    public static Complex gcd(Complex a, Complex b) {
        Complex result;
        if (eq(a, Complex.ZERO)) {
            result = b;
        } else if (eq(b, Complex.ZERO)) {
            result = a;
        } else {
            Complex r = rem(a, b);
            while (ne(r, 0)) {
                a = b;
                b = r;
                r = rem(a, b);
            }
            result = b;
        }
        return result;
    }

    /**
     * Domain extension of the remainder function onto the complex plane.
     *
     * @param a
     * @param b
     * @return a rem b
     */
    public static Complex rem(Complex a, Complex b) {
        Complex quot = div(a, b);
        quot = new Complex(floor(quot.getReal()), floor(quot.getImaginary()));
        return sub(a, mul(b, quot));
    }

    /**
     * Compute the least common multiple of two integers.
     *
     * @param a
     * @param b
     * @return lcm(a, b)
     */
    public static long lcm(long a, long b) {
        if (a == 0 || b == 0) return 0;
        return a * b / gcd(a, b);
    }

    /**
     * Domain extension of the least common multiple function onto the real line.
     *
     * @param a
     * @param b
     * @return lcm(a, b)
     */
    public static double lcm(double a, double b) {
        if (a == 0 || b == 0) return 0;
        return a * b / gcd(a, b);
    }

    /**
     * Domain extension of the least common multiple function onto the complex plane.
     *
     * @param a
     * @param b
     * @return lcm(a, b)
     */
    public static Complex lcm(Complex a, Complex b) {
        if (eq(a, Complex.ZERO) || eq(b, Complex.ZERO)) return Complex.ZERO;
        return div(mul(a, b), gcd(a, b));
    }

    /**
     * Compute the n-th Fibonacci number using Binet's formula.
     * Due to floating point precision issues, this method is only
     * accurate for n &le;= 75.
     *
     * @param a
     * @return fib(a)
     */
    public static long fib(int a) {
        double phi = 1.6180339887498948482045868343656381177203091798057628621354486227;
        double psi = -0.618033988749894848204586834365638117720309179805762862135448622;
        double sqrt5 = 2.2360679774997896964091736687312762354406183596115257242708972454;
        return Math.round((pow(phi, a) - pow(psi, a)) / sqrt5);
    }

    /**
     * Computes the Gaussian hypergeometric function (2F1) of four arguments.
     *
     * @param a
     * @param b
     * @param c
     * @param x
     * @return (2F1)(a, b, c, x)
     */
    public static double hypergeo2F1(double a, double b, double c, double x) {
        return Hypergeometric.hyp2f1(a, b, c, x);
    }

    /**
     * Computes the Confluent hypergeometric function (1F1) of three arguments.
     *
     * @param a
     * @param b
     * @param x
     * @return (1F1)(a, b, x)
     */
    public static double hypergeo1F1(double a, double b, double x) {
        return Hypergeometric.hyperg(a, b, x);
    }

    /**
     * Computes the 1F2 case of the generalised hypergeometric function.
     *
     * @param a
     * @param b
     * @param c
     * @param x
     * @return a double array the approximation and the estimated error of (1F2)(a, b, c x)
     */
    public static double hypergeo1F2(double a, double b, double c, double x) {
        Hypergeometric.DoublePtr result = new Hypergeometric.DoublePtr();
        return Hypergeometric.hypergeo1f2(a, b, c, x, result);
    }

    /**
     * Computes the 3F0 case of the generalised hypergeometric function.
     *
     * @param a
     * @param b
     * @param c
     * @param x
     * @return a double array the approximation and the estimated error of (1F2)(a, b, c x)
     */
    public static double hypergeo3F0(double a, double b, double c, double x) {
        Hypergeometric.DoublePtr result = new Hypergeometric.DoublePtr();
        return Hypergeometric.hypergeo3f0(a, b, c, x, result);
    }

    /**
     * Computes the value of the Struve function of order v at x.
     *
     * @param v order of the Struve function
     * @param x argument of the Struve function
     * @return Struve(v, x)
     */
    public static double struve(double v, double x) {
        if (v == 0.0) return Hypergeometric.struveH0(x);
        else if (v == 1.0) return Hypergeometric.struveH1(x);
        else return Hypergeometric.struve(v, x);
    }

    /**
     * Computes the value of the modified Struve function of 0th order.
     *
     * @param x argument of the Struve function
     * @return StruveL0(x)
     */
    public static double struveL0(double x) {
        return Hypergeometric.struveL0(x);
    }

    /**
     * Computes the value of the modified Struve function of 1st order.
     *
     * @param x argument of the Struve function
     * @return StruveL1(x)
     */
    public static double struveL1(double x) {
        return Hypergeometric.struveL1(x);
    }

    /**
     * Computes the logarithm of the absolute value of the gamma function of x.
     *
     * @param x
     * @return { log(|Gamma(x)|), sign(Gamma(x)) }
     * @see Maja#loggamma(double)
     */
    public static double[] logabsgamma(double x) {
        return Gamma.lgam(x);
    }

    /**
     * Compute fractional order bessel function of n and x.
     *
     * @param n
     * @param x
     * @return J_n(x)
     */
    public static double besselJv(double n, double x) {
        return Bessel.jv(n, x);
    }

    /**
     * Compute the fractional order bessel Y function of v and x.
     *
     * @param v
     * @param x
     * @return Y_n(x)
     */
    public static double besselYv(double v, double x) {
        return Bessel.yv(v, x);
    }

    /**
     * Integrate a monadic function over a finite interval [a,b] using the
     * Simpson rule. The number of intervals is given by N and the precision
     * of the final result greatly depends on this parameter.
     *
     * @param f function to integrate
     * @param a lower bound
     * @param b upper bound
     * @param N number of intervals, N=10000 tends to give a good approximation in most scenarios.
     * @return integral of f over [a,b]
     * @throws IllegalArgumentException if N is not a positive integer
     */
    public static double integrateSimpsonReal(Function<Double, Double> f, double a, double b, int N) {
        // Properly handle the configurations of a and b.
        if (a < b)
            return Integrator.finiteSimpsonRR(f, a, b, N);
        else if (a == b)
            return 0.0;
        else
            return -Integrator.finiteSimpsonRR(f, b, a, N);
    }

    /**
     * Integrate a monadic function over a finite interval [a,b] using the
     * Gauss-Legendre quadrature. The number of intervals is given by N and the precision
     * of the final result greatly depends on this parameter.
     * The computation of an integral using the Gauss-Legendre quadrature involves caching the
     * coefficients required to perform the integration depending on the value of the N parameter.
     * This means that the first call to this method will be slower than subsequent calls with the
     * same value of N. The coefficients are internally cached inside a ConcurrentHashMap.
     *
     * @param f function to integrate
     * @param a lower bound
     * @param b upper bound
     * @param N number of intervals, N=6 tends to give a good approximation in most scenarios.
     *          N must be between 1 and 30.
     * @return integral of f over [a,b]
     * @throws IllegalArgumentException if N is not a positive integer
     */
    public static double integrateGaussLegendreReal(Function<Double, Double> f, double a, double b, int N) {
        if (a < b)
            return Integrator.gaussLegendreIntegrateRR(f, a, b, N);
        else if (a == b)
            return 0.0;
        else
            return -Integrator.gaussLegendreIntegrateRR(f, b, a, N);
    }

    /**
     * Integrate a monadic function over a finite interval [a,b] using the
     * Tanh-Sinh quadrature, especially useful when singularities or infinite
     * derivatives exist at one or both endpoints. The Tanh-Sinh quadrature is
     * not as efficient as Gaussian quadrature for smooth integrands.
     *
     * @param f   function to integrate
     * @param a   lower bound
     * @param b   upper bound
     * @param N   the degree of the quadrature, usually N=6 is sufficient
     * @param eps desired precision of the result (usually 1.0e-9 is sufficient)
     * @return an array of double values, first of which is the integral of f over [a,b],
     * while the second is the estimated error.
     * @throws IllegalArgumentException if N is not a positive integer
     */
    public static double[] integrateTanhSinhReal(Function<Double, Double> f, double a, double b, int N, double eps) {
        if (a < b)
            return Integrator.finiteTanhSinhRR(f, a, b, N, eps);
        else if (a == b)
            return new double[]{0.0, 0.0};
        else {
            double[] res = Integrator.finiteTanhSinhRR(f, b, a, N, eps);
            res[0] = -res[0];
            return res;
        }
    }

    /**
     * Integrate a monadic function over a finite interval [a,b] using the
     * Simpson rule. The number of intervals is given by N and the precision
     * of the final result greatly depends on this parameter.
     *
     * @param f function to integrate
     * @param a lower bound
     * @param b upper bound
     * @param N number of intervals, N=10000 tends to give a good approximation in most scenarios.
     * @return integral of f over [a,b]
     * @throws IllegalArgumentException if N is not a positive integer
     */
    public static Complex integrateSimpsonRC(Function<Double, Complex> f, double a, double b, int N) {
        // Properly handle the configurations of a and b.
        if (a < b)
            return Integrator.finiteSimpsonRC(f, a, b, N);
        else if (a == b)
            return Complex.ZERO;
        else
            return negate(Integrator.finiteSimpsonRC(f, b, a, N));
    }

    /**
     * Integrate a monadic function over a finite interval [a,b] using the
     * Gauss-Legendre quadrature. The number of intervals is given by N and the precision
     * of the final result greatly depends on this parameter.
     * The computation of an integral using the Gauss-Legendre quadrature involves caching the
     * coefficients required to perform the integration depending on the value of the N parameter.
     * This means that the first call to this method will be slower than subsequent calls with the
     * same value of N. The coefficients are internally cached inside a ConcurrentHashMap.
     *
     * @param f function to integrate
     * @param a lower bound
     * @param b upper bound
     * @param N number of intervals, N=6 tends to give a good approximation in most scenarios.
     *          N must be between 1 and 30.
     * @return integral of f over [a,b]
     * @throws IllegalArgumentException if N is not a positive integer
     */
    public static Complex integrateGaussLegendreRC(Function<Double, Complex> f, double a, double b, int N) {
        if (a < b)
            return Integrator.gaussLegendreIntegrateRC(f, a, b, N);
        else if (a == b)
            return Complex.ZERO;
        else
            return negate(Integrator.gaussLegendreIntegrateRC(f, b, a, N));
    }

    /**
     * Integrate a monadic function over a finite interval [a,b] using the
     * Tanh-Sinh quadrature, especially useful when singularities or infinite
     * derivatives exist at one or both endpoints. The Tanh-Sinh quadrature is
     * not as efficient as Gaussian quadrature for smooth integrands.
     *
     * @param f   function to integrate
     * @param a   lower bound
     * @param b   upper bound
     * @param N   the degree of the quadrature, usually N=6 is sufficient
     * @param eps desired precision of the result (usually 1.0e-9 is sufficient)
     * @return an array of double values, first of which is the integral of f over [a,b],
     * while the second is the estimated error.
     * @throws IllegalArgumentException if N is not a positive integer
     */
    public static Complex[] integrateTanhSinhRC(Function<Double, Complex> f, double a, double b, int N, double eps) {
        if (a < b)
            return Integrator.finiteTanhSinhRC(f, a, b, N, eps);
        else if (a == b)
            return new Complex[]{Complex.ZERO, Complex.ZERO};
        else {
            Complex[] result = Integrator.finiteTanhSinhRC(f, b, a, N, eps);
            result[0] = negate(result[0]);
            return result;
        }
    }

    /**
     * Compute the binomial coefficient "n choose k".
     *
     * @param n the number of elements, n &gt; 0.
     * @param k the number of elements to choose, 0 &lt; k &lt;= n.
     * @return n! / (k! * (n-k)!)
     * @throws IllegalArgumentException if n &lt;= 0 or k &lt; 0 or k &gt; n.
     */
    public static long binomial(int n, int k) {
        if (n <= 0 || k < 0 || k > n)
            throw new IllegalArgumentException("Invalid arguments: n = " + n + ", k = " + k);
        // Naive method.
        if (k > n - k)
            k = n - k;
        long b = 1;
        for (int i = 1, m = n; i <= k; i++, m--)
            b = b * m / i;
        return b;
    }

    /**
     * Find a root of a monadic function using the Newton-Raphson method.
     *
     * @param f   the function to find a root for
     * @param df  the derivative of the function
     * @param x0  the initial guess
     * @param eps the desired precision of the result
     * @return a root of the function f within the desired precision
     * unless the iteration limit is exceeded.
     */
    public static double newtonRaphson(Function<Double, Double> f, Function<Double, Double> df, double x0, double eps) {
        return Root.newtonRaphson(f, df, x0, eps);
    }

    /**
     * Numerically finds all roots of a polynomial equation P(z) = 0.
     * There are (n+1) coefficients and n roots. The leading coefficient
     * must not be zero. Errors are returned in a boolean array for
     * each root. Success is indicated by false, i.e. no error. On return
     * each root should be checked against its flag.
     * @param coefficients
     * @return
     */
    public static boolean[] aberth(Complex[] coefficients) {
        double[] re, im;
        int n = coefficients.length;
        if (n < 1)
            throw new IllegalArgumentException("Invalid number of coefficients: " + n);
        else if (n == 2) {
            // Linear case.
            Complex a = coefficients[1];
            Complex b = coefficients[0];
            if (eq(a, Complex.ZERO)) {
                // No roots.
                return new boolean[]{ true };
            } else {
                // One root.
                Complex x = div(negate(b), a);
                coefficients[0] = x;
                return new boolean[]{ false };
            }
        }
        else if(n == 3) {
            // Quadratic case.
            Complex a = coefficients[2];
            Complex b = coefficients[1];
            Complex c = coefficients[0];
            // Discriminant.
            Complex d = sqrt(sub(mul(b, b), mul(a, mul(c, 4))));
            if (eq(d, Complex.ZERO)) {
                // One root.
                Complex x = div(negate(b), mul(a, 2));
                coefficients[0] = coefficients[1] = x;
                return new boolean[]{ false, false };
            } else {
                // Two roots.
                Complex x1 = div(add(negate(b), d), mul(a, 2));
                Complex x2 = div(sub(negate(b), d), mul(a, 2));
                coefficients[0] = x1;
                coefficients[1] = x2;
                return new boolean[]{ false, false };
            }
        }
        if (coefficients[n - 1].equals(Complex.ZERO))
            throw new IllegalArgumentException("Leading coefficient must not be zero: " + coefficients[n]);
        re = new double[n];
        im = new double[n];
        for (int i = 0; i < n; i++) {
            re[i] = coefficients[i].getReal();
            im[i] = coefficients[i].getImaginary();
        }
        boolean[] err = new boolean[n];
        Arrays.fill(err, true);
        Pzeros.aberth(re, im, err);
        for (int i = 0; i < n; i++) {
            if(abs(im[i]) > EPSILON)
                coefficients[i] = new Complex(re[i], im[i]);
            else
                coefficients[i] = new Complex(re[i], 0);
        }
        return err;
    }

    /**
     * Returns an integer b such that |B_n| &le; 2^b for all Bernoulli numbers B_n.
     * @param n
     * @return
     */
    public static int bernoulliBound(int n) {
        if (n % 2 == 1 || n <= 13) return 0; // 2^0 = 1, B_n up until 13 are |B_n| <= 1.
        // str:implode-on "," \:to-string \ceil (log2 (filter #0 (abs (:bernoulli (range 14 80)))))
        final int lut[] = { 1,3,6,10,13,17,21,25,30,34,39,44,49,55,60,66,71,77,83,89,95,102,108,115,121,128,135,141,148,155,162,170,177 };
        // lut[0] = B_14, lut[1] = B_16, lut[2] = B_18, ...
        int i = (n - 14) / 2;
        if (i < lut.length) return lut[i];
        else {
            // |B_n| < 4n!/(2pi)^n < 4(n+1)^(n+1)e^(-n)/(2pi)^n (per Stirling)
            // trivially, log2(4(n+1)^(n+1)e^(-n)/(2pi)^n)=log2(4(n+1)^(n+1))-n*log2(e)-n*log2(2pi)
            return (int) ceil((n+1)*log(n+1)/log(2)+2 - n * log2(E) - n * log2(2 * PI));
        }
    }

    /**
     * Computes the Arithmetic-Geometric mean.
     * @param a
     * @param g
     * @return the Arithmetic-Geometric mean of a and g
     */
    public static double agm(double a, double g) {
        double a1 = a;
        double g1 = g;
        while (Math.abs(a1 - g1) >= 1.0e-14) {
            double a2 = (a1 + g1) / 2.0;
            double g2 = Math.sqrt(a1 * g1);
            a1 = a2;
            g1 = g2;
        }
        return a1;
    }

    // TODO: Document.

    public static double fac2(double a) {
        // x!! = 2^(x/2) (pi/2)^((cos(pi x)-1)/4) gamma(x/2+1).
        double x = pow(2, a/2);
        double b = pow(PI_2, (cos(PI * a) - 1) / 4);
        double c = gamma((a + 1) / 2);
        return x * b * c;
    }

    public static double angerJ(double v, double z) {
        return 0.3183098861837906715377675267450287240689192914809128974953346881
                * integrateTanhSinhReal(x -> cos(v * x - z * sin(x)), 0, PI, 7, 1.0e-11)[0];
    }

    public static double weberE(double v, double z) {
        return 0.3183098861837906715377675267450287240689192914809128974953346881
                * integrateTanhSinhReal(x -> sin(v * x - z * sin(x)), 0, PI, 7, 1.0e-11)[0];
    }

    public static Complex angerJ(Complex v, Complex z) {
        return mul(0.3183098861837906715377675267450287240689192914809128974953346881,
                 integrateTanhSinhRC(x -> cos(sub(mul(v, x), mul(z, sin(x)))), 0, PI, 7, 1.0e-11)[0]);
    }

    public static Complex weberE(Complex v, Complex z) {
        return mul(0.3183098861837906715377675267450287240689192914809128974953346881,
                integrateTanhSinhRC(x -> sin(sub(mul(v, x), mul(z, sin(x)))), 0, PI, 7, 1.0e-11)[0]);
    }

    public static double lommels1(double u, double v, double z) {
        // lommels1(u,v,z) = (bessely(v,z)*int(t -> t**u*besselj(v,t), 0, z) - besselj(v,z)*int(t -> t**u*bessely(v,t), 0, z))*(pi/2)
        double by = besselYv(v, z);
        double byint = integrateTanhSinhReal(x -> pow(x, u) * besselJv(v, x), 0, z, 7, 1.0e-11)[0];
        double bj = besselJv(v, z);
        double bjint = integrateTanhSinhReal(x -> pow(x, u) * besselYv(v, x), 0, z, 7, 1.0e-11)[0];
        return (by * byint - bj * bjint) * PI_2;
    }

    public static double legendre(double n, double z) {
        // Legendre polynomials
        // P_n(z) = 2F1(-n,n+1,1,(1-z)/2)
        return hypergeo2F1(-n, n + 1, 1, (1 - z) / 2);
    }

    /**
     * Add two complex numbers together.
     *
     * @param a
     * @param b
     * @return a + b
     */
    public static Complex add(Complex a, Complex b) {
        return new Complex(a.getReal() + b.getReal(), a.getImaginary() + b.getImaginary());
    }

    /**
     * Perform dual addition.
     * @param a
     * @param b
     * @return a + b
     */
    public static Dual add(Dual a, Dual b) {
        return new Dual(a.a() + b.a(), a.b() + b.b());
    }

    /**
     * Perform real-dual addition.
     * @param a
     * @param b
     * @return a + b
     */
    public static Dual add(Dual a, double b) {
        return new Dual(a.a() + b, a.b());
    }

    /**
     * Perform real-dual addition.
     * @param a
     * @param b
     * @return a + b
     */
    public static Dual add(double a, Dual b) {
        return new Dual(b.a() + a, b.b());
    }

    public static Dual sub(Dual a, Dual b) {
        return new Dual(a.a() - b.a(), a.b() - b.b());
    }

    public static Dual sub(Dual a, double b) {
        return new Dual(a.a() - b, a.b());
    }

    public static Dual sub(double a, Dual b) {
        return new Dual(a - b.a(), -b.b());
    }

    public static Dual mul(Dual a, Dual b) {
        return new Dual(a.a() * b.a(), a.a() * b.b() + a.b() * b.a());
    }

    public static Dual div(Dual a, Dual b) {
        if(b.a() == 0)
            throw new ArithmeticException("Division by zero");
        // (a+bε)/(c+dε) = a/c + (bc-ad)/(c^2) ε
        double re = a.a() / b.a();
        double du = (a.b() * b.a() - a.a() * b.b()) / (b.a() * b.a());
        return new Dual(re, du);
    }

    public static Dual sin(Dual a) {
        // sin(a+bε) = sin(a) + cos(a)bε
        return new Dual(sin(a.a()), cos(a.a()) * a.b());
    }

    public static Dual sec(Dual a) {
        // sec(a+bε) = sec(a) + tan(a)sec(a)bε
        return new Dual(sec(a.a()), tan(a.a()) * sec(a.a()) * a.b());
    }

    public static Dual cos(Dual a) {
        // cos(a+bε) = cos(a) - sin(a)bε
        return new Dual(cos(a.a()), -sin(a.a()) * a.b());
    }

    public static Dual csc(Dual a) {
        // csc(a+bε) = csc(a) - cot(a)csc(a)bε
        return new Dual(csc(a.a()), -cot(a.a()) * csc(a.a()) * a.b());
    }

    public static Dual tan(Dual a) {
        // tan(a+bε) = tan(a) + b*sec^2(a)ε
        double re = tan(a.a());
        double ca = cos(a.a());
        double du = a.b() / (ca * ca);
        return new Dual(re, du);
    }

    public static Dual cot(Dual a) {
        // cot(a+bε) = cot(a) - b*sec^2(a)ε
        double re = tan(a.a());
        double ca = cos(a.a());
        double du = -a.b() / (ca * ca);
        return new Dual(re, du);
    }

    public static Dual exp(Dual a) {
        // exp(a+bε) = exp(a) + exp(a)bε
        double ea = exp(a.a());
        return new Dual(ea, ea * a.b());
    }

    public static Dual log(Dual a) {
        // log(a+bε) = log(a) + bε/a
        if(a.a() <= 0)
            throw new ArithmeticException("Domain error.");
        double re = log(a.a());
        double du = a.b() / a.a();
        return new Dual(re, du);
    }

    public static Dual pow(Dual a, int b) {
        if(a.a() == 0)
            throw new ArithmeticException("Domain error.");
        double re = pow(a.a(), b);
        double du = b * pow(a.a(), b - 1) * a.b();
        return new Dual(re, du);
    }

    public static Dual abs(Dual a) {
        // abs(a+bε) = |a| + sgn(a)bε
        double re = abs(a.a());
        double du = signum(a.a()) * a.b();
        return new Dual(re, du);
    }

    public static Dual asin(Dual a) {
        // asin(a+bε) = asin(a) + b/sqrt(1-a^2)ε
        double re = asin(a.a());
        double du = a.b() / sqrt(1 - a.a() * a.a());
        return new Dual(re, du);
    }

    public static Dual acos(Dual a) {
        // acos(a+bε) = acos(a) - b/sqrt(1-a^2)ε
        double re = acos(a.a());
        double du = -a.b() / sqrt(1 - a.a() * a.a());
        return new Dual(re, du);
    }

    public static Dual atan(Dual a) {
        // atan(a+bε) = atan(a) + b/(1+a^2)ε
        double re = atan(a.a());
        double du = a.b() / (1 + a.a() * a.a());
        return new Dual(re, du);
    }

    public static Dual sqrt(Dual a) {
        // sqrt(a+bε) = sqrt(a) + b/(2*sqrt(a))ε
        double re = sqrt(a.a());
        double du = a.b() / (2 * sqrt(a.a()));
        return new Dual(re, du);
    }

    public static Dual sinh(Dual a) {
        // sinh(a+bε) = sinh(a) + cosh(a)bε
        double re = sinh(a.a());
        double du = cosh(a.a()) * a.b();
        return new Dual(re, du);
    }

    public static Dual cosh(Dual a) {
        // cosh(a+bε) = cosh(a) + sinh(a)bε
        double re = cosh(a.a());
        double du = sinh(a.a()) * a.b();
        return new Dual(re, du);
    }

    public static Dual tanh(Dual a) {
        // tanh(a+bε) = tanh(a) + b*sech^2(a)ε
        double re = tanh(a.a());
        double ch = cosh(a.a());
        double du = a.b() / (ch * ch);
        return new Dual(re, du);
    }

    public static Dual gamma(Dual a) {
        // gamma'(f(x)) = f'(x)Γ(f(x))ψ(f(x)):
        // gamma(a+bε) = gamma(a) + gamma(a)ψ(a)bε
        double re = gamma(a.a());
        double du = gamma(a.a()) * digamma(a.a()) * a.b();
        return new Dual(re, du);
    }

    public static Dual digamma(Dual a) {
        // digamma(a+bε) = digamma(a) + trigamma(a)bε
        double re = digamma(a.a());
        double du = trigamma(a.a()) * a.b();
        return new Dual(re, du);
    }

    public static Dual besselI0(Dual a) {
        // i0(a+bε) = i0(a) + i1(a)bε
        double re = besselI0(a.a());
        double du = besselI1(a.a()) * a.b();
        return new Dual(re, du);
    }

    public static Dual besselJ0(Dual a) {
        // j0(a+bε) = j0(a) - j1(a)bε
        double re = besselJ0(a.a());
        double du = -besselJ1(a.a()) * a.b();
        return new Dual(re, du);
    }

    public static Dual besselK0(Dual a) {
        // k0(a+bε) = k0(a) - k1(a)bε
        double re = besselK0(a.a());
        double du = -besselK1(a.a()) * a.b();
        return new Dual(re, du);
    }

    public static Dual besselJ1(Dual a) {
        // j1(a+bε) = j1(a) + bε(j0(a) - j2(a))/2
        double re = besselJ1(a.a());
        double du = a.b() * (besselJ0(a.a()) - besselJn(2, a.a())) / 2;
        return new Dual(re, du);
    }

    public static Dual besselY0(Dual a) {
        // y0(a+bε) = y0(a) - y1(a)bε
        double re = besselY0(a.a());
        double du = -besselY1(a.a()) * a.b();
        return new Dual(re, du);
    }

    public static Dual besselK1(Dual a) {
        // k1(a+bε) = k1(a) - bε(k0(a) + k2(a))/2
        double re = besselK1(a.a());
        double du = a.b() * -(besselK0(a.a()) + besselKn(2, a.a())) / 2;
        return new Dual(re, du);
    }

    public static Dual besselJn(int n, Dual a) {
        double re = besselJn(n, a.a());
        double du = a.b() * (besselJn(n - 1, a.a()) - besselJn(n + 1, a.a())) / 2;
        return new Dual(re, du);
    }

    public static Dual besselYn(int n, Dual a) {
        double re = besselYn(n, a.a());
        double du = a.b() * (besselYn(n - 1, a.a()) - besselYn(n + 1, a.a())) / 2;
        return new Dual(re, du);
    }

    public static Dual besselYv(double n, Dual a) {
        double re = besselYv(n, a.a());
        double du = a.b() * (besselYv(n - 1, a.a()) - besselYv(n + 1, a.a())) / 2;
        return new Dual(re, du);
    }

    public static Dual besselKn(int n, Dual a) {
        double re = besselKn(n, a.a());
        double du = a.b() * -(besselKn(n - 1, a.a()) + besselKn(n + 1, a.a())) / 2;
        return new Dual(re, du);
    }

    public static Dual besselJv(double n, Dual a) {
        double re = besselJv(n, a.a());
        double du = a.b() * (besselJv(n - 1, a.a()) - besselJv(n + 1, a.a())) / 2;
        return new Dual(re, du);
    }

    public static Dual airyAi(Dual a) {
        double[] r = airy(a.a());
        return new Dual(r[0], r[1] * a.b());
    }

    public static Dual airyBi(Dual a) {
        double[] r = airy(a.a());
        return new Dual(r[2], r[3] * a.b());
    }

    public static Dual Chi(Dual a) {
        // Chi(a+bε) = Chi(a) + cosh(a)bε/a
        double re = Chi(a.a());
        double du = cosh(a.a()) * a.b() / a.a();
        return new Dual(re, du);
    }

    public static Dual Shi(Dual a) {
        // Shi(a+bε) = Shi(a) + sinh(a)bε/a
        double re = Shi(a.a());
        double du = sinh(a.a()) * a.b() / a.a();
        return new Dual(re, du);
    }

    public static Dual Ci(Dual a) {
        double re = Ci(a.a());
        double du = cos(a.a()) * a.b() / a.a();
        return new Dual(re, du);
    }

    public static Dual Si(Dual a) {
        double re = Si(a.a());
        double du = sin(a.a()) * a.b() / a.a();
        return new Dual(re, du);
    }

    public static Dual Ei(Dual a) {
        double re = Ei(a.a());
        double du = exp(a.a()) * a.b() / a.a();
        return new Dual(re, du);
    }

    public static Dual erf(Dual a) {
        double re = erf(a.a());
        double du = 2 * exp(-a.a() * a.a()) * a.b() / sqrt(PI);
        return new Dual(re, du);
    }

    public static Dual dilog(Dual a) {
        // dilog(a+bε) = dilog(a) + -b*log(1-a)ε/a
        return new Dual(dilog(a.a()), -a.b() * log(1 - a.a()) / a.a());
    }

    public static Dual acot(Dual a) {
        // acot(a+bε) = acot(a) - b/(1+a^2)ε
        double re = acot(a.a());
        double du = -a.b() / (1 + a.a() * a.a());
        return new Dual(re, du);
    }

    public static Dual acoth(Dual a) {
        // acoth(a+bε) = acoth(a) - b/(1-a^2)ε
        double re = acoth(a.a());
        double du = -a.b() / (1 - a.a() * a.a());
        return new Dual(re, du);
    }

    public static Dual asec(Dual a) {
        // asec(a+bε) = asec(a) + b/(a^2 sqrt(1-1/(a^2)))ε
        double re = asec(a.a());
        double asq = a.a() * a.a();
        double du = a.b() / (asq * sqrt(1 - 1 / asq));
        return new Dual(re, du);
    }

    public static Dual acsc(Dual a) {
        // acsc(a+bε) = acsc(a) - b/(a^2 sqrt(1-1/(a^2)))ε
        double re = acsc(a.a());
        double asq = a.a() * a.a();
        double du = -a.b() / (asq * sqrt(1 - 1 / asq));
        return new Dual(re, du);
    }

    public static Dual coth(Dual a) {
        // coth(a+bε) = coth(a) - b/sinh^2(a)ε
        double re = coth(a.a());
        double sh = sinh(a.a());
        double du = -a.b() / (sh * sh);
        return new Dual(re, du);
    }

    public static Dual sech(Dual a) {
        // sech(a+bε) = sech(a) - b*tanh(a)*sech(a)ε
        double re = sech(a.a());
        double th = tanh(a.a());
        double du = -a.b() * th * re;
        return new Dual(re, du);
    }

    public static Dual csch(Dual a) {
        // csch(a+bε) = csch(a) - b*coth(a)*csch(a)ε
        double re = csch(a.a());
        double co = coth(a.a());
        double du = -a.b() * co * re;
        return new Dual(re, du);
    }

    public static Dual asinh(Dual a) {
        // asinh(a+bε) = asinh(a) + b/sqrt(1+a^2)ε
        double re = asinh(a.a());
        double du = a.b() / sqrt(1 + a.a() * a.a());
        return new Dual(re, du);
    }

    public static Dual acosh(Dual a) {
        // acosh(a+bε) = acosh(a) + b/(sqrt(a-1)sqrt(a+1))ε
        double re = acosh(a.a());
        double du = a.b() / (sqrt(a.a() - 1) * sqrt(a.a() + 1));
        return new Dual(re, du);
    }

    public static Dual atanh(Dual a) {
        // atanh(a+bε) = atanh(a) + b/(1-a^2)ε
        double re = atanh(a.a());
        double du = a.b() / (1 - a.a() * a.a());
        return new Dual(re, du);
    }

    public static Dual asech(Dual a) {
        // asech(a+bε) = asech(a) - b/(sqrt(1/a - 1) sqrt(1/a + 1) a^2)ε
        double re = asech(a.a());
        double du = -a.b() / (sqrt(1 / a.a() - 1) * sqrt(1 / a.a() + 1) * a.a() * a.a());
        return new Dual(re, du);
    }

    public static Dual acsch(Dual a) {
        // acsch(a+bε) = acsch(a) - b/(sqrt(1+1/a^2) a^2)ε
        double re = acsch(a.a());
        double du = -a.b() / (sqrt(1 + 1 / (a.a() * a.a())) * a.a() * a.a());
        return new Dual(re, du);
    }

    public static Dual cbrt(Dual a) {
        // cbrt(a+bε) = cbrt(a) + b/(3 cbrt(a)^2)ε
        double re = cbrt(a.a());
        double du = a.b() / (3 * re * re);
        return new Dual(re, du);
    }

    /**
     * Compute the beta function of a dual and real number.
     * The derivative is taken as d/da beta(a,c).
     * @param a
     * @param c
     * @return
     */
    public static Dual beta(Dual a, double c) {
        // beta(a+bε,c) = beta(a,c) + b(digamma(a) - digamma(a+c))*beta(a,c)ε
        double re = beta(a.a(), c);
        double du = a.b() * (digamma(a.a()) - digamma(a.a() + c)) * re;
        return new Dual(re, du);
    }

    /**
     * Compute the beta function of a real and dual number.
     * The derivative is taken as d/dc beta(a,c).
     * @param a
     * @param c
     * @return
     */
    public static Dual beta(double a, Dual c) {
        // beta(a,c+bε) = beta(a,c) + b(digamma(c) - digamma(a+c))*beta(a,c)ε
        double re = beta(a, c.a());
        double du = c.b() * (digamma(c.a()) - digamma(a + c.a())) * re;
        return new Dual(re, du);
    }

    public static Dual erfc(Dual a) {
        // erfc(a+bε) = erfc(a) - 2 b exp(-a^2) / sqrt(π)ε
        double re = erfc(a.a());
        double du = -2 * a.b() * exp(-a.a() * a.a()) / 1.7724538509055160272981674833411451827975494561223871282138077898;
        return new Dual(re, du);
    }

    public static Dual erfi(Dual a) {
        // erfi(a+bε) = erfi(a) + 2 b exp(a^2) / sqrt(π)ε
        double re = erfi(a.a());
        double du = 2 * a.b() * exp(a.a() * a.a()) / 1.7724538509055160272981674833411451827975494561223871282138077898;
        return new Dual(re, du);
    }

    public static Dual fresnelC(Dual a) {
        // fresnelC(a+bε) = fresnelC(a) + b cos(π a^2/2)ε
        double re = fresnelC(a.a());
        double du = a.b() * cos(a.a() * a.a() * PI_2);
        return new Dual(re, du);
    }

    public static Dual fresnelS(Dual a) {
        // fresnelS(a+bε) = fresnelS(a) + b sin(π a^2/2)ε
        double re = fresnelS(a.a());
        double du = a.b() * sin(a.a() * a.a() * PI_2);
        return new Dual(re, du);
    }

    public static Dual li(Dual a) {
        // li(a+bε) = li(a) + b/(log(a))ε
        double re = li(a.a());
        double du = a.b() / log(a.a());
        return new Dual(re, du);
    }

    public static Dual sinc(Dual a) {
        // sinc(a+bε) = sinc(a) + b(a cos(a) - sin(a))ε/a^2
        double re = sinc(a.a());
        double du = a.b() * (a.a() * cos(a.a()) - sin(a.a())) / (a.a() * a.a());
        return new Dual(re, du);
    }

    public static Dual spence(Dual a) {
        // spence(a+bε) = spence(a) + b*(log(a)/(1-a))ε
        double re = spence(a.a());
        double du = a.b() * log(a.a()) / (1 - a.a());
        return new Dual(re, du);
    }

    public static Dual trigamma(Dual a) {
        // trigamma(a+bε) = trigamma(a) + polygamma(2,a)bε
        double re = trigamma(a.a());
        double du = polygamma(2, a.a()) * a.b();
        return new Dual(re, du);
    }

    public static Dual polygamma(double n, Dual a) {
        // polygamma(n, a+bε) = polygamma(n, a) + polygamma(n+1,a)bε
        double re = polygamma(n, a.a());
        double du = polygamma(n + 1, a.a()) * a.b();
        return new Dual(re, du);
    }

    public static Dual polylog(int n, Dual a) {
        // polylog(n, a+bε) = polylog(n, a) + polylog(n-1,a)bε/a
        double re = polylog(n, a.a());
        double du = polylog(n - 1, a.a()) * a.b() / a.a();
        return new Dual(re, du);
    }

    public static Dual lambertW0(Dual a) {
        // W(a+bε) = W(a) + bεW(a)/(a(1+W(a)))
        double re = lambertW0(a.a());
        double du = a.b() * re / (a.a() * (1 + lambertW0(re)));
        return new Dual(re, du);
    }

    public static Dual lambertWm1(Dual a) {
        // W(a+bε) = W(a) + bεW(a)/(a(1+W(a)))
        double re = lambertWm1(a.a());
        double du = a.b() * re / (a.a() * (1 + lambertWm1(re)));
        return new Dual(re, du);
    }

    public static Dual hurwitzZeta(double s, Dual a) {
        // hurwitz(s, a+bε) = hurwitz(s, a) -s hurwitz(s+1,a)bε
        // Dunno about the d/da.
        double re = hurwitzZeta(s, a.a());
        double du = -s * hurwitzZeta(s + 1, a.a()) * a.b();
        return new Dual(re, du);
    }

    /**
     * Add a complex number and a real number.
     *
     * @param a
     * @param b
     * @return a + b
     */
    public static Complex add(Complex a, double b) {
        return new Complex(a.getReal() + b, a.getImaginary());
    }

    /**
     * Add a complex number and a real number.
     *
     * @param a
     * @param b
     * @return a + b
     */
    public static Complex add(double a, Complex b) {
        return new Complex(b.getReal() + a, b.getImaginary());
    }

    /**
     * Subtract two complex numbers from each other.
     *
     * @param a
     * @param b
     * @return a - b
     */
    public static Complex sub(Complex a, Complex b) {
        return new Complex(a.getReal() - b.getReal(), a.getImaginary() - b.getImaginary());
    }

    /**
     * Subtract a complex number and a real number.
     *
     * @param a
     * @param b
     * @return a - b
     */
    public static Complex sub(Complex a, double b) {
        return new Complex(a.getReal() - b, a.getImaginary());
    }

    /**
     * Subtract a complex number and a real number.
     *
     * @param a
     * @param b
     * @return a - b
     */
    public static Complex sub(double a, Complex b) {
        return new Complex(a - b.getReal(), -b.getImaginary());
    }

    /**
     * Multiply two complex numbers.
     *
     * @param a
     * @param b
     * @return a * b
     */
    public static Complex mul(Complex a, Complex b) {
        return new Complex(a.getReal() * b.getReal() - a.getImaginary() * b.getImaginary(), a.getReal() * b.getImaginary() + a.getImaginary() * b.getReal());
    }

    /**
     * Multiply a complex number and a real number.
     *
     * @param a
     * @param b
     * @return a * b
     */
    public static Complex mul(Complex a, double b) {
        return new Complex(a.getReal() * b, a.getImaginary() * b);
    }

    /**
     * Multiply a complex number and a real number.
     *
     * @param a
     * @param b
     * @return a * b
     */
    public static Complex mul(double a, Complex b) {
        return new Complex(b.getReal() * a, b.getImaginary() * a);
    }

    /**
     * Divide two complex numbers.
     *
     * @param a
     * @param b
     * @return a / b
     */
    public static Complex div(Complex a, Complex b) {
        double d = b.getReal() * b.getReal() + b.getImaginary() * b.getImaginary();
        return new Complex((a.getReal() * b.getReal() + a.getImaginary() * b.getImaginary()) / d, (a.getImaginary() * b.getReal() - a.getReal() * b.getImaginary()) / d);
    }

    /**
     * Divide a complex number and a real number.
     *
     * @param a
     * @param b
     * @return a / b
     */
    public static Complex div(Complex a, double b) {
        return new Complex(a.getReal() / b, a.getImaginary() / b);
    }

    /**
     * Divide a complex number and a real number.
     *
     * @param a
     * @param b
     * @return a / b
     */
    public static Complex div(double a, Complex b) {
        double d = b.getReal() * b.getReal() + b.getImaginary() * b.getImaginary();
        return new Complex(a * b.getReal() / d, -a * b.getImaginary() / d);
    }

    /**
     * Compute the complex conjugate of a complex number.
     *
     * @param a
     * @return conj(a)
     */
    public static Complex conj(Complex a) {
        return new Complex(a.getReal(), -a.getImaginary());
    }

    /**
     * Compute the absolute value of a complex number.
     *
     * @param a
     * @return |a|
     */
    public static double abs(Complex a) {
        return Math.sqrt(a.getReal() * a.getReal() + a.getImaginary() * a.getImaginary());
    }

    /**
     * Compute the square root of a complex number.
     *
     * @param x
     * @return sqrt(x)
     */
    public static Complex sqrt(Complex x) {
        if (x.getImaginary() == 0)
            return new Complex(Math.sqrt(x.getReal()), 0);
        else {
            double r = abs(x);
            double t = Math.sqrt(0.5 * (r + x.getReal()));
            double u = Math.sqrt(0.5 * (r - x.getReal()));
            if (x.getImaginary() > 0)
                return new Complex(t, u);
            else
                return new Complex(t, -u);
        }
    }

    /**
     * Compute the value of the exponential function of a complex number.
     *
     * @param x
     * @return exp(x)
     */
    public static Complex exp(Complex x) {
        double r = Math.exp(x.getReal());
        return new Complex(r * Math.cos(x.getImaginary()), r * Math.sin(x.getImaginary()));
    }

    /**
     * Compute the natural logarithm of a complex number.
     *
     * @param x
     * @return ln(x)
     */
    public static Complex log(Complex x) {
        // Log(z) = ln(r) + i*theta = ln(|z|) + i*arg(z) = ln(|z|) + i*atan2(im(z), re(z))
        // There's some ridiculous nonsense involving the sign of zero here.
        // Pretend that it doesn't exist...
        return new Complex(Math.log(abs(x)), Math.atan2(x.getImaginary() == -0.0 ? 0.0 : x.getImaginary(), x.getReal() == -0.0 ? 0.0 : x.getReal()));
    }

    /**
     * Return the argument of a complex number.
     *
     * @param z
     * @return Arg(z)
     */
    public static double arg(Complex z) {
        return Math.atan2(z.getImaginary(), z.getReal());
    }

    /**
     * Compute the value of the cis function of a number.
     *
     * @param x
     * @return cos(x) + i * sin(x)
     */
    public static Complex cis(double x) {
        return new Complex(Math.cos(x), Math.sin(x));
    }

    /**
     * Compare two numbers for equality.
     *
     * @param a
     * @param b
     * @return true if a == b, false otherwise
     */
    public static boolean eq(Complex a, Complex b) {
        return a.getReal() == b.getReal() && a.getImaginary() == b.getImaginary();
    }

    /**
     * Compare two numbers for equality.
     *
     * @param a
     * @param b
     * @return true if a == b, false otherwise
     */
    public static boolean eq(double a, Complex b) {
        return b.getReal() == a && b.getImaginary() == 0;
    }

    /**
     * Compare two numbers for equality.
     *
     * @param a
     * @param b
     * @return true if a == b, false otherwise
     */
    public static boolean eq(Complex a, double b) {
        return a.getReal() == b && a.getImaginary() == 0;
    }

    /**
     * Compare two numbers for inequality.
     *
     * @param a
     * @param b
     * @return true if a != b, false otherwise
     */
    public static boolean ne(Complex a, double b) {
        return b != a.getReal() || a.getImaginary() != 0;
    }

    /**
     * Compare two numbers for inequality.
     *
     * @param a
     * @param b
     * @return true if a != b, false otherwise
     */
    public static boolean ne(double a, Complex b) {
        return a != b.getReal() || b.getImaginary() != 0;
    }

    /**
     * Compare two numbers for inequality.
     *
     * @param a
     * @param b
     * @return true if a != b, false otherwise
     */
    public static boolean ne(Complex a, Complex b) {
        return a.getReal() != b.getReal() || a.getImaginary() != b.getImaginary();
    }

    /**
     * Compare two complex numbers for equality.
     *
     * @param a
     * @param b
     * @param tol
     * @return true if a ~= b, false otherwise
     */
    public static boolean eq(Complex a, Complex b, double tol) {
        return Math.abs(a.getReal() - b.getReal()) < tol && Math.abs(a.getImaginary() - b.getImaginary()) < tol;
    }

    /**
     * Compute the n-th root of a complex number.
     *
     * @param x
     * @param deg
     * @return x^(1/deg)
     */
    public static double root(double x, int deg) {
        if (deg < 0)
            return Double.NaN;
        return Math.pow(x, 1.0 / deg);
    }

    /**
     * Find all n-th roots of a complex number.
     *
     * @param x
     * @param n
     * @return x^(1/n)
     */
    public static Complex[] root(Complex x, int n) {
        if (n < 0)
            return null;
        double magnitude = abs(x);
        double phase = Math.atan2(x.getImaginary(), x.getReal());
        double nthRootOfMagnitude = root(magnitude, n);
        Complex[] roots = new Complex[n];
        for (int i = 0; i < n; i++) {
            double theta = (phase + 2 * Math.PI * i) / n;
            roots[i] = new Complex(nthRootOfMagnitude * Math.cos(theta), nthRootOfMagnitude * Math.sin(theta));
        }
        return roots;
    }

    /**
     * Compute the principal cube root of a complex number.
     *
     * @param x
     * @return cbrt(x)
     */
    public static Complex cbrt(Complex x) {
        double magnitude = abs(x);
        double phase = Math.atan2(x.getImaginary(), x.getReal());
        double nthRootOfMagnitude = root(magnitude, 3);
        double theta = (phase + 2 * Math.PI) / 3;
        return new Complex(nthRootOfMagnitude * Math.cos(theta), nthRootOfMagnitude * Math.sin(theta));
    }

    /**
     * Compute the sine of a complex number.
     *
     * @param x
     * @return sin(x)
     */
    public static Complex sin(Complex x) {
        // sin(a+bi) = sin(a)cosh(b) + i*cos(a)sinh(b)
        return new Complex(Math.sin(x.getReal()) * Math.cosh(x.getImaginary()), Math.cos(x.getReal()) * Math.sinh(x.getImaginary()));
    }

    /**
     * Compute the cosine of a complex number.
     *
     * @param x
     * @return cos(x)
     */
    public static Complex cos(Complex x) {
        // cos(a+bi) = cos(a)cosh(b) - i*sin(a)sinh(b)
        return new Complex(Math.cos(x.getReal()) * Math.cosh(x.getImaginary()), -Math.sin(x.getReal()) * Math.sinh(x.getImaginary()));
    }

    /**
     * Compute the tangent of a complex number.
     *
     * @param x
     * @return tan(x)
     */
    public static Complex tan(Complex x) {
        // tan(a+bi) = sin(2a)+i*sinh(2b) / cos(2a)+cosh(2b)
        double d = Math.cos(2 * x.getReal()) + Math.cosh(2 * x.getImaginary());
        return new Complex(Math.sin(2 * x.getReal()) / d, Math.sinh(2 * x.getImaginary()) / d);
    }

    /**
     * Compute the cotangent of a complex number.
     *
     * @param x
     * @return cot(x)
     */
    public static Complex cot(Complex x) {
        return div(1, tan(x));
    }

    /**
     * Compute the secant of a complex number.
     *
     * @param x
     * @return
     */
    public static Complex sec(Complex x) {
        return div(1, cos(x));
    }

    /**
     * Compute the cosecant of a complex number.
     *
     * @param x
     * @return
     */
    public static Complex csc(Complex x) {
        return div(1, sin(x));
    }

    /**
     * Compute a to the power of b, where a and b are both
     * complex numbers.
     *
     * @param a
     * @param b
     * @return a^b
     */
    public static Complex pow(Complex a, Complex b) {
        // a^b = exp(b*ln(a))
        return exp(mul(b, log(a)));
    }

    /**
     * Compute a to the power of b, where a is a double precision
     * number and b is a complex number.
     *
     * @param a
     * @param b
     * @return a^b
     */
    public static Complex pow(double a, Complex b) {
        // a^b = exp(b*ln(a))
        return exp(mul(b, log(a)));
    }

    /**
     * Compute a to the power of b, where a is a complex
     *
     * @param a
     * @param b
     * @return a^b
     */
    public static Complex pow(Complex a, long b) {
        // a^b = exp(b*ln(a))
        return exp(mul(b, log(a)));
    }

    /**
     * Compute the hyperbolic sine of a complex number.
     *
     * @param a
     * @return sinh(a)
     */
    public static Complex sinh(Complex a) {
        // sinh(a+bi) = sinh(a)cos(b) + i*cosh(a)sin(b)
        return new Complex(Math.sinh(a.getReal()) * Math.cos(a.getImaginary()), Math.cosh(a.getReal()) * Math.sin(a.getImaginary()));
    }

    /**
     * Compute the hyperbolic cosine of a complex number.
     *
     * @param a
     * @return cosh(a)
     */
    public static Complex cosh(Complex a) {
        // cosh(a+bi) = cosh(a)cos(b) + i*sinh(a)sin(b)
        return new Complex(Math.cosh(a.getReal()) * Math.cos(a.getImaginary()), Math.sinh(a.getReal()) * Math.sin(a.getImaginary()));
    }

    /**
     * Compute the hyperbolic tangent of a complex number.
     *
     * @param a
     * @return tanh(a)
     */
    public static Complex tanh(Complex a) {
        // tanh(a+bi) = sinh(2a)+i*sin(2b) / cosh(2a)+cos(2b)
        double d = Math.cosh(2 * a.getReal()) + Math.cos(2 * a.getImaginary());
        return new Complex(Math.sinh(2 * a.getReal()) / d, Math.sin(2 * a.getImaginary()) / d);
    }

    /**
     * Compute the hyperbolic cotangent of a complex number.
     *
     * @param a
     * @return coth(a)
     */
    public static Complex coth(Complex a) {
        return div(1, tanh(a));
    }

    /**
     * Compute the hyperbolic secant of a complex number.
     *
     * @param a
     * @return sech(a)
     */
    public static Complex sech(Complex a) {
        return div(1, cosh(a));
    }

    /**
     * Compute the hyperbolic cosecant of a complex number.
     *
     * @param a
     * @return csch(a)
     */
    public static Complex csch(Complex a) {
        return div(1, sinh(a));
    }

    /**
     * Compute the arcus sine of a complex number.
     *
     * @param a
     * @return asin(a)
     */
    public static Complex asin(Complex a) {
        // asin(z)=1/i Ln(iz+sqrt(1-z^2))
        Complex y = sqrt(sub(1, mul(a, a)));
        Complex lnt = log(add(mul(I, a), y));
        return div(lnt, I);
    }

    /**
     * Compute the arcus cosine of a complex number.
     *
     * @param a
     * @return acos(a)
     */
    public static Complex acos(Complex a) {
        // acos(z)=1/i Ln(z+sqrt(z^2-1))
        Complex y = sqrt(sub(mul(a, a), 1));
        Complex lnt = log(add(a, y));
        return div(lnt, I);
    }

    /**
     * Compute the arcus tangent of a complex number.
     *
     * @param a
     * @return atan(a)
     */
    public static Complex atan(Complex a) {
        // atan(z)=1/2 i Ln((1-iz)/(1+iz))
        // atan(z) = (I/2)*(log(1-I*z) - log(1+I*z))
        Complex ia = mul(I, a);
        Complex y = log(sub(1, ia));
        Complex z = log(add(1, ia));
        Complex lnt = sub(y, z);
        return mul(lnt, mul(0.5, I));
    }

    /**
     * Compute the arcus cotangent of a complex number.
     *
     * @param a
     * @return acot(a)
     */
    public static Complex acot(Complex a) {
        // acot(z)=1/2i Ln((z+i)/(z-i))
        Complex y = add(a, I);
        Complex z = sub(a, I);
        Complex lnt = log(div(y, z));
        return mul(lnt, mul(0.5, I));
    }

    /**
     * Compute the arcus cosecant of a complex number.
     *
     * @param a
     * @return acsc(a)
     */
    public static Complex acsc(Complex a) {
        // asec(a) = 1/i * ln((i+sqrt(a*a-1))/a)
        Complex y = add(I, sqrt(sub(mul(a, a), 1)));
        Complex lnt = log(div(y, a));
        return div(lnt, I);
    }

    /**
     * Compute the arcus secant of a complex number.
     *
     * @param a
     * @return asec(a)
     */
    public static Complex asec(Complex a) {
        // asec(a) = 1/i * ln((1+sqrt(1-z*z))/z)
        Complex y = add(1, sqrt(sub(1, mul(a, a))));
        Complex lnt = log(div(y, a));
        return div(lnt, I);
    }

    /**
     * Compute the hyperbolic arcsine of a complex number.
     *
     * @param a
     * @return asinh(a)
     */
    public static Complex asinh(Complex a) {
        // asinh(z) = ln(z+sqrt(z^2+1))
        Complex y = sqrt(add(mul(a, a), 1));
        return log(add(a, y));
    }

    /**
     * Compute the hyperbolic arccosine of a complex number.
     *
     * @param a
     * @return acosh(a)
     */
    public static Complex acosh(Complex a) {
        // acosh(z) = sqrt(z-1)/sqrt(1-z) acos(z)
        Complex y = sqrt(sub(a, 1));
        Complex z = sqrt(sub(1, a));
        return mul(acos(a), div(y, z));
    }

    /**
     * Compute the hyperbolic arctangent of a complex number.
     *
     * @param a
     * @return atanh(a)
     */
    public static Complex atanh(Complex a) {
        // atanh(z) = 1/2 ln((1+z)/(1-z))
        Complex y = add(1, a);
        Complex z = sub(1, a);
        return div(log(div(y, z)), new Complex(0, 2));
    }

    /**
     * Compute the hyperbolic arccotangent of a complex number.
     *
     * @param a
     * @return acoth(a)
     */
    public static Complex acoth(Complex a) {
        // atanh(z) = 1/2 ln((z+1)/(z-1))
        Complex y = add(a, 1);
        Complex z = sub(a, 1);
        return div(log(div(y, z)), new Complex(0, 2));
    }

    /**
     * Compute the hyperbolic arcsecant of a complex number.
     *
     * @param a
     * @return asech(a)
     */
    public static Complex asech(Complex a) {
        // asech(z) = ln((1+sqrt(1-z*z))/z)
        Complex y = add(1, sqrt(sub(1, mul(a, a))));
        return log(div(y, a));
    }

    /**
     * Negate a complex number.
     *
     * @param x
     * @return -x
     */
    public static Complex negate(Complex x) {
        return new Complex(-x.getReal(), -x.getImaginary());
    }

    /**
     * Negate a real number.
     *
     * @param x
     * @return -x
     */
    public static double negate(double x) {
        return -x;
    }

    /**
     * Compute the hyperbolic arccosecant of a complex number.
     *
     * @param a
     * @return acsch(a)
     */
    public static Complex acsch(Complex a) {
        // asec(a) = ln((1+sqrt(a*a+1))/a)
        Complex y = add(I, sqrt(add(mul(a, a), 1)));
        return log(div(y, a));
    }

    /**
     * Compute the gamma function of a complex number.
     *
     * @param a
     * @return gamma(a)
     */
    public static Complex gamma(Complex a) {
        return Gamma.gamma(a);
    }

    /**
     * Compute the beta function of two complex numbers.
     * Defined as beta(a, b) = gamma(a) * gamma(b) / gamma(a + b)
     *
     * @param a
     * @param b
     * @return beta(a, b)
     * @see #gamma(Complex)
     */
    public static Complex beta(Complex a, Complex b) {
        return div(mul(gamma(a), gamma(b)), gamma(add(a, b)));
    }

    /**
     * Compute the value of the Airy Ai function at the specified point.
     *
     * @param x
     * @return Ai(x)
     */
    public static Complex airyAi(Complex x) {
        return Airy.airy(x)[0];
    }

    /**
     * Compute the value of the Airy Ai function's first derivative at the specified point.
     *
     * @param x
     * @return Ai'(x)
     */
    public static Complex airyAip(Complex x) {
        return Airy.airy(x)[1];
    }

    /**
     * Compute the value of the Airy Bi function at the specified point.
     *
     * @param x
     * @return Bi(x)
     */
    public static Complex airyBi(Complex x) {
        return Airy.airy(x)[2];
    }

    /**
     * Compute the value of the Airy Bi function's first derivative at the specified point.
     *
     * @param x
     * @return Bi'(x)
     */
    public static Complex airyBip(Complex x) {
        return Airy.airy(x)[3];
    }

    /**
     * Compute the value of the Airy Ai, Ai', Bi and Bi' functions at the specified point.
     *
     * @param x
     * @return a double array of length 4 containing Ai, Ai', Bi and Bi' in that order.
     */
    public static Complex[] airy(Complex x) {
        return Airy.airy(x);
    }

    /**
     * Compute the value of the exponential integral E1 at the specified point.
     * For positive values of real x, E1 and Ei relate as -E1(x) = Ei(-x).
     *
     * @param x
     * @return E1(x)
     */
    public static Complex e1(Complex x) {
        return Ei.e1(x);
    }

    /**
     * Compute the exponential integral using a formula that was revealed to me in a dream:
     * <p> Ei(z) = -E1(-z) - log(-z) + 0.5 (log(z) - log(1/z))
     *
     * @param x
     * @return Ei(x)
     */
    public static Complex Ei(Complex x) {
        return add(sub(negate(e1(negate(x))), log(negate(x))), mul(0.5, sub(log(x), log(div(1, x)))));
    }

    /**
     * Compute the value of the complementary exponential integral Ein at the specified point.
     * For all complex z, Ein and Ei relate as Ein(z) = E1(z) + EulerGamma + ln z.
     *
     * @param x
     * @return Ein(x)
     * @see #e1(Complex)
     */
    public static Complex Ein(Complex x) {
        // E1(z) = Ein(z) - ln z - EulerGamma
        // ... => Ein(z) = E1(z) + ln z + EulerGamma
        return add(e1(x), add(log(x), EULER_GAMMA));
    }

    /**
     * <p>Principal branch of the logarithm of the gamma function.
     * Defined to be log(gamma(x)) for x &gt; 0.
     * Extended to the complex plane by analytic continuation.
     * The function has a single branch cut on the negative real axis.
     *
     * <p>It is not generally true that log gamma(x) = log(gamma(x)),
     * though the real parts tend to agree. The benefit of not defining
     * loggamma as is that the latter function has a complicated branch
     * cut structure whereas loggamma is analytic except for on the negative real axis.
     *
     * <p>On the real line, loggamma is related to logabsgamma.
     *
     * @param x
     * @return
     * @see #logabsgamma(double)
     */
    public static Complex loggamma(Complex x) {
        return Gamma.loggamma(x);
    }

    /**
     * Compute the logarithm of the beta function in the complex plane.
     * Defined using complex loggamma as logbeta(a, b) = loggamma(a) - loggamma(b) - loggamma(a + b).
     *
     * @param a
     * @param b
     * @return logbeta(a, b)
     * @see #loggamma(Complex)
     */
    public static Complex logbeta(Complex a, Complex b) {
        return sub(loggamma(a), add(loggamma(b), loggamma(add(a, b))));
    }

    /**
     * Compute the lower incomplete (non-regularised) gamma function of a complex number.
     *
     * @param s
     * @param z
     * @return ligamma(s, z)
     */
    public static Complex liGamma(Complex s, Complex z) {
        // z^s e^-z sum(k=0, inf, z^k/pochhammer(s, k+1))
        if (eq(s, Complex.ZERO))
            throw new ArithmeticException("s=0 pole.");
        Complex zs = pow(z, s);
        Complex ez = exp(negate(z));
        Complex sum = Complex.ZERO;
        int maxiter = 50;
        for (int k = 0; k < maxiter; k++) {
            Complex term = div(pow(z, k), pochhammer(s, k + 1));
            sum = add(sum, term);
            if (abs(term) <= Maja.EPSILON)
                break;
        }
        return mul(mul(zs, ez), sum);
    }

    /**
     * Compute the upper incomplete (non-regularised) gamma function of a complex number.
     *
     * @param s
     * @param z
     * @return uigamma(s, z)
     */
    public static Complex uiGamma(Complex s, Complex z) {
        // ligamma(s,x) + uigamma(s,x) = gamma(s)
        // ... => uigamma(s,x) = gamma(s) - ligamma(s,x)
        return sub(gamma(s), liGamma(s, z));
    }

    /**
     * Compute the complex error function.
     *
     * @param z
     * @return erf(z)
     */
    public static Complex erf(Complex z) {
        return Erf.cerf(z);
    }

    /**
     * Compute the complex complementary error function.
     *
     * @param z
     * @return erfc(z)
     */
    public static Complex erfc(Complex z) {
        return Erf.cerfc(z);
    }

    /**
     * Compute the complex imaginary error function.
     *
     * @param z
     * @return erfi(z) = -i erf(iz)
     */
    public static Complex erfi(Complex z) {
        return Erf.cerfi(z);
    }

    /**
     * Compute the value of the complex Dawson function (D+) at z.
     *
     * @param z
     * @return D+(z)
     */
    public static Complex dawsonPlus(Complex z) {
        return Erf.cdawson(z);
    }

    /**
     * Compute the value of the complex Dawson function (D-) at z.
     *
     * @param z
     * @return D-(z)
     */
    public static Complex dawsonMinus(Complex z) {
        // sqrt(pi)/2 * exp(z*2) * erf(z)
        return mul(mul(0.8862269254527580136490837416705725913990, exp(add(z, z))), erf(z));
    }

    /**
     * Compute the Fresnel S integral on the complex plane.
     *
     * @param z
     * @return FresnelS(z)
     */
    public static Complex fresnelS(Complex z) {
        // FresnelS[z] == ((1 + I)/4) (Erf[((1 + I)/2) Sqrt[Pi] z] - I Erf[((1 - I)/2) Sqrt[Pi] z])
        double sqp = 1.7724538509055160272981674833411451827975494561223871282138077898;
        Complex t1 = div(add(1, I), 4);
        Complex t2 = erf(mul(div(add(1, I), 2), mul(sqp, z)));
        Complex t3 = mul(I, erf(mul(div(sub(1, I), 2), mul(sqp, z))));
        return mul(t1, sub(t2, t3));
    }

    /**
     * Compute the Fresnel C integral on the complex plane.
     *
     * @param z
     * @return FresnelC(z)
     */
    public static Complex fresnelC(Complex z) {
        // FresnelC[z] == ((1 - I)/4) (Erf[((1 + I)/2) Sqrt[Pi] z] + I Erf[((1 - I)/2) Sqrt[Pi] z])
        double sqp = 1.7724538509055160272981674833411451827975494561223871282138077898;
        Complex t1 = div(sub(1, I), 4);
        Complex t2 = erf(mul(div(add(1, I), 2), mul(sqp, z)));
        Complex t3 = mul(I, erf(mul(div(sub(1, I), 2), mul(sqp, z))));
        return mul(t1, add(t2, t3));
    }

    /**
     * Compute the complex digamma function.
     *
     * @param z
     * @return digamma(z)
     */
    public static Complex digamma(Complex z) {
        return Gamma.digamma(z);
    }

    /**
     * Compute the complex trigamma function.
     *
     * @param z
     * @return trigamma(z)
     */
    public static Complex trigamma(Complex z) {
        return Gamma.trigamma(z);
    }

    /**
     * Compute the logarithmic integral of x, defined as li(x) = int(1/log t, t=0..x).
     *
     * @param x
     * @return li(x)
     */
    public static double li(double x) {
        return Ei.expint(log(x));
    }

    /**
     * Compute the complex logarithmic integral of z, defined as li(z) = int(1/log t, t=0..z).
     *
     * @param z
     * @return li(z)
     */
    public static Complex li(Complex z) {
        return Ei(log(z));
    }

    /**
     * Compute the value of the complex sine integral Si(z) at z.
     *
     * @param z
     * @return Si(z)
     */
    public static Complex Si(Complex z) {
        return TrigonometricIntegral.Si(z);
    }

    /**
     * Compute the value of the complex sine integral si(z) at z.
     *
     * @param z
     * @return si(z)
     */
    public static Complex si(Complex z) {
        return TrigonometricIntegral.si(z);
    }

    /**
     * Compute the value of the complex cosine integral Ci(z) at z.
     *
     * @param z
     * @return Ci(z)
     */
    public static Complex Ci(Complex z) {
        return TrigonometricIntegral.Ci(z);
    }

    /**
     * Compute the value of the complex cosine integral Cin(z) at z.
     *
     * @param z
     * @return Cin(z)
     */
    public static Complex Cin(Complex z) {
        return TrigonometricIntegral.Cin(z);
    }

    /**
     * Compute the value of the complex hyperbolic sine integral Shi(z).
     *
     * @param z
     * @return Shi(z)
     */
    public static Complex Shi(Complex z) {
        return TrigonometricIntegral.Shi(z);
    }

    /**
     * Compute the value of the complex hyperbolic cosine integral Chi(z).
     *
     * @param z
     * @return Chi(z)
     */
    public static Complex Chi(Complex z) {
        return TrigonometricIntegral.Chi(z);
    }

    /**
     * Compute the value of the complex hyperbolic sine integral Shi(z)
     * and complex hyperbolic cosine integral Chi(z).
     *
     * @param z
     * @return { Shi(z), Chi(z) }
     */
    public static Complex[] ShiChi(Complex z) {
        return new Complex[]{Shi(z), Chi(z)};
    }

    /**
     * Compute the value of the complex generalised exponential integral E_n(z).
     * Uses A&amp;S 5.1.45 E_n(z) = z^(n-1) * uiGamma(1-n, z) to perform computation.
     *
     * @param n
     * @param z
     * @return E_n(z)
     */
    public static Complex en(Complex n, Complex z) {
        // z^(n-1) * uiGamma(1-n, z)
        return mul(pow(z, sub(n, 1)), uiGamma(sub(1, n), z));
    }

    /**
     * Compute the Hurwitz zeta function of complex arguments.
     *
     * @param s
     * @param a
     * @return zeta(s, a)
     */
    public static Complex hurwitzZeta(Complex s, Complex a) {
        return Zeta.hurwitz_zeta(s, a);
    }

    /**
     * Compute the complex dilogarithm of z.
     *
     * @param z
     * @return dilog(z)
     */
    public static Complex dilog(Complex z) {
        return Spence.dilog(z);
    }

    /**
     * Compute the Spence function of z.
     *
     * @param z
     * @return spence(z)
     */
    public static Complex spence(Complex z) {
        return Spence.spence(z);
    }

    /**
     * Compute the polylogarithm of s and z.
     *
     * @param s
     * @param z
     * @return polylog(s, z)
     * @throws ArithmeticException if the amount of numerical algorithm iterations is exceeded
     */
    public static Complex polylog(Complex s, Complex z) {
        return Spence.polylog(s, z);
    }

    /**
     * Compute the Lerch transcendent of z, s, and a.
     * A few things to note:
     * <ul>
     *     <li>The implementation may overflow the stack for particularly large, negative values of Re(a).
     *     This is being worked on.</li>
     *     <li>The precision of the output may vary. Worst-case scenario revealed during non-extensive randomised
     *     trials is around 0.006% relative error</li>
     * </ul>
     *
     * @param z
     * @param s
     * @param a
     * @return lerch(z, s, a)
     */
    public static Complex lerchPhi(Complex z, Complex s, Complex a) {
        return Zeta.lerch_phi(z, s, a);
    }

    /**
     * Compute the specified branch of the complex Lambert W function of z.
     *
     * @param z
     * @param k
     * @return lambertw(z, k)
     */
    public static Complex lambertw(Complex z, long k) {
        return Lambert.lambertW(z, k);
    }

    /**
     * Trim insignificant real/imaginary parts (below machine epsilon), round up numbers where
     * the real/imaginary part is very close to an integer.
     *
     * @param z
     * @return iround(z)
     */
    public static Complex chop(Complex z) {
        if (Math.abs(z.getReal()) < EPSILON)
            z = new Complex(0, z.getImaginary());
        if (Math.abs(z.getImaginary()) < EPSILON)
            z = new Complex(z.getReal(), 0);
        if (Math.abs(z.getReal() - Math.round(z.getReal())) < EPSILON)
            z = new Complex(Math.round(z.getReal()), z.getImaginary());
        if (Math.abs(z.getImaginary() - Math.round(z.getImaginary())) < EPSILON)
            z = new Complex(z.getReal(), Math.round(z.getImaginary()));
        return z;
    }

    /**
     * Compute the surface area of a solid of revolution created by rotating the function
     * f(x) about the x-axis. The function f(x) must be continuous and differentiable
     * on the interval [a, b]. The area is computed using the Gauss-Legendre quadrature
     * as 2pi * integral of f(x) * sqrt(1 + df(x)^2) from a to b.
     *
     * @param f  the function to rotate
     * @param df the derivative of f
     * @param a  the lower bound of the interval
     * @param b  the upper bound of the interval
     * @return the area of the solid of revolution
     */
    public static double solidArea(Function<Double, Double> f, Function<Double, Double> df, double a, double b) {
        return TWO_PI * integrateGaussLegendreReal(x -> f.apply(x) * sqrt(1 + pow(df.apply(x), 2)), a, b, 10);
    }

    /**
     * Compute the volume of a solid of revolution created by rotating the function f(x)
     * about the x-axis. The function f(x) must be continuous and differentiable on the
     * interval [a, b]. The volume is computed using the Gauss-Legendre quadrature as
     * pi * integral of f(x)^2 from a to b (the disk method).
     *
     * @param f the function to rotate
     * @param a the lower bound of the interval
     * @param b the upper bound of the interval
     * @return the volume of the solid of revolution
     */
    public static double solidVolume(Function<Double, Double> f, double a, double b) {
        return PI * integrateGaussLegendreReal(x -> pow(f.apply(x), 2), a, b, 10);
    }

    /**
     * Integrate a C -&gt; C function using the Tanh-Sinh quadrature. Performs
     * integration through a straight line contour from a to b.
     *
     * @param f   the function to integrate
     * @param a   the lower bound of the interval
     * @param b   the upper bound of the interval
     * @param n   the quadrature degree
     * @param eps the desired accuracy
     * @return the integral of f from a to b
     */
    public static Complex[] integrateTanhSinhComplex(Function<Complex, Complex> f, Complex a, Complex b, int n, double eps) {
        return Integrator.finiteTanhSinhCC(f, a, b, n, eps);
    }

    /**
     * Integrate a C -&gt; C function using the Gauss-Legendre quadrature. Performs
     * integration through a straight line contour from a to b.
     *
     * @param f the function to integrate
     * @param a the lower bound of the interval
     * @param b the upper bound of the interval
     * @param n the quadrature degree
     * @return the integral of f from a to b
     */
    public static Complex integrateGaussLegendreComplex(Function<Complex, Complex> f, Complex a, Complex b, int n) {
        return Integrator.gaussLegendreIntegrateCC(f, a, b, n);
    }

    /**
     * Compute the arc length of a curve defined by the function f(x) and its derivative
     * df(x) on the interval [a, b]. The arc length is computed using the Gauss-Legendre
     * quadrature as the integral of sqrt(1 + df(x)^2) from a to b.
     *
     * @param df the derivative of f
     * @param a  the lower bound of the interval
     * @param b  the upper bound of the interval
     * @return the arc length of the curve
     */
    public static double arcLength(Function<Double, Double> df, double a, double b) {
        return integrateGaussLegendreReal(x -> sqrt(1 + pow(df.apply(x), 2)), a, b, 10);
    }

    /**
     * Compute the Legendre F elliptic integral defined by DLMF §19.2.4
     *
     * @param phi
     * @param k
     * @return legendreF(phi, k)
     */
    public static double legendreF(double phi, double k) {
        return LegendreIntegral.legendreF(phi, k);
    }

    /**
     * Compute the Legendre E elliptic integral defined by DLMF §19.2.5
     *
     * @param phi
     * @param k
     * @return legendreE(phi, k)
     */
    public static double legendreE(double phi, double k) {
        return LegendreIntegral.legendreE(phi, k);
    }

    /**
     * Compute the Legendre D elliptic integral defined by DLMF §19.2.6
     *
     * @param phi
     * @param k
     * @return legendreD(phi, k)
     */
    public static double legendreD(double phi, double k) {
        return LegendreIntegral.legendreD(phi, k);
    }

    /**
     * Compute the Legendre Pi elliptic integral defined by DLMF §19.2.7
     *
     * @param phi
     * @param alpha
     * @param k
     * @return legendrePi(phi, alpha, k)
     */
    public static double legendrePi(double phi, double alpha, double k) {
        return LegendreIntegral.legendrePi(phi, alpha, k);
    }

    /**
     * Compute the Legendre F elliptic integral defined by DLMF §19.2.4
     *
     * @param phi
     * @param k
     * @return legendreF(phi, k)
     */
    public static Complex legendreF(Complex phi, Complex k) {
        return LegendreIntegral.legendreF(phi, k);
    }

    /**
     * Compute the Legendre E elliptic integral defined by DLMF §19.2.5
     *
     * @param phi
     * @param k
     * @return legendreE(phi, k)
     */
    public static Complex legendreE(Complex phi, Complex k) {
        return LegendreIntegral.legendreE(phi, k);
    }

    /**
     * Compute the Legendre D elliptic integral defined by DLMF §19.2.6
     *
     * @param phi
     * @param k
     * @return legendreD(phi, k)
     */
    public static Complex legendreD(Complex phi, Complex k) {
        return LegendreIntegral.legendreD(phi, k);
    }

    /**
     * Compute the Legendre Pi elliptic integral defined by DLMF §19.2.7
     *
     * @param phi
     * @param alpha
     * @param k
     * @return legendrePi(phi, alpha, k)
     */
    public static Complex legendrePi(Complex phi, Complex alpha, Complex k) {
        return LegendreIntegral.legendrePi(phi, alpha, k);
    }

    /**
     * Compute the Landau function with specified most probable value and sigma value.
     *
     * @param x     the value to evaluate the function at
     * @param mpv   the most probable value
     * @param sigma width of the distribution
     * @param norm  whether to normalize the result
     * @return landau(x, mpv, sigma, norm)
     */
    public static double landau(double x, double mpv, double sigma, boolean norm) {
        return Landau.landau(x, mpv, sigma, norm);
    }

    /**
     * Compute the Landau distribution function.
     *
     * @param x the value to evaluate the function at
     * @return landau(x)
     */
    public static double landau(double x) {
        return Landau.landauI(x);
    }

    /**
     * Compute the regularized Gamma P function.
     *
     * @param a order
     * @param x argument
     * @return gammaP(a, x)
     * @throws ArithmeticException if arguments are outside of the domain or the iteration count is exceeded.
     */
    public static double gammaP(double a, double x) {
        return Gamma.regularizedGammaP(a, x);
    }

    /**
     * Compute the regularized Gamma Q function.
     *
     * @param a order
     * @param x argument
     * @return gammaQ(a, x)
     * @throws ArithmeticException if arguments are outside of the domain or the iteration count is exceeded.
     */
    public static double gammaQ(double a, double x) {
        return Gamma.regularizedGammaQ(a, x);
    }

    /**
     * Computes quantiles for standard normal distribution N(0, 1) at probability p.
     *
     * @param p probability between 0 and 1.
     * @return quantile value
     * @throws IllegalArgumentException if p is not between 0 and 1
     */
    public static double normQuantile(double p) {
        return Landau.normQuantile(p);
    }

    /**
     * Computes quantiles for chi-squared probability distribution at probability p .
     *
     * @param p  probability between 0 and 1.
     * @param df degrees of freedom
     * @return quantile value
     * @throws IllegalArgumentException if p is not between 0 and 1
     */
    public static double chiSquaredQuantile(double p, double df) {
        return Landau.chisquareQuantile(p, df);
    }

    /**
     * Integrate a dyadic function over a finite interval [a,b] x [c,d] using the
     * Gauss-Legendre quadrature. The number of intervals is given by N and the precision
     * of the final result greatly depends on this parameter.
     * The computation of an integral using the Gauss-Legendre quadrature involves caching the
     * coefficients required to perform the integration depending on the value of the N parameter.
     * This means that the first call to this method will be slower than subsequent calls with the
     * same value of N. The coefficients are internally cached inside a ConcurrentHashMap.
     *
     * @param f function to integrate
     * @param a lower bound
     * @param b upper bound
     * @param c lower bound
     * @param d upper bound
     * @param N number of intervals, N=6 tends to give a good approximation in most scenarios.
     *          N must be between 1 and 30.
     * @return integral of f over [a,b] x [c,d]
     * @throws IllegalArgumentException if N is not a positive integer
     */
    public static double integrateGaussLegendreReal(BiFunction<Double, Double, Double> f, double a, double b, double c, double d, int N) {
        return integrateGaussLegendreReal(x -> integrateGaussLegendreReal(y -> f.apply(x, y), c, d, N), a, b, N);
    }

    /**
     * Integrate a dyadic function over a finite interval [a,b] x [c,d] using the
     * Gauss-Legendre quadrature. The number of intervals is given by N and the precision
     * of the final result greatly depends on this parameter.
     * The computation of an integral using the Gauss-Legendre quadrature involves caching the
     * coefficients required to perform the integration depending on the value of the N parameter.
     * This means that the first call to this method will be slower than subsequent calls with the
     * same value of N. The coefficients are internally cached inside a ConcurrentHashMap.
     *
     * @param f function to integrate
     * @param a lower bound
     * @param b upper bound
     * @param c lower bound
     * @param d upper bound
     * @param N number of intervals, N=6 tends to give a good approximation in most scenarios.
     *          N must be between 1 and 30.
     * @return integral of f over [a,b] x [c,d]
     * @throws IllegalArgumentException if N is not a positive integer
     */
    public static Complex integrateGaussLegendreComplex(BiFunction<Complex, Complex, Complex> f, Complex a, Complex b, Complex c, Complex d, int N) {
        return integrateGaussLegendreComplex(x -> integrateGaussLegendreComplex(y -> f.apply(x, y), c, d, N), a, b, N);
    }

    /**
     * Primality test for an integer in range [0, 2^31-1].
     */
    public static boolean isPrime(int n) {
        return Prime.is_prime_1b(n);
    }

    /**
     * Primality test for an unsigned integer in range [0, 2^64-1].
     * Notice that some values of n that will be interpreted as negative when signed
     * may yield a truthy value, as the number is reinterpreted as unsigned.
     */
    public static boolean isPrime(long n) {
        return Prime.is_prime_2_64(n);
    }

    /**
     * Find the next prime after n.
     * @throws IllegalArgumentException if n is negative
     * @throws ArithmeticException if n is too large
     */
    public static int nextPrime(int n) {
        if (n < 0)
            throw new IllegalArgumentException("n must be positive.");
        // Only even prime.
        if (n == 2) return 3;
        // 2^31-1 is a Mersenne prime.
        if (n == 2147483647)
            throw new ArithmeticException("n too large.");
        // If n is even, add 1.
        if (n % 2 == 0) n++;
        while (!isPrime(n))
            n += 2;
        return n;
    }

    // 2^64-59 is the largest 64-bit prime.
    private static final Long MAX_P = Long.parseUnsignedLong("18446744073709551557");

    /**
     * Find the next prime after (unsigned) n.
     * @throws ArithmeticException if n is too large
     */
    public static long nextPrime(long n) {
        // Only even prime.
        if (n == 2) return 3;
        if (Long.compareUnsigned(n, MAX_P) >= 0)
            throw new ArithmeticException("n too large.");
        // If n is even, add 1.
        if (n % 2 == 0) n++;
        while (!isPrime(n))
            n += 2;
        return n;
    }
}
