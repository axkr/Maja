import static org.assertj.core.api.Assertions.assertThat;
import org.hipparchus.complex.Complex;
import org.junit.jupiter.api.Test;
import rocks.palaiologos.maja.Maja;

public class TestComplex {
  @Test
  public void testArithmetic() {
    assertThat(Maja.add(new Complex(1, 2), new Complex(3, 4))).isEqualTo(new Complex(4, 6));
    assertThat(Maja.add(new Complex(1, 2), 3)).isEqualTo(new Complex(4, 2));
    assertThat(Maja.add(3, new Complex(1, 2))).isEqualTo(new Complex(4, 2));

    assertThat(Maja.sub(new Complex(1, 2), new Complex(3, 4))).isEqualTo(new Complex(-2, -2));
    assertThat(Maja.sub(new Complex(1, 2), 3)).isEqualTo(new Complex(-2, 2));
    assertThat(Maja.sub(3, new Complex(1, 2))).isEqualTo(new Complex(2, -2));

    assertThat(Maja.mul(new Complex(1, 2), new Complex(3, 4))).isEqualTo(new Complex(-5, 10));
    assertThat(Maja.mul(new Complex(1, 2), 3)).isEqualTo(new Complex(3, 6));
    assertThat(Maja.mul(3, new Complex(1, 2))).isEqualTo(new Complex(3, 6));

    assertThat(Maja.div(new Complex(1, 2), new Complex(3, 4)))
        .isEqualTo(new Complex(11.0 / 25, 2.0 / 25));
    assertThat(Maja.div(new Complex(1, 2), 3)).isEqualTo(new Complex(1.0 / 3, 2.0 / 3));
    assertThat(Maja.div(3, new Complex(1, 2))).isEqualTo(new Complex(0.6, -1.2));

    assertThat(Complex.ZERO).isEqualTo(new Complex(0, 0));
    assertThat(new Complex(1.23, 4.56).toString()).isEqualTo("(1.23, 4.56)");
    assertThat(new Complex(1.23, -4.56).toString()).isEqualTo("(1.23, -4.56)");
  }

  @Test
  public void testConjugate() {
    assertThat(Maja.conj(new Complex(1, 2))).isEqualTo(new Complex(1, -2));
  }

  @Test
  public void testAbsSqrt() {
    assertThat(Maja.abs(new Complex(3, 4))).isEqualTo(5);
    assertThat(Maja.abs(new Complex(3, -4))).isEqualTo(5);
    assertThat(Maja.sqrt(new Complex(3, 4))).isEqualTo(new Complex(2, 1));
    assertThat(Maja.sqrt(new Complex(4, 0))).isEqualTo(new Complex(2, 0));
  }

  @Test
  public void testComplexExp() {
    assertThat(Maja.exp(new Complex(0, 0))).isEqualTo(new Complex(1, 0));
    assertThat(Maja.exp(new Complex(1, 0))).isEqualTo(new Complex(Math.E, 0));
    assertThat(Maja.exp(new Complex(0, 1))).isEqualTo(new Complex(Math.cos(1), Math.sin(1)));
    assertThat(Maja.exp(new Complex(1.2, 3.4)))
        .isEqualTo(new Complex(Math.cos(3.4) * Math.exp(1.2), Math.sin(3.4) * Math.exp(1.2)));
    assertThat(Maja.exp(new Complex(1.2, -3.4)))
        .isEqualTo(new Complex(Math.cos(-3.4) * Math.exp(1.2), Math.sin(-3.4) * Math.exp(1.2)));
  }

  @Test
  public void testComplexLog() {
    assertThat(Maja.log(new Complex(1, 0))).isEqualTo(new Complex(0, 0));
    assertThat(Maja.log(new Complex(0, 1))).isEqualTo(new Complex(0, Math.PI / 2));
    assertThat(Maja.log(new Complex(0, -1))).isEqualTo(new Complex(0, -Math.PI / 2));
    assertThat(Maja.log(new Complex(1, 1)))
        .isEqualTo(new Complex(0.3465735902799727, 0.7853981633974483));
    assertThat(Maja.log(new Complex(1, -1)))
        .isEqualTo(new Complex(0.3465735902799727, -0.7853981633974483));
    assertThat(Maja.log(new Complex(1, 2)))
        .isEqualTo(new Complex(0.8047189562170503, 1.1071487177940904));
    assertThat(Maja.log(new Complex(1, -2)))
        .isEqualTo(new Complex(0.8047189562170503, -1.1071487177940904));
  }

  @Test
  public void testCmp() {
    assertThat(Maja.eq(new Complex(1, 2), new Complex(3, 4))).isFalse();
    assertThat(Maja.eq(new Complex(1, 2), new Complex(1, 2))).isTrue();

    assertThat(Maja.ne(new Complex(1, 2), new Complex(3, 4))).isTrue();
    assertThat(Maja.ne(new Complex(1, 2), new Complex(1, 2))).isFalse();
  }

  @Test
  public void testGamma() {
    assertThat(Maja.liGamma(new Complex(2, 5), new Complex(3, 6)))
        .isEqualTo(new Complex(0.003951212582048221, -0.007307084423671122));
    assertThat(Maja.liGamma(new Complex(10, 10), new Complex(10, 10)))
        .isEqualTo(new Complex(711.1040308344138, -1881.5622609709533));
    assertThat(Maja.digamma(new Complex(4, 4)))
        .isEqualTo(new Complex(1.6703598173334107, 0.8505022091860445));

    Complex[] poly = new Complex[] {new Complex(3, 0), new Complex(-4, 0), new Complex(0, 0),
        new Complex(8, 0), new Complex(-10, 0)};
    boolean[] err = Maja.aberth(poly);
    assertThat(err[0]).isFalse();
    assertThat(poly[0]).isEqualTo(new Complex(-0.7286839684871557));
    assertThat(err[1]).isFalse();
    assertThat(poly[1]).isEqualTo(new Complex(0.3758484858246309, -0.6233832239396144));
    assertThat(err[2]).isFalse();
    assertThat(poly[2]).isEqualTo(new Complex(0.3758484858246309, 0.6233832239396144));
    assertThat(err[3]).isFalse();
    assertThat(poly[3]).isEqualTo(new Complex(0.776986996837894));
    assertThat(err[4]).isTrue();
  }

  @Test
  public void testFresnel() {
    assertThat(Maja.fresnelC(new Complex(2, 5)))
        .isEqualTo(new Complex(-4.768856847987491E11, 1.221373671098564E12));
    assertThat(Maja.fresnelS(new Complex(2, 5)))
        .isEqualTo(new Complex(-1.221373671098064E12, -4.768856847992491E11));
  }

  @Test
  public void testLi() {
    assertThat(Maja.li(new Complex(1, 2)))
        .isEqualTo(new Complex(1.387678742022938, 2.508754698859233));
    assertThat(Maja.e1(new Complex(-10, 10)))
        .isEqualTo(new Complex(1576.1504265768526, 433.77763904754283));
  }

  @Test
  public void testZeta() {
    assertThat(Maja.zeta(new Complex(-2, 2)))
        .isEqualTo(new Complex(0.08638207303300506, 0.02053604281694063));
    assertThat(Maja.zeta(new Complex(20, 2)))
        .isEqualTo(new Complex(1.0000001747892098, -9.377209774992701E-7));
    assertThat(Maja.abs(Maja.sub(Maja.hurwitzZeta(new Complex(-2, 2), new Complex(1, 2)),
        new Complex(36.4631498760176, -3.5135475777949)))).isLessThan(3.8e-5);
    assertThat(Maja.abs(Maja.sub(Maja.hurwitzZeta(new Complex(5, 2), new Complex(1, 2)),
        new Complex(0.13723775153235927778, -0.11147445278404573)))).isLessThan(1e-13);
    assertThat(
        Maja.abs(Maja.sub(Maja.lerchPhi(new Complex(-2, 2), new Complex(1, 2), new Complex(2, 2)),
            new Complex(-0.539304744092355, -0.611340474677990)))).isLessThan(3.411e-7);
    assertThat(
        Maja.abs(Maja.sub(Maja.lerchPhi(new Complex(-2, 2), new Complex(1, 2), new Complex(-2, 2)),
            new Complex(-15.1536239656553, 43.4891357967102)))).isLessThan(2.19e-5);
    assertThat(
        Maja.abs(Maja.sub(Maja.lerchPhi(new Complex(0, -8), new Complex(1, -1), new Complex(2, 1)),
            new Complex(-0.0210382763183311, -0.0598747533099067)))).isLessThan(1.74e-9);
  }

  @Test
  public void testComplexTrig() {
    assertThat(Maja.abs(Maja.sub(Maja.sin(new Complex(-2.612, 3.162)),
        new Complex(-5.976332784021184, -10.172956161697219)))).isLessThan(1e-5);
  }
}
