package rocks.palaiologos.maja.expression;

import org.antlr.v4.runtime.tree.AbstractParseTreeVisitor;
import org.antlr.v4.runtime.tree.ParseTree;
import rocks.palaiologos.maja.Complex;
import rocks.palaiologos.maja.Maja;
import rocks.palaiologos.maja.matrix.ComplexMatrix;
import rocks.palaiologos.maja.matrix.DoubleMatrix;
import rocks.palaiologos.maja.matrix.Matrix;

import java.util.List;
import java.util.stream.Collectors;

public class DefaultExpressionVisitor extends AbstractParseTreeVisitor<Object> implements ExpressionVisitor<Object> {
    private Environment env;

    private Environment getEnv() {
        return env;
    }

    public DefaultExpressionVisitor(Environment env) {
        this.env = env;

        this.env.set("i", Maja.I);
        this.env.set("pi", Maja.PI);
        this.env.set("e", Maja.E);
        this.env.set("glaisher", Maja.GLAISHER_CONSTANT);
        this.env.set("catalan", Maja.CATALAN_CONSTANT);
        this.env.set("khinchin", Maja.KHINCHIN_CONSTANT);
        this.env.set("apery", Maja.APERY_CONSTANT);
        this.env.set("golden", Maja.GOLDEN_RATIO);
        this.env.set("euler_gamma", Maja.EULER_GAMMA);
        this.env.set("feigenbaum", Maja.FEIGENBAUM_CONSTANT);
        this.env.set("epsilon", Maja.EPSILON);
        this.env.set("pi2", Maja.PI_2);
        this.env.set("pi4", Maja.PI_4);
        this.env.set("tau", Maja.TWO_PI);
        this.env.set("invpi", Maja.ONE_OVER_PI);
        this.env.set("inve", Maja.ONE_OVER_E);
        this.env.set("ln2", Maja.LN2);
        this.env.set("ln10", Maja.LN10);
        this.env.set("log2e", Maja.LOG2E);
        this.env.set("mills", Maja.MILLS_CONSTANT);
        this.env.set("golomb_dickman", Maja.GOLOMB_DICKMAN_CONSTANT);
        this.env.set("deci", Maja.DECI);
        this.env.set("centi", Maja.CENTI);
        this.env.set("milli", Maja.MILLI);
        this.env.set("micro", Maja.MICRO);
        this.env.set("nano", Maja.NANO);
        this.env.set("pico", Maja.PICO);
        this.env.set("femto", Maja.FEMTO);
        this.env.set("atto", Maja.ATTO);
        this.env.set("zepto", Maja.ZEPTO);
        this.env.set("yocto", Maja.YOCTO);
        this.env.set("deca", Maja.DECA);
        this.env.set("hecto", Maja.HECTO);
        this.env.set("kilo", Maja.KILO);
        this.env.set("mega", Maja.MEGA);
        this.env.set("giga", Maja.GIGA);
        this.env.set("tera", Maja.TERA);
        this.env.set("peta", Maja.PETA);
        this.env.set("exa", Maja.EXA);
        this.env.set("zetta", Maja.ZETTA);
        this.env.set("yotta", Maja.YOTTA);

        // Regular trigonometric functions.
        this.env.set("sin", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.sin(c);
                } else if (x instanceof Double d) {
                    return Maja.sin(d);
                } else if (x instanceof Long l) {
                    return Maja.sin(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::sin);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::sin);
                } else {
                    throw new RuntimeException("Invalid argument type for sin(x).");
                }
            }
        });

        this.env.set("cos", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.cos(c);
                } else if (x instanceof Double d) {
                    return Maja.cos(d);
                } else if (x instanceof Long l) {
                    return Maja.cos(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::cos);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::cos);
                } else {
                    throw new RuntimeException("Invalid argument type for cos(x).");
                }
            }
        });

        env.set("tan", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.tan(c);
                } else if (x instanceof Double d) {
                    return Maja.tan(d);
                } else if (x instanceof Long l) {
                    return Maja.tan(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::tan);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::tan);
                } else {
                    throw new RuntimeException("Invalid argument type for tan(x).");
                }
            }
        });

        env.set("cot", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.cot(c);
                } else if (x instanceof Double d) {
                    return Maja.cot(d);
                } else if (x instanceof Long l) {
                    return Maja.cot(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::cot);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::cot);
                } else {
                    throw new RuntimeException("Invalid argument type for cot(x).");
                }
            }
        });

        env.set("sec", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.sec(c);
                } else if (x instanceof Double d) {
                    return Maja.sec(d);
                } else if (x instanceof Long l) {
                    return Maja.sec(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::sec);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::sec);
                } else {
                    throw new RuntimeException("Invalid argument type for sec(x).");
                }
            }
        });

        env.set("csc", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.csc(c);
                } else if (x instanceof Double d) {
                    return Maja.csc(d);
                } else if (x instanceof Long l) {
                    return Maja.csc(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::csc);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::csc);
                } else {
                    throw new RuntimeException("Invalid argument type for csc(x).");
                }
            }
        });

        // Inverse ("cyclometric") functions.
        this.env.set("asin", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.asin(c);
                } else if (x instanceof Double d) {
                    return Maja.asin(d);
                } else if (x instanceof Long l) {
                    return Maja.asin(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::asin);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::asin);
                } else {
                    throw new RuntimeException("Invalid argument type for asin(x).");
                }
            }
        });

        this.env.set("acos", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.acos(c);
                } else if (x instanceof Double d) {
                    return Maja.acos(d);
                } else if (x instanceof Long l) {
                    return Maja.acos(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::acos);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::acos);
                } else {
                    throw new RuntimeException("Invalid argument type for acos(x).");
                }
            }
        });

        env.set("atan", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.atan(c);
                } else if (x instanceof Double d) {
                    return Maja.atan(d);
                } else if (x instanceof Long l) {
                    return Maja.atan(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::atan);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::atan);
                } else {
                    throw new RuntimeException("Invalid argument type for atan(x).");
                }
            }
        });

        env.set("acot", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.acot(c);
                } else if (x instanceof Double d) {
                    return Maja.acot(d);
                } else if (x instanceof Long l) {
                    return Maja.acot(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::acot);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::acot);
                } else {
                    throw new RuntimeException("Invalid argument type for acot(x).");
                }
            }
        });

        env.set("asec", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.asec(c);
                } else if (x instanceof Double d) {
                    return Maja.asec(d);
                } else if (x instanceof Long l) {
                    return Maja.asec(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::asec);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::asec);
                } else {
                    throw new RuntimeException("Invalid argument type for asec(x).");
                }
            }
        });

        env.set("acsc", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.acsc(c);
                } else if (x instanceof Double d) {
                    return Maja.acsc(d);
                } else if (x instanceof Long l) {
                    return Maja.acsc(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::acsc);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::acsc);
                } else {
                    throw new RuntimeException("Invalid argument type for acsc(x).");
                }
            }
        });

        // Hyperbolic functions
        this.env.set("sinh", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.sinh(c);
                } else if (x instanceof Double d) {
                    return Maja.sinh(d);
                } else if (x instanceof Long l) {
                    return Maja.sinh(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::sinh);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::sinh);
                } else {
                    throw new RuntimeException("Invalid argument type for sinh(x).");
                }
            }
        });

        this.env.set("cosh", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.cosh(c);
                } else if (x instanceof Double d) {
                    return Maja.cosh(d);
                } else if (x instanceof Long l) {
                    return Maja.cosh(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::cosh);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::cosh);
                } else {
                    throw new RuntimeException("Invalid argument type for cosh(x).");
                }
            }
        });

        env.set("tanh", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.tanh(c);
                } else if (x instanceof Double d) {
                    return Maja.tanh(d);
                } else if (x instanceof Long l) {
                    return Maja.tanh(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::tanh);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::tanh);
                } else {
                    throw new RuntimeException("Invalid argument type for tanh(x).");
                }
            }
        });

        env.set("coth", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.coth(c);
                } else if (x instanceof Double d) {
                    return Maja.coth(d);
                } else if (x instanceof Long l) {
                    return Maja.coth(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::coth);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::coth);
                } else {
                    throw new RuntimeException("Invalid argument type for coth(x).");
                }
            }
        });

        env.set("sech", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.sech(c);
                } else if (x instanceof Double d) {
                    return Maja.sech(d);
                } else if (x instanceof Long l) {
                    return Maja.sech(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::sech);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::sech);
                } else {
                    throw new RuntimeException("Invalid argument type for sech(x).");
                }
            }
        });

        env.set("csch", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.csch(c);
                } else if (x instanceof Double d) {
                    return Maja.csch(d);
                } else if (x instanceof Long l) {
                    return Maja.csch(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::csch);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::csch);
                } else {
                    throw new RuntimeException("Invalid argument type for csch(x).");
                }
            }
        });

        // Inverse hyperbolic functions
        this.env.set("asinh", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.asinh(c);
                } else if (x instanceof Double d) {
                    return Maja.asinh(d);
                } else if (x instanceof Long l) {
                    return Maja.asinh(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::asinh);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::asinh);
                } else {
                    throw new RuntimeException("Invalid argument type for asinh(x).");
                }
            }
        });

        this.env.set("acosh", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.acosh(c);
                } else if (x instanceof Double d) {
                    return Maja.acosh(d);
                } else if (x instanceof Long l) {
                    return Maja.acosh(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::acosh);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::acosh);
                } else {
                    throw new RuntimeException("Invalid argument type for acosh(x).");
                }
            }
        });

        env.set("atanh", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.atanh(c);
                } else if (x instanceof Double d) {
                    return Maja.atanh(d);
                } else if (x instanceof Long l) {
                    return Maja.atanh(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::atanh);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::atanh);
                } else {
                    throw new RuntimeException("Invalid argument type for atanh(x).");
                }
            }
        });

        env.set("acoth", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.acoth(c);
                } else if (x instanceof Double d) {
                    return Maja.acoth(d);
                } else if (x instanceof Long l) {
                    return Maja.acoth(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::acoth);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::acoth);
                } else {
                    throw new RuntimeException("Invalid argument type for acoth(x).");
                }
            }
        });

        env.set("asech", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.asech(c);
                } else if (x instanceof Double d) {
                    return Maja.asech(d);
                } else if (x instanceof Long l) {
                    return Maja.asech(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::asech);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::asech);
                } else {
                    throw new RuntimeException("Invalid argument type for asech(x).");
                }
            }
        });

        env.set("acsch", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.acsch(c);
                } else if (x instanceof Double d) {
                    return Maja.acsch(d);
                } else if (x instanceof Long l) {
                    return Maja.acsch(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::acsch);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::acsch);
                } else {
                    throw new RuntimeException("Invalid argument type for acsch(x).");
                }
            }
        });

        // Min, max, abs, signum
        env.set("min", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x", "y");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x"), y = getEnv().get("y");
                if (x instanceof Long a && y instanceof Long b) {
                    return Math.min(a, b);
                } else if (x instanceof Double a && y instanceof Double b) {
                    return Math.min(a, b);
                } else if (x instanceof Double a && y instanceof Long b) {
                    return Math.min(a, b);
                } else if (x instanceof Long a && y instanceof Double b) {
                    return Math.min(a, b);
                } else if (x instanceof DoubleMatrix a && y instanceof DoubleMatrix b) {
                    return a.zipWith(b, Math::min);
                } else if (x instanceof DoubleMatrix a && y instanceof Long b) {
                    return a.map(e -> Math.min(e, b));
                } else if (x instanceof Long a && y instanceof DoubleMatrix b) {
                    return b.map(e -> Math.min(a, e));
                } else if (x instanceof Double a && y instanceof DoubleMatrix b) {
                    return b.map(e -> Math.min(a, e));
                } else if (x instanceof DoubleMatrix a && y instanceof Double b) {
                    return a.map(e -> Math.min(e, b));
                } else {
                    throw new RuntimeException("Invalid argument types for min(x, y).");
                }
            }
        });

        env.set("max", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x", "y");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x"), y = getEnv().get("y");
                if (x instanceof Long a && y instanceof Long b) {
                    return Math.max(a, b);
                } else if (x instanceof Double a && y instanceof Double b) {
                    return Math.max(a, b);
                } else if (x instanceof Double a && y instanceof Long b) {
                    return Math.max(a, b);
                } else if (x instanceof Long a && y instanceof Double b) {
                    return Math.max(a, b);
                } else if (x instanceof DoubleMatrix a && y instanceof DoubleMatrix b) {
                    return a.zipWith(b, Math::max);
                } else if (x instanceof DoubleMatrix a && y instanceof Long b) {
                    return a.map(e -> Math.max(e, b));
                } else if (x instanceof Long a && y instanceof DoubleMatrix b) {
                    return b.map(e -> Math.max(a, e));
                } else if (x instanceof Double a && y instanceof DoubleMatrix b) {
                    return b.map(e -> Math.max(a, e));
                } else if (x instanceof DoubleMatrix a && y instanceof Double b) {
                    return a.map(e -> Math.max(e, b));
                } else {
                    throw new RuntimeException("Invalid argument types for max(x, y).");
                }
            }
        });

        env.set("abs", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.abs(c);
                } else if (x instanceof Double d) {
                    return Maja.abs(d);
                } else if (x instanceof Long l) {
                    return Maja.abs(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::abs);
                } else if (x instanceof ComplexMatrix cm) {
                    return DoubleMatrix.into(cm.retype(Maja::abs));
                } else {
                    throw new RuntimeException("Invalid argument type for abs(x).");
                }
            }
        });

        env.set("signum", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Double d) {
                    return Maja.signum(d);
                } else if (x instanceof Long l) {
                    return Maja.signum(l);
                } else if (x instanceof Complex c) {
                    return Maja.signum(c);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::signum);
                } else if (x instanceof ComplexMatrix dm) {
                    return dm.map(Maja::signum);
                } else {
                    throw new RuntimeException("Invalid argument type for signum(x).");
                }
            }
        });

        // Exp, log
        env.set("exp", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.exp(c);
                } else if (x instanceof Double d) {
                    return Maja.exp(d);
                } else if (x instanceof Long l) {
                    return Maja.exp(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::exp);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::exp);
                } else {
                    throw new RuntimeException("Invalid argument type for exp(x).");
                }
            }
        });

        this.env.set("ln", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            private Object transform(Object x) {
                if(anyComplex(x))
                    return Maja.log((Complex) x);
                else if(allDouble(x)) {
                    try {
                        double r = Maja.log(coerceDouble(x));
                        if (isPathologic(r))
                            return Maja.log(new Complex(coerceDouble(x)));
                        else
                            return r;
                    } catch (ArithmeticException e) {
                        return Maja.log(new Complex(coerceDouble(x)));
                    }
                } else {
                    throw new RuntimeException("Invalid argument type for ln(x).");
                }
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if(x instanceof ComplexMatrix cm) {
                    return cm.map(z -> forceComplex(transform(z)));
                } else if(x instanceof DoubleMatrix dm) {
                    // Note: Will be transformed into an uniform matrix upon simplification.
                    return dm.retype(this::transform);
                } else {
                    return transform(x);
                }
            }
        });

        // TODO: arbitrary base logarithm goes here.

        env.set("log10", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Double d) {
                    return Maja.log10(d);
                } else if (x instanceof Long l) {
                    return Maja.log10(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::log10);
                } else {
                    throw new RuntimeException("Invalid argument type for log10(x).");
                }
            }
        });

        env.set("log2", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Double d) {
                    return Maja.log2(d);
                } else if (x instanceof Long l) {
                    return Maja.log2(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::log2);
                } else {
                    throw new RuntimeException("Invalid argument type for log2(x).");
                }
            }
        });

        env.set("log1p", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Double d) {
                    return Maja.log1p(d);
                } else if (x instanceof Long l) {
                    return Maja.log1p(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::log1p);
                } else {
                    throw new RuntimeException("Invalid argument type for log1p(x).");
                }
            }
        });

        this.env.set("sqrt", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            private Object transform(Object x) {
                if(anyComplex(x))
                    return Maja.sqrt((Complex) x);
                else if(allDouble(x)) {
                    try {
                        double r = Maja.sqrt(coerceDouble(x));
                        if (isPathologic(r))
                            return Maja.sqrt(new Complex(coerceDouble(x)));
                        else
                            return r;
                    } catch (ArithmeticException e) {
                        return Maja.sqrt(new Complex(coerceDouble(x)));
                    }
                } else {
                    throw new RuntimeException("Invalid argument type for sqrt(x).");
                }
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if(x instanceof ComplexMatrix cm) {
                    return cm.map(z -> forceComplex(transform(z)));
                } else if(x instanceof DoubleMatrix dm) {
                    // Note: Will be transformed into an uniform matrix upon simplification.
                    return dm.retype(this::transform);
                } else {
                    return transform(x);
                }
            }
        });

        this.env.set("cbrt", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            private Object transform(Object x) {
                if(anyComplex(x))
                    return Maja.cbrt((Complex) x);
                else if(allDouble(x)) {
                    try {
                        double r = Maja.cbrt(coerceDouble(x));
                        if (isPathologic(r))
                            return Maja.cbrt(new Complex(coerceDouble(x)));
                        else
                            return r;
                    } catch (ArithmeticException e) {
                        return Maja.cbrt(new Complex(coerceDouble(x)));
                    }
                } else {
                    throw new RuntimeException("Invalid argument type for cbrt(x).");
                }
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if(x instanceof ComplexMatrix cm) {
                    return cm.map(z -> forceComplex(transform(z)));
                } else if(x instanceof DoubleMatrix dm) {
                    // Note: Will be transformed into an uniform matrix upon simplification.
                    return dm.retype(this::transform);
                } else {
                    return transform(x);
                }
            }
        });

        env.set("hypot", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x", "y");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x"), y = getEnv().get("y");
                if (x instanceof Long lx && y instanceof Long ly) {
                    return Maja.hypot(lx, ly);
                } else if (x instanceof Double dx && y instanceof Double dy) {
                    return Maja.hypot(dx, dy);
                } else if (x instanceof DoubleMatrix dxm && y instanceof DoubleMatrix dym) {
                    return dxm.zipWith(dym, Maja::hypot);
                } else if (x instanceof Long a && y instanceof Double b) {
                    return Maja.hypot(a, b);
                } else if (x instanceof Double a && y instanceof Long b) {
                    return Maja.hypot(a, b);
                } else if (x instanceof DoubleMatrix dxm && y instanceof Double d) {
                    return dxm.map(a -> Maja.hypot(a, d));
                } else if (x instanceof Double a && y instanceof DoubleMatrix dym) {
                    return dym.map(b -> Maja.hypot(a, b));
                } else if (x instanceof DoubleMatrix dxm && y instanceof Long d) {
                    return dxm.map(a -> Maja.hypot(a, d));
                } else if (x instanceof Long a && y instanceof DoubleMatrix dym) {
                    return dym.map(b -> Maja.hypot(a, b));
                } else {
                    throw new RuntimeException("Invalid argument type for hypot(x, y).");
                }
            }
        });

        // Rounding.
        env.set("ceil", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.ceil(c);
                } else if (x instanceof Double d) {
                    return Maja.ceil(d);
                } else if (x instanceof Long l) {
                    return Maja.ceil(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::ceil);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::ceil);
                } else {
                    throw new RuntimeException("Invalid argument type for ceil(x).");
                }
            }
        });

        env.set("floor", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.floor(c);
                } else if (x instanceof Double d) {
                    return Maja.floor(d);
                } else if (x instanceof Long l) {
                    return Maja.floor(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::floor);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::floor);
                } else {
                    throw new RuntimeException("Invalid argument type for floor(x).");
                }
            }
        });

        env.set("round", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.round(c);
                } else if (x instanceof Double d) {
                    return Maja.round(d);
                } else if (x instanceof Long l) {
                    return Maja.round(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(z -> (double) Maja.round(z));
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::round);
                } else {
                    throw new RuntimeException("Invalid argument type for round(x).");
                }
            }
        });

        // Angle functions.
        env.set("atan2", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x", "y");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x"), y = getEnv().get("y");
                if (x instanceof Double d && y instanceof Double e) {
                    return Math.atan2(d, e);
                } else if (x instanceof Long c && y instanceof Long d) {
                    return Maja.atan2(c, d);
                } else if (x instanceof Long a && y instanceof Double b) {
                    return Maja.atan2(a, b);
                } else if (x instanceof Double a && y instanceof Long b) {
                    return Maja.atan2(a, b);
                } else if (x instanceof DoubleMatrix dm && y instanceof Double b) {
                    return dm.map(z -> Maja.atan2(z, b));
                } else if (x instanceof Double a && y instanceof DoubleMatrix dm) {
                    return dm.map(z -> Maja.atan2(a, z));
                } else if (x instanceof DoubleMatrix dm && y instanceof DoubleMatrix em) {
                    return dm.zipWith(em, Maja::atan2);
                } else {
                    throw new RuntimeException("Invalid argument type for atan2(x, y).");
                }
            }
        });

        env.set("sinc", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Complex c) {
                    return Maja.sinc(c);
                } else if (x instanceof Double d) {
                    return Maja.sinc(d);
                } else if (x instanceof Long l) {
                    return Maja.sinc(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::sinc);
                } else if (x instanceof ComplexMatrix cm) {
                    return cm.map(Maja::sinc);
                } else {
                    throw new RuntimeException("Invalid argument type for sinc(x).");
                }
            }
        });

        env.set("rad", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Double d) {
                    return Maja.toRadians(d);
                } else if (x instanceof Long l) {
                    return Maja.toRadians(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::toRadians);
                } else {
                    throw new RuntimeException("Invalid argument type for rad(x).");
                }
            }
        });

        env.set("deg", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Double d) {
                    return Maja.toDegrees(d);
                } else if (x instanceof Long l) {
                    return Maja.toDegrees(l);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Maja::toDegrees);
                } else {
                    throw new RuntimeException("Invalid argument type for deg(x).");
                }
            }
        });

        env.set("ulp", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Double d) {
                    return Math.ulp(d);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(Math::ulp);
                } else {
                    throw new RuntimeException("Invalid argument type for deg(x).");
                }
            }
        });

        env.set("scalb", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x", "y");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x"), y = getEnv().get("y");
                if (x instanceof Double a && y instanceof Long b) {
                    return Maja.scalb(a, b.intValue());
                } else if (x instanceof DoubleMatrix dm && y instanceof Long b) {
                    return dm.map(z -> Math.scalb(z, b.intValue()));
                } else {
                    throw new RuntimeException("Invalid argument type for scalb(x, y).");
                }
            }
        });

        env.set("copysign", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x", "y");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x"), y = getEnv().get("y");
                if (x instanceof Double a && y instanceof Long b) {
                    return Maja.copySign(a, b);
                } else if (x instanceof Long a && y instanceof Double b) {
                    return Maja.copySign(a, b);
                } else if (x instanceof Long a && y instanceof Long b) {
                    return Maja.copySign(a, b);
                } else if (x instanceof Double a && y instanceof Double b) {
                    return Maja.copySign(a, b);
                } else if (x instanceof Complex a && y instanceof Complex b) {
                    return Maja.copySign(a, b);
                } else if (x instanceof ComplexMatrix a && y instanceof ComplexMatrix b) {
                    return a.zipWith(b, Maja::copySign);
                } else if (x instanceof ComplexMatrix a && y instanceof Complex b) {
                    return a.map(z -> Maja.copySign(z, b));
                } else if (x instanceof Complex a && y instanceof ComplexMatrix b) {
                    return b.map(z -> Maja.copySign(a, z));
                } else if (x instanceof DoubleMatrix a && y instanceof DoubleMatrix b) {
                    return a.zipWith(b, Maja::copySign);
                } else if (x instanceof DoubleMatrix a && y instanceof Double b) {
                    return a.map(z -> Maja.copySign(z, b));
                } else if (x instanceof Double a && y instanceof DoubleMatrix b) {
                    return b.map(z -> Maja.copySign(a, z));
                } else if (x instanceof DoubleMatrix a && y instanceof Long b) {
                    return a.map(z -> Maja.copySign(z, b));
                } else if (x instanceof Long a && y instanceof DoubleMatrix b) {
                    return b.map(z -> Maja.copySign(a, z));
                } else {
                    throw new RuntimeException("Invalid argument type for copysign(x, y).");
                }
            }
        });

        env.set("getexp", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Double d) {
                    return Maja.getExponent(d);
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(z -> (double) Maja.getExponent(z));
                } else {
                    throw new RuntimeException("Invalid argument type for getexp(x).");
                }
            }
        });

        env.set("rand", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("begin", "end");
            }

            @Override
            public Object eval() {
                Object begin = getEnv().get("begin"), end = getEnv().get("end");
                if (begin instanceof Double a && end instanceof Double b) {
                    return Maja.random(a, b);
                } else if (begin instanceof DoubleMatrix a && end instanceof DoubleMatrix b) {
                    return a.zipWith(b, Maja::random);
                } else if (begin instanceof DoubleMatrix a && end instanceof Double b) {
                    return a.map(z -> Maja.random(z, b));
                } else if (begin instanceof Double a && end instanceof DoubleMatrix b) {
                    return b.map(z -> Maja.random(a, z));
                } else if (begin instanceof DoubleMatrix a && end instanceof Long b) {
                    return a.map(z -> Maja.random(z, b));
                } else if (begin instanceof Long a && end instanceof DoubleMatrix b) {
                    return b.map(z -> Maja.random(a, z));
                } else if (begin instanceof Long a && end instanceof Long b) {
                    return Maja.random(a, b);
                } else {
                    throw new RuntimeException("Invalid argument type for rand(begin, end).");
                }
            }
        });

        env.set("compare", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x", "y");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x"), y = getEnv().get("y");
                if (x instanceof Double a && y instanceof Double b) {
                    return Maja.compare(a, b);
                } else if (x instanceof DoubleMatrix a && y instanceof DoubleMatrix b) {
                    return a.zipWith(b, (i, j) -> (double) Maja.compare(i, j));
                } else if (x instanceof DoubleMatrix a && y instanceof Double b) {
                    return a.map(z -> (double) Maja.compare(z, b));
                } else if (x instanceof Double a && y instanceof DoubleMatrix b) {
                    return b.map(z -> (double) Maja.compare(a, z));
                } else if (x instanceof DoubleMatrix a && y instanceof Long b) {
                    return a.map(z -> (double) Maja.compare(z, b));
                } else if (x instanceof Long a && y instanceof DoubleMatrix b) {
                    return b.map(z -> (double) Maja.compare(a, z));
                } else if (x instanceof Long a && y instanceof Long b) {
                    return Maja.compare(a, b);
                } else {
                    throw new RuntimeException("Invalid argument type for compare(x, y).");
                }
            }
        });

        env.set("approx_eq", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x", "y", "tol");
            }

            @Override
            public Object eval() {
                Object a = getEnv().get("x"), b = getEnv().get("y"), tol = getEnv().get("tol");
                double dtol = tol instanceof Double d ? d : (tol instanceof Long l ? l : Double.NaN);
                if (Double.isNaN(dtol))
                    throw new RuntimeException("Invalid argument type for approx_eq(x, y, tol). Tolerance is NaN.");

                if (a instanceof Long l1 && b instanceof Long l2) {
                    return Maja.eq(l1, l2, dtol) ? 1 : 0;
                } else if (a instanceof Double d1 && b instanceof Double d2) {
                    return Maja.eq(d1, d2, dtol) ? 1 : 0;
                } else if (a instanceof Complex c1 && b instanceof Complex c2) {
                    return Maja.eq(c1, c2, dtol) ? 1 : 0;
                } else if (a instanceof Long l1 && b instanceof Double d2) {
                    return Maja.eq(l1, d2, dtol) ? 1 : 0;
                } else if (a instanceof Double d1 && b instanceof Long l2) {
                    return Maja.eq(d1, l2, dtol) ? 1 : 0;
                } else if (a instanceof Long l1 && b instanceof Complex c2) {
                    return Maja.eq(new Complex(l1), c2, dtol) ? 1 : 0;
                } else if (a instanceof Complex c1 && b instanceof Long l2) {
                    return Maja.eq(c1, new Complex(l2), dtol) ? 1 : 0;
                } else if (a instanceof Double d1 && b instanceof Complex c2) {
                    return Maja.eq(new Complex(d1), c2, dtol) ? 1 : 0;
                } else if (a instanceof Complex c1 && b instanceof Double d2) {
                    return Maja.eq(c1, new Complex(d2), dtol) ? 1 : 0;
                } else if (a instanceof DoubleMatrix dm1 && b instanceof Long l2) {
                    return dm1.map(d -> Maja.eq(d, l2, dtol) ? 1.0 : 0.0);
                } else if (a instanceof Long l1 && b instanceof DoubleMatrix dm2) {
                    return dm2.map(d -> Maja.eq(l1, d, dtol) ? 1.0 : 0.0);
                } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
                    return dm1.zipWith(dm2, (x, y) -> Maja.eq(x, y, dtol) ? 1.0 : 0.0);
                } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
                    return dm1.map(d -> Maja.eq(d, d2, dtol) ? 1.0 : 0.0);
                } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
                    return dm2.map(d -> Maja.eq(d1, d, dtol) ? 1.0 : 0.0);
                } else if (a instanceof DoubleMatrix dm1 && b instanceof Complex c2) {
                    if(c2.im() != 0) {
                        DoubleMatrix result = new DoubleMatrix(dm1.height(), dm1.width());
                        for(int i = 0; i < dm1.height(); i++) {
                            for(int j = 0; j < dm1.width(); j++) {
                                result.set(i, j, 0.0);
                            }
                        }
                        return result;
                    }
                    double c2re = c2.re();
                    return dm1.map(d -> Maja.eq(c2re, d, dtol) ? 1.0 : 0.0);
                } else if (a instanceof Complex c1 && b instanceof DoubleMatrix dm2) {
                    if(c1.im() != 0) {
                        DoubleMatrix result = new DoubleMatrix(dm2.height(), dm2.width());
                        for(int i = 0; i < dm2.height(); i++) {
                            for(int j = 0; j < dm2.width(); j++) {
                                result.set(i, j, 0.0);
                            }
                        }
                        return result;
                    }
                    double c1re = c1.re();
                    return dm2.map(d -> Maja.eq(c1re, d, dtol) ? 1.0 : 0.0);
                } else if (a instanceof ComplexMatrix dm1 && b instanceof Long l2) {
                    Complex cl2 = new Complex(l2);
                    return dm1.map(d -> Maja.eq(d, cl2, dtol) ? Complex.ONE : Complex.ZERO);
                } else if (a instanceof Long l1 && b instanceof ComplexMatrix dm2) {
                    Complex cl1 = new Complex(l1);
                    return dm2.map(d -> Maja.eq(cl1, d, dtol) ? Complex.ONE : Complex.ZERO);
                } else if (a instanceof ComplexMatrix dm1 && b instanceof ComplexMatrix dm2) {
                    return dm1.zipWith(dm2, (x, y) -> Maja.eq(x, y, dtol) ? Complex.ONE : Complex.ZERO);
                } else if (a instanceof ComplexMatrix dm1 && b instanceof Double d2) {
                    Complex cd2 = new Complex(d2);
                    return dm1.map(d -> Maja.eq(d, cd2, dtol) ? Complex.ONE : Complex.ZERO);
                } else if (a instanceof Double d1 && b instanceof ComplexMatrix dm2) {
                    Complex cd1 = new Complex(d1);
                    return dm2.map(d -> Maja.eq(cd1, d, dtol) ? Complex.ONE : Complex.ZERO);
                } else if (a instanceof ComplexMatrix dm1 && b instanceof Complex c2) {
                    return dm1.map(d -> Maja.eq(d, c2, dtol) ? Complex.ONE : Complex.ZERO);
                } else if (a instanceof Complex c1 && b instanceof ComplexMatrix dm2) {
                    return dm2.map(d -> Maja.eq(c1, d, dtol) ? Complex.ONE : Complex.ZERO);
                } else {
                    throw new RuntimeException("Invalid type for approx_eq(a, b, tol: double).");
                }
            }
        });

        env.set("is_perfect_square", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("n");
            }

            @Override
            public Object eval() {
                Object n = getEnv().get("n");
                if (n instanceof Long l) {
                    return Maja.isPerfectSquare(l) ? 1 : 0;
                } else {
                    throw new RuntimeException("Invalid type for is_perfect_square(n).");
                }
            }
        });

        env.set("linear_map", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("a_begin", "a_end", "b_begin", "b_end", "value");
            }

            @Override
            public Object eval() {
                double a_begin = coerceDouble(getEnv().get("a_begin"));
                double a_end = coerceDouble(getEnv().get("a_end"));
                double b_begin = coerceDouble(getEnv().get("b_begin"));
                double b_end = coerceDouble(getEnv().get("b_end"));
                double value = coerceDouble(getEnv().get("value"));
                return Maja.linearMap(a_begin, a_end, b_begin, b_end, value);
            }
        });

        env.set("linear_norm", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("begin", "end", "value");
            }

            @Override
            public Object eval() {
                double begin = coerceDouble(getEnv().get("begin"));
                double end = coerceDouble(getEnv().get("end"));
                double value = coerceDouble(getEnv().get("value"));
                return Maja.linearNorm(begin, end, value);
            }
        });

        env.set("linear_interpolate", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("begin", "end", "value");
            }

            @Override
            public Object eval() {
                double begin = coerceDouble(getEnv().get("begin"));
                double end = coerceDouble(getEnv().get("end"));
                double value = coerceDouble(getEnv().get("value"));
                return Maja.linearInterpolate(begin, end, value);
            }
        });

        env.set("clamp", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("begin", "end", "value");
            }

            @Override
            public Object eval() {
                Object begin = getEnv().get("begin");
                Object end = getEnv().get("end");
                Object value = getEnv().get("value");
                if(begin instanceof Long a && end instanceof Long b && value instanceof Long c) {
                    return Maja.clamp(a, b, c);
                } else {
                    return Maja.clamp(coerceDouble(begin), coerceDouble(end), coerceDouble(value));
                }
            }
        });

        env.set("is_power_of_2", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if(x instanceof Long l) {
                    return Maja.isPowerOfTwo(l);
                } else {
                    throw new RuntimeException("Invalid type for is_power_of_2(x).");
                }
            }
        });

        env.set("next_power_of_2", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if(x instanceof Long l) {
                    return Maja.isPowerOfTwo(l);
                } else {
                    throw new RuntimeException("Invalid type for next_power_of_2(x).");
                }
            }
        });

        this.env.set("fast_sin", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Double d) {
                    return (double) Maja.fastSin((float) d.doubleValue());
                } else if (x instanceof Long l) {
                    return (double) Maja.fastSin((float) l.doubleValue());
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(z -> (double) Maja.fastSin((float) z.doubleValue()));
                } else {
                    throw new RuntimeException("Invalid argument type for fast_sin(x).");
                }
            }
        });

        this.env.set("fast_cos", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if (x instanceof Double d) {
                    return (double) Maja.fastCos((float) d.doubleValue());
                } else if (x instanceof Long l) {
                    return (double) Maja.fastCos((float) l.doubleValue());
                } else if (x instanceof DoubleMatrix dm) {
                    return dm.map(z -> (double) Maja.fastCos((float) z.doubleValue()));
                } else {
                    throw new RuntimeException("Invalid argument type for fast_cos(x).");
                }
            }
        });

        // Integer functions.
        this.env.set("isqrt", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if(x instanceof Long l) {
                    return (long) Maja.isqrt(l);
                } else if(x instanceof Double d) {
                    return (long) Maja.isqrt(d.longValue());
                } else if(x instanceof DoubleMatrix dm) {
                    return dm.map(z -> (double) Maja.isqrt(z.longValue()));
                } else {
                    throw new RuntimeException("Invalid argument type for isqrt(x).");
                }
            }
        });

        this.env.set("ilog10", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if(x instanceof Long l) {
                    return (long) Maja.ilog10(Math.toIntExact(l));
                } else if(x instanceof Double d) {
                    return (long) Maja.ilog10(d.intValue());
                } else if(x instanceof DoubleMatrix dm) {
                    return dm.map(z -> (double) Maja.ilog10(z.intValue()));
                } else {
                    throw new RuntimeException("Invalid argument type for ilog10(x).");
                }
            }
        });
        this.env.set("iexp", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x", "y");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x"), y = getEnv().get("y");
                long ix, iy;
                if (x instanceof Long l) ix = l; else if (x instanceof Double d) ix = d.longValue(); else throw new RuntimeException("Invalid argument type for iexp(x, y).");
                if (y instanceof Long l) iy = l; else if (y instanceof Double d) iy = d.longValue(); else throw new RuntimeException("Invalid argument type for iexp(x, y).");
                return Maja.ipow(ix, iy);
            }
        });
        this.env.set("icbrt", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if(x instanceof Long l) {
                    return Maja.icbrt(l);
                } else if(x instanceof Double d) {
                    return Maja.icbrt(d.longValue());
                } else if(x instanceof DoubleMatrix dm) {
                    return dm.map(z -> (double) Maja.icbrt(z.longValue()));
                } else {
                    throw new RuntimeException("Invalid argument type for icbrt(x).");
                }
            }
        });

        this.env.set("airy_ai", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            private Object transform(Object x) {
                if(anyComplex(x))
                    return Maja.airyAi((Complex) x);
                else if(allDouble(x)) {
                    try {
                        double r = Maja.airyAi(coerceDouble(x));
                        if (isPathologic(r))
                            return Maja.airyAi(new Complex(coerceDouble(x)));
                        else
                            return r;
                    } catch (ArithmeticException e) {
                        return Maja.airyAi(new Complex(coerceDouble(x)));
                    }
                } else {
                    throw new RuntimeException("Invalid argument type for airy_ai(x).");
                }
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if(x instanceof ComplexMatrix cm) {
                    return cm.map(z -> forceComplex(transform(z)));
                } else if(x instanceof DoubleMatrix dm) {
                    // Note: Will be transformed into an uniform matrix upon simplification.
                    return dm.retype(this::transform);
                } else {
                    return transform(x);
                }
            }
        });

        this.env.set("airy_bi", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            private Object transform(Object x) {
                if(anyComplex(x))
                    return Maja.airyBi((Complex) x);
                else if(allDouble(x)) {
                    try {
                        double r = Maja.airyBi(coerceDouble(x));
                        if (isPathologic(r))
                            return Maja.airyBi(new Complex(coerceDouble(x)));
                        else
                            return r;
                    } catch (ArithmeticException e) {
                        return Maja.airyBi(new Complex(coerceDouble(x)));
                    }
                } else {
                    throw new RuntimeException("Invalid argument type for airy_bi(x).");
                }
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if(x instanceof ComplexMatrix cm) {
                    return cm.map(z -> forceComplex(transform(z)));
                } else if(x instanceof DoubleMatrix dm) {
                    // Note: Will be transformed into an uniform matrix upon simplification.
                    return dm.retype(this::transform);
                } else {
                    return transform(x);
                }
            }
        });

        this.env.set("airy_aip", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            private Object transform(Object x) {
                if(anyComplex(x))
                    return Maja.airyAip((Complex) x);
                else if(allDouble(x)) {
                    try {
                        double r = Maja.airyAip(coerceDouble(x));
                        if (isPathologic(r))
                            return Maja.airyAip(new Complex(coerceDouble(x)));
                        else
                            return r;
                    } catch (ArithmeticException e) {
                        return Maja.airyAip(new Complex(coerceDouble(x)));
                    }
                } else {
                    throw new RuntimeException("Invalid argument type for airy_aip(x).");
                }
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if(x instanceof ComplexMatrix cm) {
                    return cm.map(z -> forceComplex(transform(z)));
                } else if(x instanceof DoubleMatrix dm) {
                    // Note: Will be transformed into an uniform matrix upon simplification.
                    return dm.retype(this::transform);
                } else {
                    return transform(x);
                }
            }
        });

        this.env.set("airy_bip", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            private Object transform(Object x) {
                if(anyComplex(x))
                    return Maja.airyBip((Complex) x);
                else if(allDouble(x)) {
                    try {
                        double r = Maja.airyBip(coerceDouble(x));
                        if (isPathologic(r))
                            return Maja.airyBip(new Complex(coerceDouble(x)));
                        else
                            return r;
                    } catch (ArithmeticException e) {
                        return Maja.airyBip(new Complex(coerceDouble(x)));
                    }
                } else {
                    throw new RuntimeException("Invalid argument type for airy_bip(x).");
                }
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if(x instanceof ComplexMatrix cm) {
                    return cm.map(z -> forceComplex(transform(z)));
                } else if(x instanceof DoubleMatrix dm) {
                    // Note: Will be transformed into an uniform matrix upon simplification.
                    return dm.retype(this::transform);
                } else {
                    return transform(x);
                }
            }
        });

        this.env.set("gamma", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            private Object transform(Object x) {
                if(anyComplex(x))
                    return Maja.gamma((Complex) x);
                else if(allDouble(x)) {
                    try {
                        double r = Maja.gamma(coerceDouble(x));
                        if (isPathologic(r))
                            return Maja.gamma(new Complex(coerceDouble(x)));
                        else
                            return r;
                    } catch (ArithmeticException e) {
                        return Maja.gamma(new Complex(coerceDouble(x)));
                    }
                } else {
                    throw new RuntimeException("Invalid argument type for gamma(x).");
                }
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if(x instanceof ComplexMatrix cm) {
                    return cm.map(z -> forceComplex(transform(z)));
                } else if(x instanceof DoubleMatrix dm) {
                    // Note: Will be transformed into an uniform matrix upon simplification.
                    return dm.retype(this::transform);
                } else {
                    return transform(x);
                }
            }
        });

        this.env.set("loggamma", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            private Object transform(Object x) {
                if(anyComplex(x))
                    return Maja.loggamma((Complex) x);
                else if(allDouble(x)) {
                    try {
                        double r = Maja.loggamma(coerceDouble(x));
                        if (isPathologic(r))
                            return Maja.loggamma(new Complex(coerceDouble(x)));
                        else
                            return r;
                    } catch (ArithmeticException e) {
                        return Maja.loggamma(new Complex(coerceDouble(x)));
                    }
                } else {
                    throw new RuntimeException("Invalid argument type for loggamma(x).");
                }
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if(x instanceof ComplexMatrix cm) {
                    return cm.map(z -> forceComplex(transform(z)));
                } else if(x instanceof DoubleMatrix dm) {
                    // Note: Will be transformed into an uniform matrix upon simplification.
                    return dm.retype(this::transform);
                } else {
                    return transform(x);
                }
            }
        });

        this.env.set("digamma", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            private Object transform(Object x) {
                if(anyComplex(x))
                    return Maja.digamma((Complex) x);
                else if(allDouble(x)) {
                    try {
                        double r = Maja.digamma(coerceDouble(x));
                        if (isPathologic(r))
                            return Maja.digamma(new Complex(coerceDouble(x)));
                        else
                            return r;
                    } catch (ArithmeticException e) {
                        return Maja.digamma(new Complex(coerceDouble(x)));
                    }
                } else {
                    throw new RuntimeException("Invalid argument type for digamma(x).");
                }
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if(x instanceof ComplexMatrix cm) {
                    return cm.map(z -> forceComplex(transform(z)));
                } else if(x instanceof DoubleMatrix dm) {
                    // Note: Will be transformed into an uniform matrix upon simplification.
                    return dm.retype(this::transform);
                } else {
                    return transform(x);
                }
            }
        });

        this.env.set("trigamma", new ExpressionFunction() {
            @Override
            public List<String> params() {
                return List.of("x");
            }

            private Object transform(Object x) {
                if(anyComplex(x))
                    return Maja.trigamma((Complex) x);
                else if(allDouble(x)) {
                    try {
                        double r = Maja.trigamma(coerceDouble(x));
                        if (isPathologic(r))
                            return Maja.trigamma(new Complex(coerceDouble(x)));
                        else
                            return r;
                    } catch (ArithmeticException e) {
                        return Maja.trigamma(new Complex(coerceDouble(x)));
                    }
                } else {
                    throw new RuntimeException("Invalid argument type for trigamma(x).");
                }
            }

            @Override
            public Object eval() {
                Object x = getEnv().get("x");
                if(x instanceof ComplexMatrix cm) {
                    return cm.map(z -> forceComplex(transform(z)));
                } else if(x instanceof DoubleMatrix dm) {
                    // Note: Will be transformed into an uniform matrix upon simplification.
                    return dm.retype(this::transform);
                } else {
                    return transform(x);
                }
            }
        });
    }

    private static Complex forceComplex(Object d) {
        if (d instanceof Complex c) {
            return c;
        } else if (d instanceof Double) {
            return new Complex((double) d);
        } else if (d instanceof Long) {
            return new Complex((long) d);
        } else {
            throw new RuntimeException("Could not coerce the value to a complex number.");
        }
    }

    private static boolean isPathologic(double d) {
        return Double.isNaN(d) || Double.isInfinite(d);
    }

    private static boolean anyComplex(Object... objs) {
        for (Object obj : objs) {
            if (obj instanceof Complex) {
                return true;
            }
        }
        return false;
    }

    private static boolean allDouble(Object... objs) {
        for (Object obj : objs) {
            if (!(obj instanceof Double) && !(obj instanceof Long)) {
                return false;
            }
        }
        return true;
    }

    private static double coerceDouble(Object obj) {
        if (obj instanceof Double d) {
            return d;
        } else if (obj instanceof Long l) {
            return l;
        } else {
            return Double.NaN;
        }
    }

    @Override
    public List<Object> visitMain(ExpressionParser.MainContext ctx) {
        return ctx.toplevel().stream().map(this::visit).toList();
    }

    @Override
    public Object visitToplevel(ExpressionParser.ToplevelContext ctx) {
        return visit(ctx.getChild(0));
    }

    @Override
    public Object visitBlock(ExpressionParser.BlockContext ctx) {
        Environment old = env; env = env.createChild();
        Object last = null;
        for (int i = 0; i < ctx.toplevel().size(); i++) {
            if (ctx.toplevel(i).expression() != null) {
                return visit(ctx.toplevel(i).expression());
            } else {
                last = visit(ctx.toplevel(i));
            }
        }
        env = old;
        return last;
    }

    @Override
    public Object visitSimpleAssignment(ExpressionParser.SimpleAssignmentContext ctx) {
        Object value = visit(ctx.expression());
        env.set(ctx.ID().getText(), value);
        return value;
    }

    @Override
    public Object visitSimpleFunctionDeclaration(ExpressionParser.SimpleFunctionDeclarationContext ctx) {
        List<String> args = ctx.ID().stream().skip(1).map(ParseTree::getText).toList();
        var f = new ExpressionFunction() {
            @Override
            public List<String> params() {
                return args;
            }

            @Override
            public Object eval() {
                return visit(ctx.expression());
            }
        };
        env.set(ctx.ID(0).getText(), f);
        return f;
    }

    @Override
    public Object visitFunctionDeclaration(ExpressionParser.FunctionDeclarationContext ctx) {
        List<String> args = ctx.ID().stream().skip(1).map(ParseTree::getText).toList();
        var f = new ExpressionFunction() {
            @Override
            public List<String> params() {
                return args;
            }

            @Override
            public Object eval() {
                try {
                    Object res = visit(ctx.block());
                    if (res == null)
                        throw new RuntimeException("Function " + ctx.ID(0).getText() + " did not return a value.");
                    return res;
                } catch (ReturnError err) {
                    return err.value;
                }
            }
        };
        env.set(ctx.ID(0).getText(), f);
        return f;
    }

    private boolean truthy(Object cond) {
        boolean b;
        if (cond instanceof Long l) {
            b = l != 0;
        } else if (cond instanceof Double d) {
            b = d != 0;
        } else if (cond instanceof Complex c) {
            b = Maja.ne(c, Complex.ZERO);
        } else {
            throw new RuntimeException("Invalid condition type.");
        }
        return b;
    }

    @Override
    public Object visitIf(ExpressionParser.IfContext ctx) {
        Object cond = visit(ctx.expression());
        if (truthy(cond)) {
            return visit(ctx.block(0));
        } else if (ctx.block().size() > 1) {
            return visit(ctx.block(1));
        }
        return null;
    }

    @Override
    public Object visitSimpleIf(ExpressionParser.SimpleIfContext ctx) {
        Object cond = visit(ctx.expression(0));
        if (truthy(cond)) {
            return visit(ctx.expression(1));
        } else if (ctx.expression().size() > 2) {
            return visit(ctx.expression(2));
        }
        return null;
    }

    @Override
    public Object visitWhile(ExpressionParser.WhileContext ctx) {
        Object r = null;
        while(truthy(visit(ctx.expression())))
            r = visit(ctx.block());
        return r;
    }

    private static Object determineStep(Object begin, Object end) {
        if (begin instanceof Long l1 && end instanceof Long l2) {
            return l1 < l2 ? 1L : -1L;
        } else if (begin instanceof Double d1 && end instanceof Double d2) {
            return d1 < d2 ? 1.0 : -1.0;
        } else if (begin instanceof Long l1 && end instanceof Double d2) {
            return l1 < d2 ? 1.0 : -1.0;
        } else if (begin instanceof Double d1 && end instanceof Long l2) {
            return d1 < l2 ? 1.0 : -1.0;
        } else {
            throw new RuntimeException("Invalid type for for loop.");
        }
    }

    private static boolean mayFor(Object cur, Object end) {
        if (cur instanceof Long l1 && end instanceof Long l2) {
            return l1 < l2;
        } else if (cur instanceof Double d1 && end instanceof Double d2) {
            return d1 < d2;
        } else if (cur instanceof Long l1 && end instanceof Double d2) {
            return l1 < d2;
        } else if (cur instanceof Double d1 && end instanceof Long l2) {
            return d1 < l2;
        } else {
            throw new RuntimeException("Invalid type for for loop.");
        }
    }

    private static Object addFor(Object cur, Object step) {
        if (cur instanceof Long l1 && step instanceof Long l2) {
            return l1 + l2;
        } else if (cur instanceof Double d1 && step instanceof Double d2) {
            return d1 + d2;
        } else if (cur instanceof Long l1 && step instanceof Double d2) {
            return l1 + d2;
        } else if (cur instanceof Double d1 && step instanceof Long l2) {
            return d1 + l2;
        } else {
            throw new RuntimeException("Invalid type for for loop.");
        }
    }

    @Override
    public Object visitFor(ExpressionParser.ForContext ctx) {
        String objname = ctx.ID().getText();
        Object initial = visit(ctx.expression(0));
        Object end = visit(ctx.expression(1));
        Object step = ctx.expression().size() > 2 ? visit(ctx.expression(2)) : determineStep(initial, end);
        Environment old = env; env = env.createChild();
        env.set(objname, initial);
        Object r = null;
        while(mayFor(env.get(objname), end)) {
            r = visit(ctx.block());
            env.set(objname, addFor(env.get(objname), step));
        }
        env = old;
        return r;
    }

    @Override
    public Object visitReturn(ExpressionParser.ReturnContext ctx) {
        throw new ReturnError(visit(ctx.expression()));
    }

    @Override
    public Object visitExprGcd(ExpressionParser.ExprGcdContext ctx) {
        Object a = visit(ctx.expression(0)), b = visit(ctx.expression(1));
        if (a instanceof Long l1 && b instanceof Long l2) {
            return Maja.gcd(l1, l2);
        } else if (a instanceof Double d1 && b instanceof Double d2) {
            return Maja.gcd(d1, d2);
        } else if (a instanceof Complex c1 && b instanceof Complex c2) {
            return Maja.gcd(c1, c2);
        } else if (a instanceof Long l1 && b instanceof Double d2) {
            return Maja.gcd(l1, d2);
        } else if (a instanceof Double d1 && b instanceof Long l2) {
            return Maja.gcd(d1, l2);
        } else if (a instanceof Long l1 && b instanceof Complex c2) {
            return Maja.gcd(new Complex(l1), c2);
        } else if (a instanceof Complex c1 && b instanceof Long l2) {
            return Maja.gcd(c1, new Complex(l2));
        } else if (a instanceof Double d1 && b instanceof Complex c2) {
            return Maja.gcd(new Complex(d1), c2);
        } else if (a instanceof Complex c1 && b instanceof Double d2) {
            return Maja.gcd(c1, new Complex(d2));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Long l2) {
            return dm1.map(d -> Maja.gcd(d, l2));
        } else if (a instanceof Long l1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.gcd(l1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
            return dm1.zipWith(dm2, Maja::gcd);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> Maja.gcd(d, d2));
        } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.gcd(d1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Complex c2) {
            return ComplexMatrix.into(dm1.retype(Complex::new).map(c -> Maja.gcd(c, c2)));
        } else if (a instanceof Complex c1 && b instanceof DoubleMatrix dm2) {
            return ComplexMatrix.into(dm2.retype(Complex::new).map(c -> Maja.gcd(c1, c)));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Long l2) {
            Complex cl2 = new Complex(l2);
            return dm1.map(d -> Maja.gcd(d, cl2));
        } else if (a instanceof Long l1 && b instanceof ComplexMatrix dm2) {
            Complex cl1 = new Complex(l1);
            return dm2.map(d -> Maja.gcd(cl1, d));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof ComplexMatrix dm2) {
            return dm1.zipWith(dm2, Maja::gcd);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Double d2) {
            Complex cd2 = new Complex(d2);
            return dm1.map(d -> Maja.gcd(d, cd2));
        } else if (a instanceof Double d1 && b instanceof ComplexMatrix dm2) {
            Complex cd1 = new Complex(d1);
            return dm2.map(d -> Maja.gcd(cd1, d));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Complex c2) {
            return dm1.map(d -> Maja.gcd(d, c2));
        } else if (a instanceof Complex c1 && b instanceof ComplexMatrix dm2) {
            return dm2.map(d -> Maja.gcd(c1, d));
        } else {
            throw new RuntimeException("Invalid type for gcd.");
        }
    }

    @Override
    public Object visitExprIndex(ExpressionParser.ExprIndexContext ctx) {
        Object a = visit(ctx.expression(0));
        List b = ctx.expression().stream().skip(1).map(this::visit).collect(Collectors.toList());
        boolean allLong = b.stream().allMatch(o -> o instanceof Long);
        if (a instanceof DoubleMatrix dm && b.size() > 2 && allLong) {
            return dm.get((int) b.get(0), (int) b.get(1));
        } else if (a instanceof ComplexMatrix cm && b.size() > 2 && allLong) {
            return cm.get((int) b.get(0), (int) b.get(1));
        } else {
            throw new RuntimeException("Invalid type for index.");
        }
    }

    @Override
    public Object visitExprNeg(ExpressionParser.ExprNegContext ctx) {
        Object a = visit(ctx.expression());
        if (a instanceof Long l) {
            return -l;
        } else if (a instanceof Double d) {
            return -d;
        } else if (a instanceof Complex c) {
            return Maja.negate(c);
        } else if (a instanceof DoubleMatrix dm) {
            return dm.map(d -> -d);
        } else if (a instanceof ComplexMatrix cm) {
            return cm.map(Maja::negate);
        } else {
            throw new RuntimeException("Invalid type for -.");
        }
    }

    @Override
    public Object visitExprNot(ExpressionParser.ExprNotContext ctx) {
        Object a = visit(ctx.expression());
        if (a instanceof Long l) {
            return (l == 0) ? 1 : 0;
        } else if (a instanceof Double d) {
            return (d == 0) ? 1 : 0;
        } else if (a instanceof Complex c) {
            return Maja.eq(c, Complex.ZERO) ? 1 : 0;
        } else if (a instanceof DoubleMatrix dm) {
            return dm.map(d -> (d == 0) ? 1.0 : 0.0);
        } else if (a instanceof ComplexMatrix cm) {
            return cm.map(c -> Maja.eq(c, Complex.ZERO) ? Complex.ONE : Complex.ZERO);
        } else {
            throw new RuntimeException("Invalid type for ~.");
        }
    }

    @Override
    public Object visitExprPos(ExpressionParser.ExprPosContext ctx) {
        return visit(ctx.expression());
    }

    @Override
    public Object visitExprDiv(ExpressionParser.ExprDivContext ctx) {
        Object a = visit(ctx.expression(0)), b = visit(ctx.expression(1));
        if (a instanceof Long l1 && b instanceof Long l2) {
            return Maja.div(l1, l2);
        } else if (a instanceof Double d1 && b instanceof Double d2) {
            return Maja.div(d1, d2);
        } else if (a instanceof Complex c1 && b instanceof Complex c2) {
            return Maja.div(c1, c2);
        } else if (a instanceof Long l1 && b instanceof Double d2) {
            return Maja.div(l1, d2);
        } else if (a instanceof Double d1 && b instanceof Long l2) {
            return Maja.div(d1, l2);
        } else if (a instanceof Long l1 && b instanceof Complex c2) {
            return Maja.div(new Complex(l1), c2);
        } else if (a instanceof Complex c1 && b instanceof Long l2) {
            return Maja.div(c1, new Complex(l2));
        } else if (a instanceof Double d1 && b instanceof Complex c2) {
            return Maja.div(new Complex(d1), c2);
        } else if (a instanceof Complex c1 && b instanceof Double d2) {
            return Maja.div(c1, new Complex(d2));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Long l2) {
            return dm1.map(d -> Maja.div(d, l2));
        } else if (a instanceof Long l1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.div(l1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
            return dm1.zipWith(dm2, Maja::div);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> Maja.div(d, d2));
        } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.div(d1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Complex c2) {
            return ComplexMatrix.into(dm1.retype(Complex::new).map(c -> Maja.div(c, c2)));
        } else if (a instanceof Complex c1 && b instanceof DoubleMatrix dm2) {
            return ComplexMatrix.into(dm2.retype(Complex::new).map(c -> Maja.div(c1, c)));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Long l2) {
            Complex cl2 = new Complex(l2);
            return dm1.map(d -> Maja.div(d, cl2));
        } else if (a instanceof Long l1 && b instanceof ComplexMatrix dm2) {
            Complex cl1 = new Complex(l1);
            return dm2.map(d -> Maja.div(cl1, d));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof ComplexMatrix dm2) {
            return dm1.zipWith(dm2, Maja::div);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Double d2) {
            Complex cd2 = new Complex(d2);
            return dm1.map(d -> Maja.div(d, cd2));
        } else if (a instanceof Double d1 && b instanceof ComplexMatrix dm2) {
            Complex cd1 = new Complex(d1);
            return dm2.map(d -> Maja.div(cd1, d));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Complex c2) {
            return dm1.map(d -> Maja.div(d, c2));
        } else if (a instanceof Complex c1 && b instanceof ComplexMatrix dm2) {
            return dm2.map(d -> Maja.div(c1, d));
        } else {
            throw new RuntimeException("Invalid type for /.");
        }
    }

    @Override
    public Object visitExprOr(ExpressionParser.ExprOrContext ctx) {
        Object a = visit(ctx.expression(0)), b = visit(ctx.expression(1));
        if (a instanceof Long l1 && b instanceof Long l2) {
            return (l1 != 0 || l2 != 0) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Double d2) {
            return (d1 != 0 || d2 != 0) ? 1 : 0;
        } else if (a instanceof Complex c1 && b instanceof Complex c2) {
            return (Maja.ne(c1, Complex.ZERO) || Maja.ne(c2, Complex.ZERO)) ? 1 : 0;
        } else if (a instanceof Long l1 && b instanceof Double d2) {
            return (l1 != 0 || d2 != 0) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Long l2) {
            return (d1 != 0 || l2 != 0) ? 1 : 0;
        } else if (a instanceof Long l1 && b instanceof Complex c2) {
            return (Maja.ne(new Complex(l1), Complex.ZERO) || Maja.ne(c2, Complex.ZERO)) ? 1 : 0;
        } else if (a instanceof Complex c1 && b instanceof Long l2) {
            return (Maja.ne(c1, Complex.ZERO) || Maja.ne(new Complex(l2), Complex.ZERO)) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Complex c2) {
            return (Maja.ne(new Complex(d1), Complex.ZERO) || Maja.ne(c2, Complex.ZERO)) ? 1 : 0;
        } else if (a instanceof Complex c1 && b instanceof Double d2) {
            return (Maja.ne(c1, Complex.ZERO) || Maja.ne(new Complex(d2), Complex.ZERO)) ? 1 : 0;
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Long l2) {
            return dm1.map(d -> (d != 0 || l2 != 0) ? 1.0 : 0.0);
        } else if (a instanceof Long l1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> (l1 != 0 || d != 0) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
            return dm1.zipWith(dm2, (x, y) -> (x != 0 || y != 0) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> (d != 0 || d2 != 0) ? 1.0 : 0.0);
        } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> (d1 != 0 || d != 0) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Complex c2) {
            return dm1.map(d -> (Maja.ne(c2, Complex.ZERO) || d != 0) ? 1.0 : 0.0);
        } else if (a instanceof Complex c1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> (Maja.ne(c1, Complex.ZERO) || d != 0) ? 1.0 : 0.0);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Long l2) {
            return dm1.map(d -> (Maja.ne(d, Complex.ZERO) || l2 != 0) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof Long l1 && b instanceof ComplexMatrix dm2) {
            return dm2.map(d -> (Maja.ne(d, Complex.ZERO) || l1 != 0) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof ComplexMatrix dm2) {
            return dm1.zipWith(dm2, (x, y) -> (Maja.ne(x, Complex.ZERO) || Maja.ne(y, Complex.ZERO)) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> (Maja.ne(d, Complex.ZERO) || d2 != 0) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof Double d1 && b instanceof ComplexMatrix dm2) {
            return dm2.map(d -> (Maja.ne(d, Complex.ZERO) || d1 != 0) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Complex c2) {
            return dm1.map(d -> (Maja.ne(d, Complex.ZERO) || Maja.ne(c2, Complex.ZERO)) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof Complex c1 && b instanceof ComplexMatrix dm2) {
            return dm2.map(d -> (Maja.ne(d, Complex.ZERO) || Maja.ne(c1, Complex.ZERO)) ? Complex.ONE : Complex.ZERO);
        } else {
            throw new RuntimeException("Invalid type for ||.");
        }
    }

    @Override
    public Object visitExprSub(ExpressionParser.ExprSubContext ctx) {
        Object a = visit(ctx.expression(0)), b = visit(ctx.expression(1));
        if (a instanceof Long l1 && b instanceof Long l2) {
            return Maja.sub(l1, l2);
        } else if (a instanceof Double d1 && b instanceof Double d2) {
            return Maja.sub(d1, d2);
        } else if (a instanceof Complex c1 && b instanceof Complex c2) {
            return Maja.sub(c1, c2);
        } else if (a instanceof Long l1 && b instanceof Double d2) {
            return Maja.sub(l1, d2);
        } else if (a instanceof Double d1 && b instanceof Long l2) {
            return Maja.sub(d1, l2);
        } else if (a instanceof Long l1 && b instanceof Complex c2) {
            return Maja.sub(new Complex(l1), c2);
        } else if (a instanceof Complex c1 && b instanceof Long l2) {
            return Maja.sub(c1, new Complex(l2));
        } else if (a instanceof Double d1 && b instanceof Complex c2) {
            return Maja.sub(new Complex(d1), c2);
        } else if (a instanceof Complex c1 && b instanceof Double d2) {
            return Maja.sub(c1, new Complex(d2));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Long l2) {
            return dm1.map(d -> Maja.sub(d, l2));
        } else if (a instanceof Long l1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.sub(l1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
            return dm1.zipWith(dm2, Maja::sub);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> Maja.sub(d, d2));
        } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.sub(d1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Complex c2) {
            return ComplexMatrix.into(dm1.retype(Complex::new).map(c -> Maja.sub(c, c2)));
        } else if (a instanceof Complex c1 && b instanceof DoubleMatrix dm2) {
            return ComplexMatrix.into(dm2.retype(Complex::new).map(c -> Maja.sub(c1, c)));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Long l2) {
            Complex cl2 = new Complex(l2);
            return dm1.map(d -> Maja.sub(d, cl2));
        } else if (a instanceof Long l1 && b instanceof ComplexMatrix dm2) {
            Complex cl1 = new Complex(l1);
            return dm2.map(d -> Maja.sub(cl1, d));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof ComplexMatrix dm2) {
            return dm1.zipWith(dm2, Maja::sub);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Double d2) {
            Complex cd2 = new Complex(d2);
            return dm1.map(d -> Maja.sub(d, cd2));
        } else if (a instanceof Double d1 && b instanceof ComplexMatrix dm2) {
            Complex cd1 = new Complex(d1);
            return dm2.map(d -> Maja.sub(cd1, d));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Complex c2) {
            return dm1.map(d -> Maja.sub(d, c2));
        } else if (a instanceof Complex c1 && b instanceof ComplexMatrix dm2) {
            return dm2.map(d -> Maja.sub(c1, d));
        } else {
            throw new RuntimeException("Invalid type for -.");
        }
    }

    @Override
    public Object visitExprMul(ExpressionParser.ExprMulContext ctx) {
        Object a = visit(ctx.expression(0)), b = visit(ctx.expression(1));
        if (a instanceof Long l1 && b instanceof Long l2) {
            return Maja.mul(l1, l2);
        } else if (a instanceof Double d1 && b instanceof Double d2) {
            return Maja.mul(d1, d2);
        } else if (a instanceof Complex c1 && b instanceof Complex c2) {
            return Maja.mul(c1, c2);
        } else if (a instanceof Long l1 && b instanceof Double d2) {
            return Maja.mul(l1, d2);
        } else if (a instanceof Double d1 && b instanceof Long l2) {
            return Maja.mul(d1, l2);
        } else if (a instanceof Long l1 && b instanceof Complex c2) {
            return Maja.mul(new Complex(l1), c2);
        } else if (a instanceof Complex c1 && b instanceof Long l2) {
            return Maja.mul(c1, new Complex(l2));
        } else if (a instanceof Double d1 && b instanceof Complex c2) {
            return Maja.mul(new Complex(d1), c2);
        } else if (a instanceof Complex c1 && b instanceof Double d2) {
            return Maja.mul(c1, new Complex(d2));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Long l2) {
            return dm1.map(d -> Maja.mul(d, l2));
        } else if (a instanceof Long l1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.mul(l1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
            return dm1.zipWith(dm2, Maja::mul);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> Maja.mul(d, d2));
        } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.mul(d1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Complex c2) {
            return ComplexMatrix.into(dm1.retype(Complex::new).map(c -> Maja.mul(c, c2)));
        } else if (a instanceof Complex c1 && b instanceof DoubleMatrix dm2) {
            return ComplexMatrix.into(dm2.retype(Complex::new).map(c -> Maja.mul(c1, c)));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Long l2) {
            Complex cl2 = new Complex(l2);
            return dm1.map(d -> Maja.mul(d, cl2));
        } else if (a instanceof Long l1 && b instanceof ComplexMatrix dm2) {
            Complex cl1 = new Complex(l1);
            return dm2.map(d -> Maja.mul(cl1, d));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof ComplexMatrix dm2) {
            return dm1.zipWith(dm2, Maja::mul);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Double d2) {
            Complex cd2 = new Complex(d2);
            return dm1.map(d -> Maja.mul(d, cd2));
        } else if (a instanceof Double d1 && b instanceof ComplexMatrix dm2) {
            Complex cd1 = new Complex(d1);
            return dm2.map(d -> Maja.mul(cd1, d));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Complex c2) {
            return dm1.map(d -> Maja.mul(d, c2));
        } else if (a instanceof Complex c1 && b instanceof ComplexMatrix dm2) {
            return dm2.map(d -> Maja.mul(c1, d));
        } else {
            throw new RuntimeException("Invalid type for *.");
        }
    }

    @Override
    public Object visitExprGe(ExpressionParser.ExprGeContext ctx) {
        Object a = visit(ctx.expression(0)), b = visit(ctx.expression(1));
        if (a instanceof Long l1 && b instanceof Long l2) {
            return Maja.ge(l1, l2) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Double d2) {
            return Maja.ge(d1, d2) ? 1 : 0;
        } else if (a instanceof Long l1 && b instanceof Double d2) {
            return Maja.ge(l1, d2) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Long l2) {
            return Maja.ge(d1, l2) ? 1 : 0;
        } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
            return dm1.zipWith(dm2, (d1, d2) -> Maja.ge(d1, d2) ? 1.0 : 0.0);
        } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.ge(d1, d) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> Maja.ge(d, d2) ? 1.0 : 0.0);
        } else if (a instanceof Long d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.ge(d1, d) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Long d2) {
            return dm1.map(d -> Maja.ge(d, d2) ? 1.0 : 0.0);
        } else {
            throw new RuntimeException("Invalid type for >=.");
        }
    }

    @Override
    public Object visitExprMod(ExpressionParser.ExprModContext ctx) {
        Object a = visit(ctx.expression(0)), b = visit(ctx.expression(1));
        if (a instanceof Long l1 && b instanceof Long l2) {
            return Maja.mod(l1, l2);
        } else if (a instanceof Double d1 && b instanceof Double d2) {
            return Maja.mod(d1, d2);
        } else if (a instanceof Long l1 && b instanceof Double d2) {
            return Maja.mod(l1, d2);
        } else if (a instanceof Double d1 && b instanceof Long l2) {
            return Maja.mod(d1, l2);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
            return dm1.zipWith(dm2, Maja::mod);
        } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.mod(d1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> Maja.mod(d, d2));
        } else if (a instanceof Long d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.mod(d1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Long d2) {
            return dm1.map(d -> Maja.mod(d, d2));
        } else {
            throw new RuntimeException("Invalid type for mod.");
        }
    }

    @Override
    public Object visitExprLe(ExpressionParser.ExprLeContext ctx) {
        Object a = visit(ctx.expression(0)), b = visit(ctx.expression(1));
        if (a instanceof Long l1 && b instanceof Long l2) {
            return Maja.le(l1, l2) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Double d2) {
            return Maja.le(d1, d2) ? 1 : 0;
        } else if (a instanceof Long l1 && b instanceof Double d2) {
            return Maja.le(l1, d2) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Long l2) {
            return Maja.le(d1, l2) ? 1 : 0;
        } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
            return dm1.zipWith(dm2, (d1, d2) -> Maja.le(d1, d2) ? 1.0 : 0.0);
        } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.le(d1, d) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> Maja.le(d, d2) ? 1.0 : 0.0);
        } else if (a instanceof Long d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.le(d1, d) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Long d2) {
            return dm1.map(d -> Maja.le(d, d2) ? 1.0 : 0.0);
        } else {
            throw new RuntimeException("Invalid type for <=.");
        }
    }

    @Override
    public Object visitExprParen(ExpressionParser.ExprParenContext ctx) {
        return visit(ctx.expression());
    }

    @Override
    public Object visitExprInt(ExpressionParser.ExprIntContext ctx) {
        return Long.parseLong(ctx.INT().getText());
    }

    @Override
    public Object visitMatrix(ExpressionParser.MatrixContext ctx) {
        return ctx.expression().stream().map(this::visit).collect(Collectors.toList());
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    @Override
    public Object visitExprMatrix(ExpressionParser.ExprMatrixContext ctx) {
        List constructor = ctx.matrix().stream().map(this::visit).map(List.class::cast).collect(Collectors.toList());
        if (constructor.stream().allMatch(l -> ((List) l).stream().allMatch(o -> o instanceof Double || o instanceof Long))) {
            return DoubleMatrix.into(new Matrix(constructor).retype(x -> x instanceof Long ? ((Long) x).doubleValue() : x));
        } else if (constructor.stream().allMatch(l -> ((List) l).stream().allMatch(o -> o instanceof Complex))) {
            return new ComplexMatrix((List<List<Complex>>) constructor);
        } else {
            return new Matrix(constructor);
        }
    }

    @Override
    public Object visitExprGt(ExpressionParser.ExprGtContext ctx) {
        Object a = visit(ctx.expression(0)), b = visit(ctx.expression(1));
        if (a instanceof Long l1 && b instanceof Long l2) {
            return Maja.gt(l1, l2) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Double d2) {
            return Maja.gt(d1, d2) ? 1 : 0;
        } else if (a instanceof Long l1 && b instanceof Double d2) {
            return Maja.gt(l1, d2) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Long l2) {
            return Maja.gt(d1, l2) ? 1 : 0;
        } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
            return dm1.zipWith(dm2, (d1, d2) -> Maja.gt(d1, d2) ? 1.0 : 0.0);
        } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.gt(d1, d) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> Maja.gt(d, d2) ? 1.0 : 0.0);
        } else if (a instanceof Long d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.gt(d1, d) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Long d2) {
            return dm1.map(d -> Maja.gt(d, d2) ? 1.0 : 0.0);
        } else {
            throw new RuntimeException("Invalid type for >.");
        }
    }

    @Override
    public Object visitExprEq(ExpressionParser.ExprEqContext ctx) {
        Object a = visit(ctx.expression(0)), b = visit(ctx.expression(1));
        if (a instanceof Long l1 && b instanceof Long l2) {
            return Maja.eq(l1, l2) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Double d2) {
            return Maja.eq(d1, d2) ? 1 : 0;
        } else if (a instanceof Complex c1 && b instanceof Complex c2) {
            return Maja.eq(c1, c2) ? 1 : 0;
        } else if (a instanceof Long l1 && b instanceof Double d2) {
            return Maja.eq(l1, d2) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Long l2) {
            return Maja.eq(d1, l2) ? 1 : 0;
        } else if (a instanceof Long l1 && b instanceof Complex c2) {
            return Maja.eq(new Complex(l1), c2) ? 1 : 0;
        } else if (a instanceof Complex c1 && b instanceof Long l2) {
            return Maja.eq(c1, new Complex(l2)) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Complex c2) {
            return Maja.eq(new Complex(d1), c2) ? 1 : 0;
        } else if (a instanceof Complex c1 && b instanceof Double d2) {
            return Maja.eq(c1, new Complex(d2)) ? 1 : 0;
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Long l2) {
            return dm1.map(d -> Maja.eq(d, l2) ? 1.0 : 0.0);
        } else if (a instanceof Long l1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.eq(l1, d) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
            return dm1.zipWith(dm2, (x, y) -> Maja.eq(x, y) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> Maja.eq(d, d2) ? 1.0 : 0.0);
        } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.eq(d1, d) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Complex c2) {
            return dm1.map(d -> Maja.eq(c2, d) ? 1.0 : 0.0);
        } else if (a instanceof Complex c1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.eq(c1, d) ? 1.0 : 0.0);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Long l2) {
            Complex cl2 = new Complex(l2);
            return dm1.map(d -> Maja.eq(d, cl2) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof Long l1 && b instanceof ComplexMatrix dm2) {
            Complex cl1 = new Complex(l1);
            return dm2.map(d -> Maja.eq(cl1, d) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof ComplexMatrix dm2) {
            return dm1.zipWith(dm2, (x, y) -> Maja.eq(x, y) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Double d2) {
            Complex cd2 = new Complex(d2);
            return dm1.map(d -> Maja.eq(d, cd2) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof Double d1 && b instanceof ComplexMatrix dm2) {
            Complex cd1 = new Complex(d1);
            return dm2.map(d -> Maja.eq(cd1, d) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Complex c2) {
            return dm1.map(d -> Maja.eq(d, c2) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof Complex c1 && b instanceof ComplexMatrix dm2) {
            return dm2.map(d -> Maja.eq(c1, d) ? Complex.ONE : Complex.ZERO);
        } else {
            throw new RuntimeException("Invalid type for ==.");
        }
    }

    @Override
    public Object visitExprAnd(ExpressionParser.ExprAndContext ctx) {
        Object a = visit(ctx.expression(0)), b = visit(ctx.expression(1));
        if (a instanceof Long l1 && b instanceof Long l2) {
            return (l1 != 0 && l2 != 0) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Double d2) {
            return (d1 != 0 && d2 != 0) ? 1 : 0;
        } else if (a instanceof Complex c1 && b instanceof Complex c2) {
            return (Maja.ne(c1, Complex.ZERO) && Maja.ne(c2, Complex.ZERO)) ? 1 : 0;
        } else if (a instanceof Long l1 && b instanceof Double d2) {
            return (l1 != 0 && d2 != 0) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Long l2) {
            return (d1 != 0 && l2 != 0) ? 1 : 0;
        } else if (a instanceof Long l1 && b instanceof Complex c2) {
            return (Maja.ne(new Complex(l1), Complex.ZERO) && Maja.ne(c2, Complex.ZERO)) ? 1 : 0;
        } else if (a instanceof Complex c1 && b instanceof Long l2) {
            return (Maja.ne(c1, Complex.ZERO) && Maja.ne(new Complex(l2), Complex.ZERO)) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Complex c2) {
            return (Maja.ne(new Complex(d1), Complex.ZERO) && Maja.ne(c2, Complex.ZERO)) ? 1 : 0;
        } else if (a instanceof Complex c1 && b instanceof Double d2) {
            return (Maja.ne(c1, Complex.ZERO) && Maja.ne(new Complex(d2), Complex.ZERO)) ? 1 : 0;
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Long l2) {
            return dm1.map(d -> (d != 0 && l2 != 0) ? 1.0 : 0.0);
        } else if (a instanceof Long l1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> (l1 != 0 && d != 0) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
            return dm1.zipWith(dm2, (x, y) -> (x != 0 && y != 0) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> (d != 0 && d2 != 0) ? 1.0 : 0.0);
        } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> (d1 != 0 && d != 0) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Complex c2) {
            return dm1.map(d -> (Maja.ne(c2, Complex.ZERO) && d != 0) ? 1.0 : 0.0);
        } else if (a instanceof Complex c1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> (Maja.ne(c1, Complex.ZERO) && d != 0) ? 1.0 : 0.0);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Long l2) {
            return dm1.map(d -> (Maja.ne(d, Complex.ZERO) && l2 != 0) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof Long l1 && b instanceof ComplexMatrix dm2) {
            return dm2.map(d -> (Maja.ne(d, Complex.ZERO) && l1 != 0) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof ComplexMatrix dm2) {
            return dm1.zipWith(dm2, (x, y) -> (Maja.ne(x, Complex.ZERO) && Maja.ne(y, Complex.ZERO)) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> (Maja.ne(d, Complex.ZERO) && d2 != 0) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof Double d1 && b instanceof ComplexMatrix dm2) {
            return dm2.map(d -> (Maja.ne(d, Complex.ZERO) && d1 != 0) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Complex c2) {
            return dm1.map(d -> (Maja.ne(d, Complex.ZERO) && Maja.ne(c2, Complex.ZERO)) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof Complex c1 && b instanceof ComplexMatrix dm2) {
            return dm2.map(d -> (Maja.ne(d, Complex.ZERO) && Maja.ne(c1, Complex.ZERO)) ? Complex.ONE : Complex.ZERO);
        } else {
            throw new RuntimeException("Invalid type for &&.");
        }
    }

    @Override
    public Object visitExprPow(ExpressionParser.ExprPowContext ctx) {
        Object a = visit(ctx.expression(0)), b = visit(ctx.expression(1));
        if (a instanceof Long l1 && b instanceof Long l2) {
            return Maja.pow(l1, l2);
        } else if (a instanceof Double d1 && b instanceof Double d2) {
            return Maja.pow(d1, d2);
        } else if (a instanceof Complex c1 && b instanceof Complex c2) {
            return Maja.pow(c1, c2);
        } else if (a instanceof Long l1 && b instanceof Double d2) {
            return Maja.pow(l1, d2);
        } else if (a instanceof Double d1 && b instanceof Long l2) {
            return Maja.pow(d1, l2);
        } else if (a instanceof Long l1 && b instanceof Complex c2) {
            return Maja.pow(new Complex(l1), c2);
        } else if (a instanceof Complex c1 && b instanceof Long l2) {
            return Maja.pow(c1, new Complex(l2));
        } else if (a instanceof Double d1 && b instanceof Complex c2) {
            return Maja.pow(new Complex(d1), c2);
        } else if (a instanceof Complex c1 && b instanceof Double d2) {
            return Maja.pow(c1, new Complex(d2));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Long l2) {
            return dm1.map(d -> Maja.pow(d, l2));
        } else if (a instanceof Long l1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.pow(l1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
            return dm1.zipWith(dm2, Maja::pow);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> Maja.pow(d, d2));
        } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.pow(d1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Complex c2) {
            return ComplexMatrix.into(dm1.retype(Complex::new).map(c -> Maja.pow(c, c2)));
        } else if (a instanceof Complex c1 && b instanceof DoubleMatrix dm2) {
            return ComplexMatrix.into(dm2.retype(Complex::new).map(c -> Maja.pow(c1, c)));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Long l2) {
            Complex cl2 = new Complex(l2);
            return dm1.map(d -> Maja.pow(d, cl2));
        } else if (a instanceof Long l1 && b instanceof ComplexMatrix dm2) {
            Complex cl1 = new Complex(l1);
            return dm2.map(d -> Maja.pow(cl1, d));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof ComplexMatrix dm2) {
            return dm1.zipWith(dm2, Maja::pow);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Double d2) {
            Complex cd2 = new Complex(d2);
            return dm1.map(d -> Maja.pow(d, cd2));
        } else if (a instanceof Double d1 && b instanceof ComplexMatrix dm2) {
            Complex cd1 = new Complex(d1);
            return dm2.map(d -> Maja.pow(cd1, d));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Complex c2) {
            return dm1.map(d -> Maja.pow(d, c2));
        } else if (a instanceof Complex c1 && b instanceof ComplexMatrix dm2) {
            return dm2.map(d -> Maja.pow(c1, d));
        } else {
            throw new RuntimeException("Invalid type for **.");
        }
    }

    @Override
    public Object visitExprFunctionCall(ExpressionParser.ExprFunctionCallContext ctx) {
        String function = ctx.ID().getText();
        Object callable = env.get(function);
        if (callable instanceof ExpressionFunction e) {
            if(e.params().size() != ctx.expression().size())
                throw new RuntimeException("Invalid number of arguments for function " + function + ".");
            Environment old = env;
            env = env.createChild();
            for(int i = 0; i < e.params().size(); i++)
                env.set(e.params().get(i), visit(ctx.expression(i)));
            Object result = e.eval();
            env = old;
            return result;
        } else {
            throw new RuntimeException("Invalid function call.");
        }
    }

    @Override
    public Object visitExprLcm(ExpressionParser.ExprLcmContext ctx) {
        Object a = visit(ctx.expression(0)), b = visit(ctx.expression(1));
        if (a instanceof Long l1 && b instanceof Long l2) {
            return Maja.lcm(l1, l2);
        } else if (a instanceof Double d1 && b instanceof Double d2) {
            return Maja.lcm(d1, d2);
        } else if (a instanceof Complex c1 && b instanceof Complex c2) {
            return Maja.lcm(c1, c2);
        } else if (a instanceof Long l1 && b instanceof Double d2) {
            return Maja.lcm(l1, d2);
        } else if (a instanceof Double d1 && b instanceof Long l2) {
            return Maja.lcm(d1, l2);
        } else if (a instanceof Long l1 && b instanceof Complex c2) {
            return Maja.lcm(new Complex(l1), c2);
        } else if (a instanceof Complex c1 && b instanceof Long l2) {
            return Maja.lcm(c1, new Complex(l2));
        } else if (a instanceof Double d1 && b instanceof Complex c2) {
            return Maja.lcm(new Complex(d1), c2);
        } else if (a instanceof Complex c1 && b instanceof Double d2) {
            return Maja.lcm(c1, new Complex(d2));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Long l2) {
            return dm1.map(d -> Maja.lcm(d, l2));
        } else if (a instanceof Long l1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.lcm(l1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
            return dm1.zipWith(dm2, Maja::lcm);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> Maja.lcm(d, d2));
        } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.lcm(d1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Complex c2) {
            return ComplexMatrix.into(dm1.retype(Complex::new).map(c -> Maja.lcm(c, c2)));
        } else if (a instanceof Complex c1 && b instanceof DoubleMatrix dm2) {
            return ComplexMatrix.into(dm2.retype(Complex::new).map(c -> Maja.lcm(c1, c)));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Long l2) {
            Complex cl2 = new Complex(l2);
            return dm1.map(d -> Maja.lcm(d, cl2));
        } else if (a instanceof Long l1 && b instanceof ComplexMatrix dm2) {
            Complex cl1 = new Complex(l1);
            return dm2.map(d -> Maja.lcm(cl1, d));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof ComplexMatrix dm2) {
            return dm1.zipWith(dm2, Maja::lcm);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Double d2) {
            Complex cd2 = new Complex(d2);
            return dm1.map(d -> Maja.lcm(d, cd2));
        } else if (a instanceof Double d1 && b instanceof ComplexMatrix dm2) {
            Complex cd1 = new Complex(d1);
            return dm2.map(d -> Maja.lcm(cd1, d));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Complex c2) {
            return dm1.map(d -> Maja.lcm(d, c2));
        } else if (a instanceof Complex c1 && b instanceof ComplexMatrix dm2) {
            return dm2.map(d -> Maja.lcm(c1, d));
        } else {
            throw new RuntimeException("Invalid type for lcm.");
        }
    }

    @Override
    public Object visitExprLt(ExpressionParser.ExprLtContext ctx) {
        Object a = visit(ctx.expression(0)), b = visit(ctx.expression(1));
        if (a instanceof Long l1 && b instanceof Long l2) {
            return Maja.lt(l1, l2) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Double d2) {
            return Maja.lt(d1, d2) ? 1 : 0;
        } else if (a instanceof Long l1 && b instanceof Double d2) {
            return Maja.lt(l1, d2) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Long l2) {
            return Maja.lt(d1, l2) ? 1 : 0;
        } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
            return dm1.zipWith(dm2, (d1, d2) -> Maja.lt(d1, d2) ? 1.0 : 0.0);
        } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.lt(d1, d) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> Maja.lt(d, d2) ? 1.0 : 0.0);
        } else if (a instanceof Long d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.lt(d1, d) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Long d2) {
            return dm1.map(d -> Maja.lt(d, d2) ? 1.0 : 0.0);
        } else {
            throw new RuntimeException("Invalid type for <.");
        }
    }

    @Override
    public Object visitExprRem(ExpressionParser.ExprRemContext ctx) {
        Object a = visit(ctx.expression(0)), b = visit(ctx.expression(1));
        if (a instanceof Long l1 && b instanceof Long l2) {
            return Maja.rem(l1, l2);
        } else if (a instanceof Double d1 && b instanceof Double d2) {
            return Maja.rem(d1, d2);
        } else if (a instanceof Complex c1 && b instanceof Complex c2) {
            return Maja.rem(c1, c2);
        } else if (a instanceof Long l1 && b instanceof Double d2) {
            return Maja.rem(l1, d2);
        } else if (a instanceof Double d1 && b instanceof Long l2) {
            return Maja.rem(d1, l2);
        } else if (a instanceof Long l1 && b instanceof Complex c2) {
            return Maja.rem(new Complex(l1), c2);
        } else if (a instanceof Complex c1 && b instanceof Long l2) {
            return Maja.rem(c1, new Complex(l2));
        } else if (a instanceof Double d1 && b instanceof Complex c2) {
            return Maja.rem(new Complex(d1), c2);
        } else if (a instanceof Complex c1 && b instanceof Double d2) {
            return Maja.rem(c1, new Complex(d2));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Long l2) {
            return dm1.map(d -> Maja.rem(d, l2));
        } else if (a instanceof Long l1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.rem(l1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
            return dm1.zipWith(dm2, Maja::rem);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> Maja.rem(d, d2));
        } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.rem(d1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Complex c2) {
            return ComplexMatrix.into(dm1.retype(Complex::new).map(c -> Maja.rem(c, c2)));
        } else if (a instanceof Complex c1 && b instanceof DoubleMatrix dm2) {
            return ComplexMatrix.into(dm2.retype(Complex::new).map(c -> Maja.rem(c1, c)));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Long l2) {
            Complex cl2 = new Complex(l2);
            return dm1.map(d -> Maja.rem(d, cl2));
        } else if (a instanceof Long l1 && b instanceof ComplexMatrix dm2) {
            Complex cl1 = new Complex(l1);
            return dm2.map(d -> Maja.rem(cl1, d));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof ComplexMatrix dm2) {
            return dm1.zipWith(dm2, Maja::rem);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Double d2) {
            Complex cd2 = new Complex(d2);
            return dm1.map(d -> Maja.rem(d, cd2));
        } else if (a instanceof Double d1 && b instanceof ComplexMatrix dm2) {
            Complex cd1 = new Complex(d1);
            return dm2.map(d -> Maja.rem(cd1, d));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Complex c2) {
            return dm1.map(d -> Maja.rem(d, c2));
        } else if (a instanceof Complex c1 && b instanceof ComplexMatrix dm2) {
            return dm2.map(d -> Maja.rem(c1, d));
        } else {
            throw new RuntimeException("Invalid type for rem.");
        }
    }

    @Override
    public Object visitExprNeq(ExpressionParser.ExprNeqContext ctx) {
        Object a = visit(ctx.expression(0)), b = visit(ctx.expression(1));
        if (a instanceof Long l1 && b instanceof Long l2) {
            return Maja.ne(l1, l2) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Double d2) {
            return Maja.ne(d1, d2) ? 1 : 0;
        } else if (a instanceof Complex c1 && b instanceof Complex c2) {
            return Maja.ne(c1, c2) ? 1 : 0;
        } else if (a instanceof Long l1 && b instanceof Double d2) {
            return Maja.ne(l1, d2) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Long l2) {
            return Maja.ne(d1, l2) ? 1 : 0;
        } else if (a instanceof Long l1 && b instanceof Complex c2) {
            return Maja.ne(new Complex(l1), c2) ? 1 : 0;
        } else if (a instanceof Complex c1 && b instanceof Long l2) {
            return Maja.ne(c1, new Complex(l2)) ? 1 : 0;
        } else if (a instanceof Double d1 && b instanceof Complex c2) {
            return Maja.ne(new Complex(d1), c2) ? 1 : 0;
        } else if (a instanceof Complex c1 && b instanceof Double d2) {
            return Maja.ne(c1, new Complex(d2)) ? 1 : 0;
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Long l2) {
            return dm1.map(d -> Maja.ne(d, l2) ? 1.0 : 0.0);
        } else if (a instanceof Long l1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.ne(l1, d) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
            return dm1.zipWith(dm2, (x, y) -> Maja.ne(x, y) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> Maja.ne(d, d2) ? 1.0 : 0.0);
        } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.ne(d1, d) ? 1.0 : 0.0);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Complex c2) {
            return dm1.map(d -> Maja.ne(c2, d) ? 1.0 : 0.0);
        } else if (a instanceof Complex c1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.ne(c1, d) ? 1.0 : 0.0);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Long l2) {
            Complex cl2 = new Complex(l2);
            return dm1.map(d -> Maja.ne(d, cl2) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof Long l1 && b instanceof ComplexMatrix dm2) {
            Complex cl1 = new Complex(l1);
            return dm2.map(d -> Maja.ne(cl1, d) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof ComplexMatrix dm2) {
            return dm1.zipWith(dm2, (x, y) -> Maja.ne(x, y) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Double d2) {
            Complex cd2 = new Complex(d2);
            return dm1.map(d -> Maja.ne(d, cd2) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof Double d1 && b instanceof ComplexMatrix dm2) {
            Complex cd1 = new Complex(d1);
            return dm2.map(d -> Maja.ne(cd1, d) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Complex c2) {
            return dm1.map(d -> Maja.ne(d, c2) ? Complex.ONE : Complex.ZERO);
        } else if (a instanceof Complex c1 && b instanceof ComplexMatrix dm2) {
            return dm2.map(d -> Maja.ne(c1, d) ? Complex.ONE : Complex.ZERO);
        } else {
            throw new RuntimeException("Invalid type for !=.");
        }
    }

    @Override
    public Object visitExprReal(ExpressionParser.ExprRealContext ctx) {
        return Double.parseDouble(ctx.REAL().getText());
    }

    @Override
    public Object visitExprVariable(ExpressionParser.ExprVariableContext ctx) {
        return env.get(ctx.ID().getText());
    }

    @Override
    public Object visitExprAdd(ExpressionParser.ExprAddContext ctx) {
        Object a = visit(ctx.expression(0)), b = visit(ctx.expression(1));
        if (a instanceof Long l1 && b instanceof Long l2) {
            return Maja.add(l1, l2);
        } else if (a instanceof Double d1 && b instanceof Double d2) {
            return Maja.add(d1, d2);
        } else if (a instanceof Complex c1 && b instanceof Complex c2) {
            return Maja.add(c1, c2);
        } else if (a instanceof Long l1 && b instanceof Double d2) {
            return Maja.add(l1, d2);
        } else if (a instanceof Double d1 && b instanceof Long l2) {
            return Maja.add(d1, l2);
        } else if (a instanceof Long l1 && b instanceof Complex c2) {
            return Maja.add(new Complex(l1), c2);
        } else if (a instanceof Complex c1 && b instanceof Long l2) {
            return Maja.add(c1, new Complex(l2));
        } else if (a instanceof Double d1 && b instanceof Complex c2) {
            return Maja.add(new Complex(d1), c2);
        } else if (a instanceof Complex c1 && b instanceof Double d2) {
            return Maja.add(c1, new Complex(d2));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Long l2) {
            return dm1.map(d -> Maja.add(d, l2));
        } else if (a instanceof Long l1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.add(l1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof DoubleMatrix dm2) {
            return dm1.zipWith(dm2, Maja::add);
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Double d2) {
            return dm1.map(d -> Maja.add(d, d2));
        } else if (a instanceof Double d1 && b instanceof DoubleMatrix dm2) {
            return dm2.map(d -> Maja.add(d1, d));
        } else if (a instanceof DoubleMatrix dm1 && b instanceof Complex c2) {
            return ComplexMatrix.into(dm1.retype(Complex::new).map(c -> Maja.add(c, c2)));
        } else if (a instanceof Complex c1 && b instanceof DoubleMatrix dm2) {
            return ComplexMatrix.into(dm2.retype(Complex::new).map(c -> Maja.add(c1, c)));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Long l2) {
            Complex cl2 = new Complex(l2);
            return dm1.map(d -> Maja.add(d, cl2));
        } else if (a instanceof Long l1 && b instanceof ComplexMatrix dm2) {
            Complex cl1 = new Complex(l1);
            return dm2.map(d -> Maja.add(cl1, d));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof ComplexMatrix dm2) {
            return dm1.zipWith(dm2, Maja::add);
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Double d2) {
            Complex cd2 = new Complex(d2);
            return dm1.map(d -> Maja.add(d, cd2));
        } else if (a instanceof Double d1 && b instanceof ComplexMatrix dm2) {
            Complex cd1 = new Complex(d1);
            return dm2.map(d -> Maja.add(cd1, d));
        } else if (a instanceof ComplexMatrix dm1 && b instanceof Complex c2) {
            return dm1.map(d -> Maja.add(d, c2));
        } else if (a instanceof Complex c1 && b instanceof ComplexMatrix dm2) {
            return dm2.map(d -> Maja.add(c1, d));
        } else {
            throw new RuntimeException("Invalid type for +.");
        }
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private static Object simplify(Object input, boolean recurrent) {
        // As a rule of thumb:
        //  - If we have a generic matrix that is not a double or a complex matrix, see if we can cast it.
        //  - If we have a complex value, if it has im=0 then cast it to a double, otherwise leave it as a complex.
        //  - If we have a double value, if it has no fractional part then cast it to a long, otherwise leave it as a double.
        //  - If we have a long value, leave it as a long.

        // This generally breaks cases like log(-5). These functions need to be adjusted to call the complex overload in case
        // the function returns a value that is not normal (NaN, +Infinity, -Infinity) or throws an exception.

        if (!(input instanceof DoubleMatrix) && !(input instanceof ComplexMatrix) && input instanceof Matrix mat) {
            if(mat.ravel().stream().allMatch(e -> e instanceof Double)) {
                return DoubleMatrix.into(mat);
            } else if(mat.ravel().stream().allMatch(e -> e instanceof Complex)) {
                return DoubleMatrix.into(mat);
            } else if(mat.ravel().stream().anyMatch(e -> e instanceof Complex)
                    && mat.ravel().stream().allMatch(e -> e instanceof Complex || e instanceof Double || e instanceof Long)) {
                return ComplexMatrix.into(mat.retype(x -> {
                    if (x instanceof Complex) return x;
                    else if (x instanceof Double d) return new Complex(d);
                    else return new Complex(((Long) x).doubleValue());
                }));
            } else if(mat.ravel().stream().anyMatch(e -> e instanceof Double)
                    && mat.ravel().stream().allMatch(e -> e instanceof Double || e instanceof Long)) {
                return DoubleMatrix.into(mat.retype(x -> {
                    if (x instanceof Double) return x;
                    else return ((Long) x).doubleValue();
                }));
            } else {
                if(!recurrent)
                    return simplify(mat.map(x -> simplify(x, false)), true);
                else
                    return mat;
            }
        } else if (input instanceof Complex c) {
            if (c.im() == 0) {
                return c.re();
            } else {
                return c;
            }
        } else if (input instanceof Double d) {
            if (d == Math.floor(d)) {
                return (long) d.doubleValue();
            } else {
                return d;
            }
        } else {
            return input;
        }
    }

    @Override
    public Object visit(ParseTree tree) {
        return simplify(tree.accept(this), false);
    }
}
