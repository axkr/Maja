package rocks.palaiologos.maja;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

class PrimeConstantPool {
    public static short[] bases_2;
    public static long[] bases;

    static {
        // Read the lines of the file "maja_bases.txt" from classpath
        // and parse them into the arrays bases_2 and bases.
        try {
            InputStream in = PrimeConstantPool.class.getClassLoader().getResourceAsStream("maja_bases.txt");
            byte[] data = new byte[1024 * 1024];
            int n = in.read(data);
            String[] lines = new String(data, 0, n, StandardCharsets.US_ASCII).split(System.lineSeparator());
            String[] bases_2_str = lines[0].split(",");
            short[] bases_2 = new short[bases_2_str.length];
            for (int i = 0; i < bases_2_str.length; ++i) {
                bases_2[i] = Short.parseShort(bases_2_str[i]);
            }
            PrimeConstantPool.bases_2 = bases_2;
            String[] bases_str = lines[1].split(",");
            long[] bases = new long[bases_str.length];
            for (int i = 0; i < bases_str.length; ++i) {
                bases[i] = Long.parseUnsignedLong(bases_str[i]);
            }
            PrimeConstantPool.bases = bases;
            in.close();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
}
