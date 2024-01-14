package rocks.palaiologos.maja;

class Prime {
    // Fast Primality Testing for Integers That Fit into a Machine Word
    // Michal Forisek and Jakub Jancina
    // All codes published under https://people.ksp.sk/~misof/primes/
    // are available under the CC BY-NC 4.0 Int'l license.
    // [https://creativecommons.org/licenses/by-nc/4.0/legalcode]

    public static long pow(long a, long e, long n) {
        long ans = 1;
        while (e > 0) {
            if (e % 2 == 1) ans = (ans * a) % n;
            a = (a * a) % n;
            e >>= 1;
        }
        return ans;
    }
    
    public static boolean is_SPRP(long a, long n) {
        int d = 0;
        long t = n - 1;
        while (t % 2 == 0) {
            ++d;
            t >>= 1;
        }
        long x = pow(a, t, n);
        if (x == 1) return true;
        for (int i = 0; i < d; ++i) {
            if (x == n - 1) return true;
            x = (x * x) % n;
            if (x == 1) return false;
        }
        return false;
    }
    public static boolean is_prime_1b(int x) {
        if (x < 2) return false;
        if (x == 2 || x == 3 || x == 5 || x == 7) return true;
        if (x % 2 == 0 || x % 3 == 0 || x % 5 == 0 || x % 7 == 0) return false;
        if (x < 121) return true;
        return is_SPRP(PrimeConstantPool.bases_2[(int) (((0xAFF7B4L * x) >> 7) & 1023)], x);
    }

    static class MontgomeryResult {
        public long u, v;
        public MontgomeryResult(long u, long v) {
            this.u = u; this.v = v;
        }
    }

    public static void xbinGCD(long a, long b, MontgomeryResult r)  {
        long alpha, beta, u, v;
        u = 1; v = 0;
        alpha = a; beta = b;

        while (Long.compareUnsigned(a, 0) > 0) {
            a = a >>> 1;
            if ((u & 1) == 0) {
                u = u >>> 1; v = v >>> 1;
            } else {
                u = ((u ^ beta) >>> 1) + (u & beta);
                v = (v >>> 1) + alpha;
            }
        }

        r.u = u; r.v = v;
    }


    public static long modul64(long x, long y, long z) {
        long i, t;

        for (i = 1; i <= 64; i++) {
            t = x >>> 63;
            x = (x << 1) | (y >>> 63);
            y = y << 1;
            if (Long.compareUnsigned(x | t, z) >= 0) {
                x = x - z;
                y = y + 1;
            }
        }
        return x;
    }

    public static void mulul64(long u, long v, MontgomeryResult r) {
        long u0, u1, v0, v1, k, t;
        long w0, w1, w2;

        u1 = u >>> 32; u0 = u & 0xFFFFFFFFL;
        v1 = v >>> 32; v0 = v & 0xFFFFFFFFL;

        t = u0*v0;
        w0 = t & 0xFFFFFFFFL;
        k = t >>> 32;

        t = u1*v0 + k;
        w1 = t & 0xFFFFFFFFL;
        w2 = t >>> 32;

        t = u0*v1 + w1;
        k = t >>> 32;

        long wlo = (t << 32) + w0;
        long whi = u1*v1 + w2 + k;

        r.u = wlo; r.v = whi;
    }

    public static long montmul(long abar, long bbar, long m, long mprime, MontgomeryResult r) {
        long thi, tlo, tm, tmmhi, tmmlo, uhi, ulo; boolean ov;

        mulul64(abar, bbar, r); thi = r.v; tlo = r.u;

        tm = tlo*mprime;

        mulul64(tm, m, r); tmmhi = r.v; tmmlo = r.u;

        ulo = tlo + tmmlo;
        uhi = thi + tmmhi;
        if (Long.compareUnsigned(ulo, tlo) < 0) uhi = uhi + 1;

        ov = (Long.compareUnsigned(uhi, thi) < 0) || ((Long.compareUnsigned(uhi, thi) == 0) && (Long.compareUnsigned(ulo, tlo) < 0));

        ulo = uhi;
        uhi = 0;

        if (ov || Long.compareUnsigned(ulo, m) >= 0)
            ulo = ulo - m;

        return ulo;
    }

    public static long mulmodMont(long baseM,long e,long modul,long pv,long oneM, MontgomeryResult r) {
        long ans = oneM;
        while(Long.compareUnsigned(e, 0) > 0) {
            if((e&1) == 1)
                ans = montmul(baseM,ans,modul,pv,r);
            baseM = montmul(baseM,baseM,modul,pv,r);
            e>>>=1;
        }
        return ans;
    }

    public static boolean is_SPRP(long base,long modul, MontgomeryResult r) {
        if(Long.compareUnsigned(base,modul)>=0) base=Long.remainderUnsigned(base, modul);
        long pu,pv;
        xbinGCD(1l<<63l,modul,r);
        pu = r.u; pv = r.v;
        long baseM = modul64(base,0,modul);
        long oneM = modul64(1,0,modul);
        long moneM = modul - oneM;
        long e = modul-1;
        while(Long.compareUnsigned(e&1,0)==0) e>>>=1;
        long t = mulmodMont(baseM,e,modul,pv,oneM,r);
        if(t==oneM) return true;
        while(Long.compareUnsigned(e,modul-1)<0) {
            if(Long.compareUnsigned(t,moneM)==0) return true;
            t = montmul(t,t,modul,pv,r);
            if(Long.compareUnsigned(t,oneM)==0) return false;
            e<<=1;
        }
        return false;
    }

    public static int hashh(long x) {
        x = ((x >>> 32) ^ x) * 0x45d9f3b3335b369l;
        x = ((x >>> 32) ^ x) * 0x3335b36945d9f3bl;
        x = ((x >>> 32) ^ x);
        return (int) (x & 262143l);
    }

    public static boolean is_prime_2_64(long a) {
        MontgomeryResult r = new MontgomeryResult(0, 0);
        if (Long.compareUnsigned(a,2)==0 || Long.compareUnsigned(a,3)==0 || Long.compareUnsigned(a,5)==0 || Long.compareUnsigned(a,7)==0) return true;
        if (Long.remainderUnsigned(a,2)==0 || Long.remainderUnsigned(a,3)==0 || Long.remainderUnsigned(a,5)==0 || Long.remainderUnsigned(a,7)==0) return false;
        if (Long.compareUnsigned(a,121)<0) return Long.compareUnsigned(a,1) > 0;
        return is_SPRP(2,a,r) && is_SPRP(PrimeConstantPool.bases[hashh(a)],a,r);
    }
}
