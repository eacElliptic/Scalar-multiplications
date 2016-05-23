package fr.univtln.sac_glv;

import java.math.BigInteger;
import java.security.SecureRandom;

public class SGLVTiming {

    private int l;
    private AffPoint[] PP;
    private BigInteger a, b, p, phi, efp, refp, beta, aa, bb, Na;


    public SGLVTiming() {
        init();
    }


    private void init() {
        BigInteger q, X, Y;

        PP = new AffPoint[2];
        p = BigInteger.ZERO;
        p = p.setBit(256); /* p = 2^256 */
        q = new BigInteger("1539");
        p = p.subtract(q); /* p = 2^256-1539 */
        /* Curve is y^2 = x^3 + 5 */
        a = BigInteger.ZERO;
        b = new BigInteger("5");
        X = new BigInteger("66043340678279369258981193985450715448372246692524118399919799979198175326756");
        Y = new BigInteger("52931614173969837860927627620185828149756459753631369661231760168280520196499");
        PP[0] = new AffPoint(X, Y, p);  /* Point P */
        /* 5 is a generator of Fp and beta is an element of order 3 */
        beta = new BigInteger("5");
        beta = beta.modPow(p.subtract(BigInteger.ONE).divide(new BigInteger("3")), p);
        /* order of P */
        efp = new BigInteger("115792089237316195423570985008687907852920869663551026044856920202296625078603");
        /* ceil part of log_2 (efp) /2 +1  for GLV-SAC */
        l = 129;
        /* integer part of square root of efp */
        refp = new BigInteger("340282366920938463463374607431768211455");
        /* root of  phi^2+phi+1 = 0 mod #E(Fp) */
        phi = new BigInteger("86115113571596370384202877098459685148505633832452548870570058599719872579953");

        initGLV(phi, efp, refp);
        
        if (PP[0].isOnCurve(a, b)) 
            System.out.println("\nP OK");
        else 
            System.out.println("\nP not OK");
        
    }


    public void initGLV(BigInteger lambda, BigInteger n, BigInteger rootn) {
        BigInteger u, v, x1, x2, y1, y2, q, r, x, y;
        u = n;
        v = lambda;

        y = BigInteger.ZERO;
        x1 = BigInteger.ONE;
        y1 = BigInteger.ZERO;
        x2 = BigInteger.ZERO;
        y2 = BigInteger.ONE;
        while (u.compareTo(rootn) == 1) {
            q = v.divide(u);
            r = v.subtract(q.multiply(u));
            x = x2.subtract(q.multiply(x1));
            y = y2.subtract(q.multiply(y1));
            v = u;
            u = r;
            x2 = x1;
            x1 = x;
            y2 = y1;
            y1 = y;
        }
        aa = u;
        bb = y.negate();
        Na = aa.multiply(aa).add(bb.multiply(bb)).subtract(aa.multiply(bb));
        aa = aa.subtract(bb);
    }


    public void go(int nbiter) {

        int i;
        SGLVScalar k;
        BigInteger ZZ;
        Point Q = new Point(p);

        SecureRandom rg = new SecureRandom();
        
        System.out.println("Running...");
        System.out.println("SGLV 256 bits benchmark");

        PP[1] = new AffPoint(p);
        
        long at, bt, diff1, diff2;

        diff1 = 0;
        diff2 = 0;

        for (i = 0; i < nbiter; i++) {
            PP[1].X = PP[0].X.multiply(beta).mod(p); /* Point PHI_P */
            PP[1].Y = PP[0].Y;
            PP[1].ADD(PP[0]); /* Point P+PHI_P */

            at = System.currentTimeMillis();
            
            k = new SGLVScalar(256, aa, bb, Na, efp, rg, l); 
            
            bt = System.currentTimeMillis();
            diff1 += (bt - at);
            
            at = System.currentTimeMillis();
            
            Q = Point.PointFromGLVScalar(k, PP);
            
            bt = System.currentTimeMillis();
            diff2 += (bt - at);
            
            ZZ = Q.Z.modInverse(p);
            /* P is transformed in affine coordinates for next iteration */
            PP[0].X = Q.X.multiply(ZZ).multiply(ZZ).mod(p);
            PP[0].Y = Q.Y.multiply(ZZ).multiply(ZZ).multiply(ZZ).mod(p);
        }
        
        if (Q.isOnCurve(a, b)){ 
            System.out.println("\nMULT: kP OK");
            System.out.println("Time in milliseconds: " + diff1 + "," + diff2 + "," + (diff1+diff2));
        }
        else 
            System.out.println("\nMULT: kP not OK");
        
    }


    public boolean testpoint(Point Q_, BigInteger a, BigInteger b) {
        if (Q_.isOnCurve(a, b)) {
            System.out.println("kP is on the curve !\n\n");
            return true;
        }
        System.out.println("Error : kP is not on the curve !");
        return false;
    }

}
