package fr.univtln.eac;

import java.math.BigInteger;
import java.security.SecureRandom;

import static java.math.BigInteger.ONE;
import static java.math.BigInteger.ZERO;

public class EACTiming {

    private BigInteger beta, a, b;
    private Point P;

    public EACTiming() {
        init();
    }

    private void init() {
        BigInteger p, q, X, Y;

        p = ZERO;
        p = p.setBit(358);        
        q = new BigInteger("36855");
        p = p.subtract(q); /* p = 2^358-36855 */
        /* Curve is y^2 = x^3 + 17 */
        a = new BigInteger("0");
        b = new BigInteger("17");
        /* 2 is a generator of Fp and beta is an element of order 3 */
        beta = new BigInteger("2");
        beta = beta.modPow(p.subtract(ONE).divide(new BigInteger("3")), p);

        X = new BigInteger("568554169514108215082878628903438431409037009413741676480372717420637861491120753156713065337030704866093074");
        Y = new BigInteger("125804816918392729528454862544298396901812597725830867877741173065843595271902914550421712710898294754082775");
        /* P represents point P itself and the point Phi_P */
        /* see Point.java for class definition */
        P = new Point(X, Y, p);

        if (P.isOnCurve(0, a, b)) 
            System.out.println("\nP OK");
        else  
            System.out.println("\nP not OK");

    }

    public void go(int nbiter) {
        int i, j, len;
        byte[] eac;
        BigInteger ZZ;
        SecureRandom rg = new SecureRandom();

        len = 256;
        eac = new byte[len];
        
        System.out.println("Running...");
        System.out.println("EAC 358 bits benchmark");

        long at, bt, diff1;

        diff1 = 0;
        
        for (i = 0; i < nbiter; i++) {
            P.Z = ONE;
            /* we compute PHI_P */
            P.PQ[1] = P.PQ[0].multiply(beta).mod(Point.p);
            P.PQ[3] = P.PQ[2];

            for (j = 0; j < len; j++) 
                eac[j] = (byte) rg.nextInt(2);
            
            at = System.currentTimeMillis();
            /* coordinates of P are updated so that at the end of the method, */
            /* P contains ((k-t)P, kP) for some integer t */
            P.PointFromEAC(eac);
            
            bt = System.currentTimeMillis();
            diff1 += (bt - at);

            ZZ = P.Z.modInverse(Point.p);
            /* P is transformed in affine coordinates for next iteration */
            P.PQ[0] = P.PQ[0].multiply(ZZ).multiply(ZZ).mod(Point.p);
            P.PQ[2] = P.PQ[2].multiply(ZZ).multiply(ZZ).multiply(ZZ).mod(Point.p);
        }
        
        if (P.isOnCurve(1, a, b)) {
            System.out.println("\nMULT: kP OK");
            System.out.println("Time in milliseconds: " + diff1);
        }
        else {
            System.out.println("\nMULT: kP not OK");
            System.out.println("Time in milliseconds: " + diff1);
        }
    }

    public boolean testpoint(Point Q_, BigInteger a, BigInteger b) {
        if (Q_.isOnCurve(1, a, b)) {
            System.out.println("kP is on the curve !\n\n");
            System.out.println(Q_.toString(1, 16));
            return true;
        }
        System.out.println("Error : kP is not on the curve !");
        return false;
    }

}

