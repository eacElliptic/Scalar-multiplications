package fr.android.bri.simple_glv;

import android.app.Activity;
import android.os.Bundle;
import android.os.Debug;
import android.view.View;
import android.widget.SeekBar;
import android.widget.TextView;
import android.widget.Toast;

import java.math.BigInteger;
import java.security.SecureRandom;

public class MainActivity extends Activity implements SeekBar.OnSeekBarChangeListener {

    private SeekBar bar;
    private TextView iteration, state;
    private int nbiter = 1;

    private AffPoint[] PP;
    private BigInteger a, b, p, phi, efp, refp, beta, aa, bb, Na;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        bar = (SeekBar) findViewById(R.id.seekBar1);
        bar.setOnSeekBarChangeListener(this);
        iteration = (TextView) findViewById(R.id.nbiter);
        state = (TextView) findViewById(R.id.state);

        init();
    }

    @Override
    public void onProgressChanged(SeekBar seekBar, int progress, boolean fromUser) {
        nbiter = (1 << progress);
        iteration.setText(String.valueOf(progress));
    }



    public void init() {
        // Note: N, T and c are always equal to ONE, so we don't include them in the code.
        BigInteger q, X, Y;

        PP = new AffPoint[4];
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
		/* integer part of square root of efp */
        refp = new BigInteger("340282366920938463463374607431768211455");
		/* root of  phi^2+phi+1 = 0 mod #E(Fp) */
        phi = new BigInteger("86115113571596370384202877098459685148505633832452548870570058599719872579953");

        initGLV(phi, efp, refp);

        if (PP[0].isOnCurve(a, b))
            Toast.makeText(this, "P OK", Toast.LENGTH_SHORT).show();
        else
            Toast.makeText(this, "P not OK", Toast.LENGTH_SHORT).show();
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

    public void go(View v) {
        int i;
        GLVScalar k;
        BigInteger ZZ;

        SecureRandom rg = new SecureRandom();
        Point Q = new Point(p);

        PP[1] = new AffPoint(p);
        PP[3] = new AffPoint(p);

        long at, bt, diff1, diff2;

        diff1=0;
        diff2=0;

        for (i = 0; i < nbiter; i++) {
            PP[2] = PP[0].endo(beta); /* Point PHI_P */
            PP[3].X = PP[2].X;
            PP[3].Y = PP[2].Y;
            PP[3].ADD(PP[0]); /* Point P+PHI_P */
            PP[1].X = PP[0].X;
            PP[1].Y = PP[0].Y.negate();
            PP[1].ADD(PP[2]); /* Point -P+PHI_P */

            at = System.currentTimeMillis();
            k = new GLVScalar(256, aa, bb, Na, efp, rg);
            bt = System.currentTimeMillis();
            diff1+=(bt-at);

            at = System.currentTimeMillis();
            Q = Point.PointFromGLVScalar(k, PP);
            bt = System.currentTimeMillis();
            diff2+=(bt-at);

            ZZ = Q.Z.modInverse(p);
            /* P is transformed in affine coordinates for next iteration */
            PP[0].X = Q.X.multiply(ZZ).multiply(ZZ).mod(p);
            PP[0].Y = Q.Y.multiply(ZZ).multiply(ZZ).multiply(ZZ).mod(p);
        }

        if (Q.isOnCurve(a, b))
            state.setText("Completed, MULT: kP OK\nTime in milliseconds: "  + diff1 + "," + diff2 + "," + (diff1+diff2));
        else
            state.setText("Completed, MULT: kP not OK");
    }



    @Override
    public void onStartTrackingTouch(SeekBar seekBar) {
    }

    @Override
    public void onStopTrackingTouch(SeekBar seekBar) {
    }

}
