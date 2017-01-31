package com.example.fang.sac_glv;

import android.app.Activity;
import android.os.Bundle;
import android.os.Debug;
import android.view.View;
import android.widget.SeekBar;
import android.widget.TextView;
import android.widget.Toast;

import java.math.BigInteger;
import java.security.SecureRandom;

public class MainActivity extends Activity implements SeekBar.OnSeekBarChangeListener  {

    private SeekBar bar;
    private TextView iteration, state;
    private int nbiter = 1;

    private ExtAffPoint[] PP;
    private BigInteger curve_a, curve_d, p, phi, efp, refp, beta, aa, bb, Na;



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
        BigInteger q, X, Y;

        PP = new ExtAffPoint[2];
        p = BigInteger.ZERO;
        p = p.setBit(256);
        q = new BigInteger("43443");
        p = p.subtract(q); /* p = 2^256-43443 */

//        Curve is : -(x^2) + y^2 = 1 + (x^2)*(y^2)
        curve_a = new BigInteger("-1");
        curve_d = BigInteger.ONE;

        X = new BigInteger("22174792725782664025844270156401847138026886831119649840138701280259288954423");
        Y = new BigInteger("87737099222887081277295330229931648697203685991117339269318803967832062780143");
        PP[0] = new ExtAffPoint(X, Y, p);  /* ExtAffPoint P */

//        beta is an element of order 4
        beta = new BigInteger("9058966984510276008007633023999327412385577717625298889713350828587238576126");
		/* order of P */
        efp = new BigInteger("14474011154664524427946373126085988481582509411733023890661985475889632981401");
		/* integer part of square root of efp */
        refp = new BigInteger("120307984584002255772516886238812528463");
		/* root of  phi^2 + 1 = 0 mod #E(Fp) */
        phi = new BigInteger("10749260569431236026102217475317958236773166001345671471970756507155106163337");

        initGLV(phi, efp, refp);

        if (PP[0].isOnCurve(curve_a, curve_d))
            Toast.makeText(this, "P OK", Toast.LENGTH_LONG).show();
        else
            Toast.makeText(this, "P non OK", Toast.LENGTH_LONG).show();
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
        Na = aa.multiply(aa).add(bb.multiply(bb));
    }

    public void go(View v) {
        int i, scal_size, scal_l;
        SGLVScalar k;
        BigInteger ZZ, k3;

        scal_size = 253;    // size of base point order.
        scal_l = 128;       // (ceil part of log_2 (efp)/2) + 1
        SecureRandom rg = new SecureRandom();
        ExtProjPoint Q = new ExtProjPoint(p);


        long at, bt, diff1, diff2;

        diff1=0;
        diff2=0;


        //        Debug.startMethodTracing("ted_glvtrace");
        Debug.startMethodTracing("sac_glv_ted", 1950*1024*1024);

        for (i = 0; i < nbiter; i++) {

            PP[1] = PP[0].endo(beta); /* ExtAffPoint PHI_P */
            PP[1].ADD(PP[0], curve_a, curve_d); /* ExtAffPoint P+PHI_P */



            k3 = (new BigInteger(scal_size, rg)).mod(efp);


//            at = System.currentTimeMillis();

            k = new SGLVScalar(k3, aa, bb, Na, scal_l);

//            bt = System.currentTimeMillis();
//            diff1+=(bt-at);



//            at = System.currentTimeMillis();

            Q = ExtProjPoint.PointFromGLVScalar(k, PP);

//            bt = System.currentTimeMillis();
//            diff2+=(bt-at);



            ZZ = Q.Z.modInverse(p);
            /* P is transformed in affine coordinates for next iteration */
            PP[0].X = Q.X.multiply(ZZ).mod(p);
            PP[0].Y = Q.Y.multiply(ZZ).mod(p);
            PP[0].T = PP[0].X.multiply(PP[0].Y).mod(p);
        }


        Debug.stopMethodTracing();

        if (Q.isOnCurve(curve_a, curve_d))
            state.setText("MULT : kP OK\n\nTime in milliseconds : "  + diff1 + " , " + diff2 + " , " + (diff1+diff2));
        else
            state.setText(R.string.glv_com_ko);
    }



    @Override
    public void onStartTrackingTouch(SeekBar seekBar) {
    }

    @Override
    public void onStopTrackingTouch(SeekBar seekBar) {
    }

}
