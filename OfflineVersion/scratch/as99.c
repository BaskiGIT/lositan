#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* START OF AS 99 */
/* Subroutine */ int jnsn( xbar,  sd, rb1,  bb2, 
	 itype,  gamma,  delta,  xlam,  xi, 
	 ifault)
	 float  *xbar,  *sd,  *rb1,  *bb2, *gamma,  *delta,  *xlam,
	  *xi;
	  int *itype, *ifault;
{
    /* Initialized data */

    static float tol = .01f;
    static float zero = 0.f;
    static float quart = .25f;
    static float half = .5f;
    static float one = 1.f;
    static float two = 2.f;
    static float three = 3.f;
    static float four = 4.f;

    /* System generated locals */
    float r__1;
    double d__1, d__2;

    /* Builtin functions */

    /* Local variables */
    float u, w, x, y;
    int sbfit();
    int fault;
    float b1, b2;
    int sufit();


/*        ALGORITHM AS 99  APPL. STATIST. (1976) VOL.25, P.180 */

/*        FINDS TYPE AND PARAMETERS OF A JOHNSON CURVE */
/*        WITH GIVEN FIRST FOUR MOMENTS */




    *ifault = 1;
    if (*sd < zero) {
	return 0;
    }
    *ifault = 0;
    *xi = zero;
    *xlam = zero;
    *gamma = zero;
    *delta = zero;
    if (*sd > zero) {
	goto L10;
    }
    *itype = 5;
    *xi = *xbar;
    return 0;
L10:
    b1 = *rb1 * *rb1;
    b2 = *bb2;
    fault = 0;

/*        TEST WHETHER LOGNORMAL (OR NORMAL) REQUESTED */

    if (b2 >= zero) {
	goto L30;
    }
L20:
    if (fabs(*rb1) <= tol) {
	goto L70;
    }
    goto L80;

/*        TEST FOR POSITION RELATIVE TO BOUNDARY LINE */

L30:
    if (b2 > b1 + tol + one) {
	goto L60;
    }
    if (b2 < b1 + one) {
	goto L50;
    }

/*        ST DISTRIBUTION */

L40:
    *itype = 5;
    r__1 = one - four / (b1 + four);
    y = half + half * sqrt(r__1);
    if (*rb1 > zero) {
	y = one - y;
    }
    r__1 = y * (one - y);
    x = *sd / sqrt(r__1);
    *xi = *xbar - y * x;
    *xlam = *xi + x;
    *delta = y;
    return 0;
L50:
    *ifault = 2;
    return 0;
L60:
    r__1 = b2 - three;
    if (fabs(*rb1) > tol || fabs(r__1) > tol) {
	goto L80;
    }

/*        NORMAL DISTRIBUTION */

L70:
    *itype = 4;
    *delta = one / *sd;
    *gamma = -(*xbar) / *sd;
    return 0;

/*        TEST FOR POSITION RELATIVE TO LOGNORMAL LINE */

L80:
    x = half * b1 + one;
    r__1 = quart * b1 + one;
    y = fabs(*rb1) * sqrt(r__1);
    d__1 = (double) (x + y);
    d__2 = (double) (one / three);
    u = pow(d__1, d__2);
    w = u + one / u - one;
    u = w * w * (three + w * (two + w)) - three;
    if (b2 < zero || fault) {
	b2 = u;
    }
    x = u - b2;
    if (fabs(x) > tol) {
	goto L90;
    }

/*        LOGNORMAL (SL) DISTRIBUTION */

    *itype = 1;
    *xlam = (*rb1 < 0.0 ? -one : one);
    u = *xlam * *xbar;
    r__1 = log(w);
    x = one / sqrt(r__1);
    *delta = x;
    r__1 = w * (w - one) / (*sd * *sd);
    y = half * x * log(r__1);
    *gamma = y;
    r__1 = (half / x - y) / x;
    *xi = *xlam * (u - exp(r__1));
    return 0;

/*        SB OR SU DISTRIBUTION */

L90:
    if (x > zero) {
	goto L100;
    }
    *itype = 2;
    sufit(xbar, sd, rb1, &b2, gamma, delta, xlam, xi);
    return 0;
L100:
    *itype = 3;
    sbfit(xbar, sd, rb1, &b2, gamma, delta, xlam, xi, &fault);
    if (! fault) {
	return 0;
    }

/*        FAILURE - TRY TO FIT APPROXIMATE RESULT */

    *ifault = 3;
    if (b2 > b1 + two) {
	goto L20;
    }
    goto L40;
} /* jnsn_ */


/* Subroutine */ int sufit(xbar, sd, rb1, b2,gamma, delta, xlam, xi)
	float *xbar, *sd, *rb1, *b2, *gamma, *delta, *xlam, *xi;
{
    /* Initialized data */

    static float tol = .01f;
    static float zero = 0.f;
    static float one = 1.f;
    static float two = 2.f;
    static float three = 3.f;
    static float four = 4.f;
    static float six = 6.f;
    static float seven = 7.f;
    static float eight = 8.f;
    static float nine = 9.f;
    static float ten = 10.f;
    static float sixten = 16.f;
    static float half = .5f;
    static float one5 = 1.5f;
    static float two8 = 2.8f;

    /* System generated locals */
    float r__1, r__2, r__3;

    /* Builtin functions */

    /* Local variables */
    float a, b, v, w, x, y, z, b1, b3, w1, wm1;


/*        ALGORITHM AS 99.1  APPL. STATIST. (1976) VOL.25, P.180 */

/*        FINDS PARAMETERS OF JOHNSON SU CURVE WITH */
/*        GIVEN FIRST FOUR MOMENTS */




    b1 = *rb1 * *rb1;
    b3 = *b2 - three;

/*        W IS FIRST ESTIMATE OF EXP(DELTA ** (-2)) */

    r__1 = two * *b2 - two8 * b1 - two;
    w = sqrt(r__1);
    r__1 = w - one;
    w = sqrt(r__1);
    if (fabs(*rb1) > tol) {
	goto L10;
    }

/*        SYMMETRICAL CASE - RESULTS ARE KNOWN */

    y = zero;
    goto L20;

/*        JOHNSON ITERATION (USING Y FOR HIS M) */

L10:
    w1 = w + one;
    wm1 = w - one;
    z = w1 * b3;
    v = w * (six + w * (three + w));
    a = eight * (wm1 * (three + w * (seven + v)) - z);
    b = sixten * (wm1 * (six + v) - b3);
    r__1 = a * a - two * b * (wm1 * (three + w * (nine + w * (ten + v))) - 
	    two * w1 * z);
    y = (sqrt(r__1) - a) / b;
/* Computing 2nd power */
    r__1 = four * (w + two) * y + three * w1 * w1;
/* Computing 3rd power */
    r__2 = two * y + w1, r__3 = r__2;
    z = y * wm1 * (r__1 * r__1) / (two * (r__3 * (r__2 * r__2)));
    v = w * w;
    r__1 = one - two * (one5 - *b2 + b1 * (*b2 - one5 - v * (one + half * v)) 
	    / z);
    w = sqrt(r__1);
    r__1 = w - one;
    w = sqrt(r__1);
    r__1 = b1 - z;
    if (fabs(r__1) > tol) {
	goto L10;
    }

/*        END OF ITERATION */

    y /= w;
    r__1 = y + one;
    r__2 = sqrt(y) + sqrt(r__1);
    y = log(r__2);
    if (*rb1 > zero) {
	y = -y;
    }
L20:
    r__1 = one / log(w);
    x = sqrt(r__1);
    *delta = x;
    *gamma = y * x;
    y = exp(y);
    z = y * y;
    r__1 = half * (w - one) * (half * w * (z + one / z) + one);
    x = *sd / sqrt(r__1);
    *xlam = x;
    *xi = half * sqrt(w) * (y - one / y) * x + *xbar;
    return 0;
} /* sufit_ */

/* Subroutine */ int sbfit(xbar, sigma, rtb1, b2, gamma, delta, xlam, xi, fault)
	float *xbar, *sigma, *rtb1, *b2, *gamma, *delta, *xlam, *xi;
	int *fault;
{
    /* Initialized data */

    static float tt = 1e-4f;
    static float tol = .01f;
    static int limit = 50;
    static float zero = 0.f;
    static float one = 1.f;
    static float two = 2.f;
    static float three = 3.f;
    static float four = 4.f;
    static float six = 6.f;
    static float half = .5f;
    static float quart = .25f;
    static float one5 = 1.5f;
    static float a1 = .0124f;
    static float a2 = .0623f;
    static float a3 = .4043f;
    static float a4 = .408f;
    static float a5 = .479f;
    static float a6 = .485f;
    static float a7 = .5291f;
    static float a8 = .5955f;
    static float a9 = .626f;
    static float a10 = .64f;
    static float a11 = .7077f;
    static float a12 = .7466f;
    static float a13 = .8f;
    static float a14 = .9281f;
    static float a15 = 1.0614f;
    static float a16 = 1.25f;
    static float a17 = 1.7973f;
    static float a18 = 1.8f;
    static float a19 = 2.163f;
    static float a20 = 2.5f;
    static float a21 = 8.5245f;
    static float a22 = 11.346f;

    /* System generated locals */
    float r__1;
    double d__1, d__2, d__3, d__4;

    /* Builtin functions */

    /* Local variables */
    float rbet, d, e, f, g;
    int j, k, m;
    float s, t, u, w, x, y, deriv[4], b1, h2, h3, h4, dd[4], h2a, h2b, rb1;
    int neg;
    int moma();
    float hmu[6], bet2;


/*        ALGORITHM AS 99.2  APPL. STATIST. (1976) VOL.25, P.180 */

/*        FINDS PARAMETERS OF JOHNSON SB CURVE WITH */
/*        GIVEN FIRST FOUR MOMENTS */




    rb1 = fabs(*rtb1);
    b1 = rb1 * rb1;
    neg = *rtb1 < zero;

/*        GET D AS FIRST ESTIMATE OF DELTA */

    e = b1 + one;
    x = half * b1 + one;
    r__1 = quart * b1 + one;
    y = fabs(rb1) * sqrt(r__1);
    d__1 = (double) (x + y);
    d__2 = (double) (one / three);
    u = pow(d__1, d__2);
    w = u + one / u - one;
    f = w * w * (three + w * (two + w)) - three;
    e = (*b2 - e) / (f - e);
    if (fabs(rb1) > tol) {
	goto L5;
    }
    f = two;
    goto L20;
L5:
    r__1 = log(w);
    d = one / sqrt(r__1);
    if (d < a10) {
	goto L10;
    }
    f = two - a21 / (d * (d * (d - a19) + a22));
    goto L20;
L10:
    f = a16 * d;
L20:
    f = e * f + one;
    if (f < a18) {
	goto L25;
    }
    d__1 = (double) (three - f);
    d__2 = (double) (-a5);
    d = (a9 * f - a4) * pow(d__1, d__2);
    goto L30;
L25:
    d = a13 * (f - one);

/*        GET G AS FIRST ESTIMATE OF GAMMA */

L30:
    g = zero;
    if (b1 < tt) {
	goto L70;
    }
    if (d > one) {
	goto L40;
    }
    d__1 = (double) d;
    d__2 = (double) a17;
    d__3 = (double) b1;
    d__4 = (double) a6;
    g = (a12 * pow(d__1, d__2) + a8) * pow(d__3, d__4);
    goto L70;
L40:
    if (d <= a20) {
	goto L50;
    }
    u = a1;
    y = a7;
    goto L60;
L50:
    u = a2;
    y = a3;
L60:
    d__1 = (double) b1;
    d__2 = (double) (u * d + y);
    g = pow(d__1, d__2) * (a14 + d * (a15 * d - a11));
L70:
    m = 0;

/*        MAIN ITERATION STARTS HERE */

L80:
    ++m;
    *fault = m > limit;
    if (*fault) {
	return 0;
    }

/*        GET FIRST SIX MOMENTS FOR LATEST G AND D VALUES */

    moma(&g, &d, hmu, fault);
    if (*fault) {
	return 0;
    }
    s = hmu[0] * hmu[0];
    h2 = hmu[1] - s;
    *fault = h2 <= zero;
    if (*fault) {
	return 0;
    }
    t = sqrt(h2);
    h2a = t * h2;
    h2b = h2 * h2;
    h3 = hmu[2] - hmu[0] * (three * hmu[1] - two * s);
    rbet = h3 / h2a;
    h4 = hmu[3] - hmu[0] * (four * hmu[2] - hmu[0] * (six * hmu[1] - three * 
	    s));
    bet2 = h4 / h2b;
    w = g * d;
    u = d * d;

/*        GET DERIVATIVES */

    for (j = 1; j <= 2; ++j) {
	for (k = 1; k <= 4; ++k) {
	    t = (float) k;
	    if (j == 1) {
		goto L90;
	    }
	    s = ((w - t) * (hmu[k - 1] - hmu[k]) + (t + one) * (hmu[k] - hmu[
		    k + 1])) / u;
	    goto L100;
L90:
	    s = hmu[k] - hmu[k - 1];
L100:
	    dd[k - 1] = t * s / d;
/* L110: */
	}
	t = two * hmu[0] * dd[0];
	s = hmu[0] * dd[1];
	y = dd[1] - t;
	deriv[j - 1] = (dd[2] - three * (s + hmu[1] * dd[0] - t * hmu[0]) - 
		one5 * h3 * y / h2) / h2a;
	deriv[j + 1] = (dd[3] - four * (dd[2] * hmu[0] + dd[0] * hmu[2]) + 
		six * (hmu[1] * t + hmu[0] * (s - t * hmu[0])) - two * h4 * y 
		/ h2) / h2b;
/* L120: */
    }
    t = one / (deriv[0] * deriv[3] - deriv[1] * deriv[2]);
    u = (deriv[3] * (rbet - rb1) - deriv[1] * (bet2 - *b2)) * t;
    y = (deriv[0] * (bet2 - *b2) - deriv[2] * (rbet - rb1)) * t;

/*        FORM NEW ESTIMATES OF G AND D */

    g -= u;
    if (b1 == zero || g < zero) {
	g = zero;
    }
    d -= y;
    if (fabs(u) > tt || fabs(y) > tt) {
	goto L80;
    }

/*        END OF ITERATION */

    *delta = d;
    *xlam = *sigma / sqrt(h2);
    if (neg) {
	goto L130;
    }
    *gamma = g;
    goto L140;
L130:
    *gamma = -g;
    hmu[0] = one - hmu[0];
L140:
    *xi = *xbar - *xlam * hmu[0];
    return 0;
} /* sbfit_ */


/* Subroutine */ int moma(g, d, a, fault)
float *g, *d,*a;
int *fault;
{
    /* Initialized data */

    static float zz = 1e-5f;
    static float vv = 1e-8f;
    static int limit = 500;
    static float rttwo = 1.414213562f;
    static float rrtpi = .5641895835f;
    static float expa = 80.f;
    static float expb = 23.7f;
    static float zero = 0.f;
    static float quart = .25f;
    static float half = .5f;
    static float p75 = .75f;
    static float one = 1.f;
    static float two = 2.f;
    static float three = 3.f;

    /* System generated locals */
    float r__1;

    /* Builtin functions */

    /* Local variables */
    float b[6], c[6], e, f, h;
    int i, k;
    int l;
    int m;
    float p, q, r, s, t, u, v, w, x, y, z, aa, ab;


/*        ALGORITHM AS 99.3  APPL. STATIST. (1976) VOL.25, P.180 */

/*        EVALUATES FIRST SIX MOMENTS OF A JOHNSON */
/*        SB DISTRIBUTION, USING GOODWIN METHOD */


    /* Parameter adjustments */
    --a;

    /* Function Body */

/*        RTTWO IS SQRT(2.0) */
/*        RRTPI IS RECIPROCAL OF SQRT(PI) */
/*        EXPA IS A VALUE SUCH THAT EXP(EXPA) DOES NOT QUITE */
/*          CAUSE OVERFLOW */
/*        EXPB IS A VALUE SUCH THAT 1.0 + EXP(-EXPB) MAY BE */
/*          TAKEN TO BE 1.0 */



    *fault = 0;
    for (i = 1; i <= 6; ++i) {
/* L10: */
	c[i - 1] = zero;
    }
    w = *g / *d;

/*        TRIAL VALUE OF H */

    if (w > expa) {
	goto L140;
    }
    e = exp(w) + one;
    r = rttwo / *d;
    h = p75;
    if (*d < three) {
	h = quart * *d;
    }
    k = 1;
    goto L40;

/*        START OF OUTER LOOP */

L20:
    ++k;
    if (k > limit) {
	goto L140;
    }
    for (i = 1; i <= 6; ++i) {
/* L30: */
	c[i - 1] = a[i];
    }

/*        NO CONVERGENCE YET - TRY SMALLER H */

    h = half * h;
L40:
    t = w;
    u = t;
    y = h * h;
    x = two * y;
    a[1] = one / e;
    for (i = 2; i <= 6; ++i) {
/* L50: */
	a[i] = a[i - 1] / e;
    }
    v = y;
    f = r * h;
    m = 0;

/*        START OF INNER LOOP */
/*        TO EVALUATE INFINITE SERIES */

L60:
    ++m;
    if (m > limit) {
	goto L140;
    }
    for (i = 1; i <= 6; ++i) {
/* L70: */
	b[i - 1] = a[i];
    }
    u -= f;
    z = one;
    if (u > -expb) {
	z = exp(u) + z;
    }
    t += f;
    l = t > expb;
    if (! l) {
	s = exp(t) + one;
    }
    r__1 = -v;
    p = exp(r__1);
    q = p;
    for (i = 1; i <= 6; ++i) {
	aa = a[i];
	p /= z;
	ab = aa;
	aa += p;
	if (aa == ab) {
	    goto L100;
	}
	if (l) {
	    goto L80;
	}
	q /= s;
	ab = aa;
	aa += q;
	l = aa == ab;
L80:
	a[i] = aa;
/* L90: */
    }
L100:
    y += x;
    v += y;
    for (i = 1; i <= 6; ++i) {
	if (a[i] == zero) {
	    goto L140;
	}
	r__1 = (a[i] - b[i - 1]) / a[i];
	if (fabs(r__1) > vv) {
	    goto L60;
	}
/* L110: */
    }

/*        END OF INNER LOOP */

    v = rrtpi * h;
    for (i = 1; i <= 6; ++i) {
/* L120: */
	a[i] = v * a[i];
    }
    for (i = 1; i <= 6; ++i) {
	if (a[i] == zero) {
	    goto L140;
	}
	r__1 = (a[i] - c[i - 1]) / a[i];
	if (fabs(r__1) > zz) {
	    goto L20;
	}
/* L130: */
    }

/*        END OF OUTER LOOP */

    return 0;
L140:
    *fault = 1;
    return 0;
} /* mom */

