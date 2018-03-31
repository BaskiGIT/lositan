#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* START OF AS 100 */
/* Subroutine */ int varj(ajv, snv, itype, gamma, 
	delta, xlam, xi, ifault)
	float *ajv, *snv, *gamma, *delta, *xlam, *xi;
	int *itype, *ifault;
{
    /* Initialized data */

    static float zero = 0.f;
    static float half = .5f;
    static float one = 1.f;

    /* System generated locals */
    float r__1;

    /* Builtin functions */

    /* Local variables */
    float v, w;


/*        ALGORITHM AS 100.1  APPL. STATIST. (1976) VOL.25, P.190 */

/*        CONVERTS A STANDARD NORMAL VARIATE (SNV) TO A */
/*        JOHNSON VARIATE (AJV) */




    *ajv = zero;
    *ifault = 1;
    if (*itype < 1 || *itype > 4) {
	return 0;
    }
    *ifault = 0;
    switch (*itype) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
    }

/*        SL DISTRIBUTION */

L10:
    r__1 = (*xlam * *snv - *gamma) / *delta;
    *ajv = *xlam * exp(r__1) + *xi;
    return 0;

/*        SU DISTRIBUTION */

L20:
    r__1 = (*snv - *gamma) / *delta;
    w = exp(r__1);
    w = half * (w - one / w);
    *ajv = *xlam * w + *xi;
    return 0;

/*        SB DISTRIBUTION */

L30:
    w = (*snv - *gamma) / *delta;
    r__1 = -fabs(w);
    v = exp(r__1);
    v = (one - v) / (one + v);
    *ajv = half * *xlam * ((w < 0.0 ? -v : v) + one) + *xi;
    return 0;

/*        NORMAL DISTRIBUTION */

L40:
    *ajv = (*snv - *gamma) / *delta;
    return 0;
} /* varj_ */


/* Subroutine */ int varn(ajv, snv, itype, gamma, 
	delta, xlam, xi, ifault)
	float *ajv, *snv, *gamma, *delta, *xlam, *xi;
	int *itype, *ifault;
{
    /* Initialized data */

    static float zero = 0.f;
    static float half = .5f;
    static float one = 1.f;
    static float c = -63.f;

    /* System generated locals */
    float r__1;

    /* Builtin functions */

    /* Local variables */
    float v, w;


/*        ALGORITHM AS 100.2  APPL. STATIST. (1976) VOL.25, P.190 */

/*        CONVERTS A JOHNSON VARIATE (AJV) TO A */
/*        STANDARD NORMAL VARIATE (SNV) */




    *snv = zero;
    *ifault = 1;
    if (*itype < 1 || *itype > 4) {
	return 0;
    }
    *ifault = 0;
    switch (*itype) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
    }

/*        SL DISTRIBUTION */

L10:
    w = *xlam * (*ajv - *xi);
    if (w <= zero) {
	goto L15;
    }
    *snv = *xlam * (log(w) * *delta + *gamma);
    return 0;
L15:
    *ifault = 2;
    return 0;

/*        SU DISTRIBUTION */

L20:
    w = (*ajv - *xi) / *xlam;
    if (w > c) {
	goto L23;
    }
    w = -half / w;
    goto L27;
L23:
    r__1 = w * w + one;
    w = sqrt(r__1) + w;
L27:
    *snv = log(w) * *delta + *gamma;
    return 0;

/*        SB DISTRIBUTION */

L30:
    w = *ajv - *xi;
    v = *xlam - w;
    if (w <= zero || v <= zero) {
	goto L35;
    }
    r__1 = w / v;
    *snv = log(r__1) * *delta + *gamma;
    return 0;
L35:
    *ifault = 2;
    return 0;

/*        NORMAL DISTRIBUTION */

L40:
    *snv = *delta * *ajv + *gamma;
    return 0;
} /* varn_ */

