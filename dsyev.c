#include "blaswrap.h"
#include "f2c.h"
int dsyev_(char *jobz, char *uplo, integer *n, doublereal *a,
integer *lda, doublereal *w, doublereal *work, integer *lwork, 
integer *info)
{
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__0 = 0;
    static doublereal c_b17 = 1.;
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;
    double sqrt(doublereal);
    static integer inde;
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax;
    static integer lopt;
    extern int dscal_(integer *, doublereal *, doublereal *, 
    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *);
    static integer iinfo;
    static logical lower, wantz;
    static integer nb;
    extern doublereal dlamch_(char *);
    static integer iscale;
    extern int dlascl_(char *, integer *, integer *, 
    doublereal *, doublereal *, integer *, integer *, doublereal *, 
    integer *, integer *);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
    integer *, integer *, ftnlen, ftnlen);
    extern int xerbla_(char *, integer *);
    static doublereal bignum;
    static integer indtau;
    extern int dsterf_(integer *, doublereal *, doublereal *,
     integer *);
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
    integer *, doublereal *);
    static integer indwrk;
    extern int dorgtr_(char *, integer *, doublereal *, 
    integer *, doublereal *, doublereal *, integer *, integer *), dsteqr_(char *, integer *, doublereal *, doublereal *, 
    doublereal *, integer *, doublereal *, integer *), 
    dsytrd_(char *, integer *, doublereal *, integer *, doublereal *, 
    doublereal *, doublereal *, doublereal *, integer *, integer *);
    static integer llwork;
    static doublereal smlnum;
    static integer lwkopt;
    static logical lquery;
    static doublereal eps;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --w;
    --work;
    wantz = lsame_(jobz, "V");
    lower = lsame_(uplo, "L");
    lquery = *lwork == -1;
 
    *info = 0;
    if (! (wantz || lsame_(jobz, "N"))) {
 
*info = -1;
     
} else if (! (lower || lsame_(uplo, "U"))) {
 
*info = -2;
     
} else if (*n < 0) {
 
*info = -3;
     
} else if (*lda < max(1,*n)) {
 
*info = -5;
     
} 
else  {
i__1 = 1, i__2 = *n * 3 - 1;
if (*lwork < max(i__1,i__2) && ! lquery) {
 
    *info = -8;
 
}
     
}
 
    if (*info == 0) {
 
nb = ilaenv_(&c__1, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
 (ftnlen)1);
i__1 = 1, i__2 = (nb + 2) * *n;
lwkopt = max(i__1,i__2);
work[1] = (doublereal) lwkopt;
     
}
 
    if (*info != 0) {
 
i__1 = -(*info);
xerbla_("DSYEV ", &i__1);
return 0;
     
} else if (lquery) {
 
return 0;
     
}
    if (*n == 0) {
 
    work[1] = 1.;
    return 0;
     
}
 
    if (*n == 1) {
 
w[1] = a_ref(1, 1);
work[1] = 3.;
if (wantz) {
 
    a_ref(1, 1) = 1.;
 
}
return 0;
     
}
    safmin = dlamch_("Safe minimum");
    eps = dlamch_("Precision");
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);
    anrm = dlansy_("M", uplo, n, &a[a_offset], lda, &work[1]);
    iscale = 0;
    if (anrm > 0. && anrm < rmin) {
 
iscale = 1;
sigma = rmin / anrm;
     
} else if (anrm > rmax) {
 
iscale = 1;
sigma = rmax / anrm;
     
}
    if (iscale == 1) {
 
dlascl_(uplo, &c__0, &c__0, &c_b17, &sigma, n, n, &a[a_offset], lda, 
info);
     
}
    inde = 1;
    indtau = inde + *n;
    indwrk = indtau + *n;
    llwork = *lwork - indwrk + 1;
    dsytrd_(uplo, n, &a[a_offset], lda, &w[1], &work[inde], &work[indtau], &
    work[indwrk], &llwork, &iinfo);
    lopt = (integer) ((*n << 1) + work[indwrk]);
    if (! wantz) {
 
dsterf_(n, &w[1], &work[inde], info);
     
} else {
 
dorgtr_(uplo, n, &a[a_offset], lda, &work[indtau], &work[indwrk], &
llwork, &iinfo);
dsteqr_(jobz, n, &w[1], &work[inde], &a[a_offset], lda, &work[indtau],
 info);
     
}
    if (iscale == 1) {
 
if (*info == 0) {
 
    imax = *n;
 
} else {
 
    imax = *info - 1;
 
}
d__1 = 1. / sigma;
dscal_(&imax, &d__1, &w[1], &c__1);
     
}
    work[1] = (doublereal) lwkopt;
 
    return 0;
}
