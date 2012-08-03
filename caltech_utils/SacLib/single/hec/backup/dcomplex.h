/* double-precision complex arithmetic subroutines */

typedef struct { double re; double im;}  dcomplex;

extern dcomplex dcmult(),dconj(),dcdiv(),dcsum(),dcdiff();
extern dcomplex dscamult(),dscadiv(),zsqrt();
extern double dmodu();
