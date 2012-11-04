#include "randomext.h"


extern void randomext_standard_normal_init(VALUE cRandom);
extern void randomext_gamma_init(VALUE cRandom);
extern void randomext_binomial_init(VALUE cRandom);
extern void randomext_poisson_init(VALUE cRandom);
extern void randomext_hypergeometric_init(VALUE cRandom);

void Init_randomext_native()
{
  VALUE cRandom = rb_const_get(rb_cObject, rb_intern("Random"));
  
  randomext_standard_normal_init(cRandom);
  randomext_gamma_init(cRandom);
  randomext_binomial_init(cRandom);
  randomext_poisson_init(cRandom);
  randomext_hypergeometric_init(cRandom);
}
