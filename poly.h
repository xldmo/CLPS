#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "config.h"
#include "params.h"

typedef struct {
  int32_t coeffs[N];
} poly;

#define poly_reduce DILITHIUM_NAMESPACE(poly_reduce)
void poly_reduce(poly *a);
#define poly_caddq DILITHIUM_NAMESPACE(poly_caddq)
void poly_caddq(poly *a);

#define poly_add DILITHIUM_NAMESPACE(poly_add)
void poly_add(poly *c, const poly *a, const poly *b);
#define poly_sub DILITHIUM_NAMESPACE(poly_sub)
void poly_sub(poly *c, const poly *a, const poly *b);
#define poly_shiftl DILITHIUM_NAMESPACE(poly_shiftl)
void poly_shiftl(poly *a);

#define poly_ntt DILITHIUM_NAMESPACE(poly_ntt)
void poly_ntt(poly *a);
#define poly_invntt_tomont DILITHIUM_NAMESPACE(poly_invntt_tomont)
void poly_invntt_tomont(poly *a);
#define poly_pointwise_montgomery DILITHIUM_NAMESPACE(poly_pointwise_montgomery)
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b);

#define poly_power2round DILITHIUM_NAMESPACE(poly_power2round)
void poly_power2round(poly *a1, poly *a0, const poly *a);
#define poly_decompose DILITHIUM_NAMESPACE(poly_decompose)
void poly_decompose(poly *a1, poly *a0, const poly *a);
#define sclls_poly_make_hw DILITHIUM_NAMESPACE(sclls_poly_make_hw)
void sclls_poly_make_hw(poly *hw, const poly *a1, const poly *b1, const poly *c1);
/* CLPS Makeh_w: hw = d1 - a1 - b1 - c1 (d1 = HighBits(a+b+c)) */
#define clps_poly_make_hw3 CLPS_NAMESPACE(clps_poly_make_hw3)
void clps_poly_make_hw3(poly *hw, const poly *a1, const poly *b1,
                        const poly *c1, const poly *d1);
#define poly_make_hint DILITHIUM_NAMESPACE(poly_make_hint)
unsigned int poly_make_hint(poly *h, const poly *a0, const poly *a1);
#define poly_use_hint DILITHIUM_NAMESPACE(poly_use_hint)
void poly_use_hint(poly *b, const poly *a, const poly *h);

#define poly_chknorm DILITHIUM_NAMESPACE(poly_chknorm)
int poly_chknorm(const poly *a, int32_t B);
#define poly_uniform DILITHIUM_NAMESPACE(poly_uniform)
void poly_uniform(poly *a,
                  const uint8_t seed[SEEDBYTES],
                  uint16_t nonce);
#define poly_uniform_eta DILITHIUM_NAMESPACE(poly_uniform_eta)
void poly_uniform_eta(poly *a,
                      const uint8_t seed[CRHBYTES],
                      uint16_t nonce);
#define poly_uniform_gamma1 DILITHIUM_NAMESPACE(poly_uniform_gamma1)
void poly_uniform_gamma1(poly *a,
                         const uint8_t seed[CRHBYTES],
                         uint16_t nonce);
#define poly_challenge DILITHIUM_NAMESPACE(poly_challenge)
void poly_challenge(poly *c, const uint8_t seed[CTILDEBYTES]);

#define polyeta_pack DILITHIUM_NAMESPACE(polyeta_pack)
void polyeta_pack(uint8_t *r, const poly *a);
#define polyeta_unpack DILITHIUM_NAMESPACE(polyeta_unpack)
void polyeta_unpack(poly *r, const uint8_t *a);

#define polyt1_pack DILITHIUM_NAMESPACE(polyt1_pack)
void polyt1_pack(uint8_t *r, const poly *a);
#define polyt1_unpack DILITHIUM_NAMESPACE(polyt1_unpack)
void polyt1_unpack(poly *r, const uint8_t *a);

#define sclls_polyt_pack DILITHIUM_NAMESPACE(sclls_polyt_pack)
void sclls_polyt_pack(uint8_t *r, const poly *a);
#define sclls_polyt_unpack DILITHIUM_NAMESPACE(sclls_polyt_unpack)
void sclls_polyt_unpack(poly *r, const uint8_t *a);

#define polyt0_pack DILITHIUM_NAMESPACE(polyt0_pack)
void polyt0_pack(uint8_t *r, const poly *a);
#define polyt0_unpack DILITHIUM_NAMESPACE(polyt0_unpack)
void polyt0_unpack(poly *r, const uint8_t *a);

#define sclls_polyz0_pack DILITHIUM_NAMESPACE(sclls_polyz0_pack)
void sclls_polyz0_pack(uint8_t *r, const poly *a);
#define sclls_polyz0_unpack DILITHIUM_NAMESPACE(sclls_polyz0_unpack)
void sclls_polyz0_unpack(poly *r, const uint8_t *a);

/* CLPS z_0_A: 22-bit per coefficient, offset CLPS_Z0A_BOUND */
#define clps_polyz0a_pack CLPS_NAMESPACE(clps_polyz0a_pack)
void clps_polyz0a_pack(uint8_t *r, const poly *a);
#define clps_polyz0a_unpack CLPS_NAMESPACE(clps_polyz0a_unpack)
void clps_polyz0a_unpack(poly *r, const uint8_t *a);

#define polyz_pack DILITHIUM_NAMESPACE(polyz_pack)
void polyz_pack(uint8_t *r, const poly *a);
#define polyz_unpack DILITHIUM_NAMESPACE(polyz_unpack)
void polyz_unpack(poly *r, const uint8_t *a);

#define sclls_polyhw_pack DILITHIUM_NAMESPACE(sclls_polyhw_pack)
void sclls_polyhw_pack(uint8_t *r, const poly *a);
#define sclls_polyhw_unpack DILITHIUM_NAMESPACE(sclls_polyhw_unpack)
void sclls_polyhw_unpack(poly *a, const uint8_t *r);

#define polyw1_pack DILITHIUM_NAMESPACE(polyw1_pack)
void polyw1_pack(uint8_t *r, const poly *a);
#define sclls_polyw1_unpack DILITHIUM_NAMESPACE(sclls_polyw1_unpack)
void sclls_polyw1_unpack(poly *r, const uint8_t *a);

#endif
