#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "rounding.h"
#include "symmetric.h"

#ifdef DBENCH
#include "test/cpucycles.h"
extern const uint64_t timing_overhead;
extern uint64_t *tred, *tadd, *tmul, *tround, *tsample, *tpack;
#define DBENCH_START() uint64_t time = cpucycles()
#define DBENCH_STOP(t) t += cpucycles() - time - timing_overhead
#else
#define DBENCH_START()
#define DBENCH_STOP(t)
#endif

/*************************************************
* Name:        poly_reduce
*
* Description: Inplace reduction of all coefficients of polynomial to
*              representative in [-6283008,6283008].
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_reduce(poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    a->coeffs[i] = reduce32(a->coeffs[i]);

  DBENCH_STOP(*tred);
}

/*************************************************
* Name:        poly_caddq
*
* Description: For all coefficients of in/out polynomial add Q if
*              coefficient is negative.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_caddq(poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    a->coeffs[i] = caddq(a->coeffs[i]);

  DBENCH_STOP(*tred);
}

/*************************************************
* Name:        poly_add
*
* Description: Add polynomials. No modular reduction is performed.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first summand
*              - const poly *b: pointer to second summand
**************************************************/
void poly_add(poly *c, const poly *a, const poly *b)  {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    c->coeffs[i] = a->coeffs[i] + b->coeffs[i];

  DBENCH_STOP(*tadd);
}

/*************************************************
* Name:        poly_sub
*
* Description: Subtract polynomials. No modular reduction is
*              performed.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial to be
*                               subtraced from first input polynomial
**************************************************/
void poly_sub(poly *c, const poly *a, const poly *b) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    c->coeffs[i] = a->coeffs[i] - b->coeffs[i];

  DBENCH_STOP(*tadd);
}

/*************************************************
* Name:        poly_shiftl
*
* Description: Multiply polynomial by 2^D without modular reduction. Assumes
*              input coefficients to be less than 2^{31-D} in absolute value.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_shiftl(poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    a->coeffs[i] <<= D;

  DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        poly_ntt
*
* Description: Inplace forward NTT. Coefficients can grow by
*              8*Q in absolute value.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_ntt(poly *a) {
  DBENCH_START();

  ntt(a->coeffs);

  DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        poly_invntt_tomont
*
* Description: Inplace inverse NTT and multiplication by 2^{32}.
*              Input coefficients need to be less than Q in absolute
*              value and output coefficients are again bounded by Q.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_invntt_tomont(poly *a) {
  DBENCH_START();

  invntt_tomont(a->coeffs);

  DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        poly_pointwise_montgomery
*
* Description: Pointwise multiplication of polynomials in NTT domain
*              representation and multiplication of resulting polynomial
*              by 2^{-32}.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    c->coeffs[i] = montgomery_reduce((int64_t)a->coeffs[i] * b->coeffs[i]);

  DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        poly_power2round
*
* Description: For all coefficients c of the input polynomial,
*              compute c0, c1 such that c mod Q = c1*2^D + c0
*              with -2^{D-1} < c0 <= 2^{D-1}. Assumes coefficients to be
*              standard representatives.
*
* Arguments:   - poly *a1: pointer to output polynomial with coefficients c1
*              - poly *a0: pointer to output polynomial with coefficients c0
*              - const poly *a: pointer to input polynomial
**************************************************/
void poly_power2round(poly *a1, poly *a0, const poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    a1->coeffs[i] = power2round(&a0->coeffs[i], a->coeffs[i]);

  DBENCH_STOP(*tround);
}

/*************************************************
* Name:        poly_decompose
*
* Description: For all coefficients c of the input polynomial,
*              compute high and low bits c0, c1 such c mod Q = c1*ALPHA + c0
*              with -ALPHA/2 < c0 <= ALPHA/2 except c1 = (Q-1)/ALPHA where we
*              set c1 = 0 and -ALPHA/2 <= c0 = c mod Q - Q < 0.
*              Assumes coefficients to be standard representatives.
*
* Arguments:   - poly *a1: pointer to output polynomial with coefficients c1
*              - poly *a0: pointer to output polynomial with coefficients c0
*              - const poly *a: pointer to input polynomial
**************************************************/
void poly_decompose(poly *a1, poly *a0, const poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    a1->coeffs[i] = decompose(&a0->coeffs[i], a->coeffs[i]);

  DBENCH_STOP(*tround);
}

/*************************************************
* Name:        sclls_poly_make_hw
*
* Description: For each coefficient i, set hw[i] = sclls_make_hw(a1[i],b1[i],c1[i])
*              i.e. (c1 - a1 - b1) as unsigned int per coefficient.
*
* Arguments:   - poly *hw: pointer to output polynomial
*              - const poly *a1: first decomposed high part (w1-style)
*              - const poly *b1: second decomposed high part
*              - const poly *c1: third decomposed high part
**************************************************/
void sclls_poly_make_hw(poly *hw, const poly *a1, const poly *b1, const poly *c1) {
  unsigned int i;
  DBENCH_START();
  for(i = 0; i < N; ++i) {
    hw->coeffs[i] = sclls_make_hw(a1->coeffs[i], b1->coeffs[i], c1->coeffs[i]);
  }
  DBENCH_STOP(*tround);
}

/*************************************************
* Name:        clps_poly_make_hw3
*
* Description: CLPS Makeh_w on one polynomial:
*              hw[i] = d1[i] - a1[i] - b1[i] - c1[i]
*              where a1,b1,c1 = HighBits of the three summands and
*              d1 = HighBits(a+b+c) (all decomposed with alpha = 2*GAMMA2).
*
* Arguments:   - poly *hw: pointer to output polynomial
*              - const poly *a1: first decomposed high part
*              - const poly *b1: second decomposed high part
*              - const poly *c1: third decomposed high part
*              - const poly *d1: high part of the sum
**************************************************/
void clps_poly_make_hw3(poly *hw, const poly *a1, const poly *b1,
                        const poly *c1, const poly *d1) {
  unsigned int i;
  DBENCH_START();
  for(i = 0; i < N; ++i) {
    hw->coeffs[i] = d1->coeffs[i] - a1->coeffs[i] - b1->coeffs[i] - c1->coeffs[i];
  }
  DBENCH_STOP(*tround);
}

/*************************************************
* Name:        poly_make_hint
*
* Description: Compute hint polynomial. The coefficients of which indicate
*              whether the low bits of the corresponding coefficient of
*              the input polynomial overflow into the high bits.
*
* Arguments:   - poly *h: pointer to output hint polynomial
*              - const poly *a0: pointer to low part of input polynomial
*              - const poly *a1: pointer to high part of input polynomial
*
* Returns number of 1 bits.
**************************************************/
unsigned int poly_make_hint(poly *h, const poly *a0, const poly *a1) {
  unsigned int i, s = 0;
  DBENCH_START();

  for(i = 0; i < N; ++i) {
    h->coeffs[i] = make_hint(a0->coeffs[i], a1->coeffs[i]);
    s += h->coeffs[i];
  }

  DBENCH_STOP(*tround);
  return s;
}

/*************************************************
* Name:        poly_use_hint
*
* Description: Use hint polynomial to correct the high bits of a polynomial.
*
* Arguments:   - poly *b: pointer to output polynomial with corrected high bits
*              - const poly *a: pointer to input polynomial
*              - const poly *h: pointer to input hint polynomial
**************************************************/
void poly_use_hint(poly *b, const poly *a, const poly *h) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    b->coeffs[i] = use_hint(a->coeffs[i], h->coeffs[i]);

  DBENCH_STOP(*tround);
}

/*************************************************
* Name:        poly_chknorm
*
* Description: Check infinity norm of polynomial against given bound.
*              Assumes input coefficients were reduced by reduce32().
*
* Arguments:   - const poly *a: pointer to polynomial
*              - int32_t B: norm bound
*
* Returns 0 if norm is strictly smaller than B <= (Q-1)/8 and 1 otherwise.
**************************************************/
int poly_chknorm(const poly *a, int32_t B) {
  unsigned int i;
  int32_t t;
  DBENCH_START();

  /* B must fit comfortably so that the absolute-value trick below
   * (t = a - (t & 2a)) does not overflow int32. The trick only needs
   * |a| <= INT32_MAX (so 2a doesn't overflow), but conservatively we
   * require B <= (Q-1)/2. Larger B values (e.g. MODE 3 z0 bound
   * 2*(GAMMA1 - BETA) = 1048184, which is between (Q-1)/8=1047552
   * and (Q-1)/2=4190208) are still safe to check, so only reject if
   * B exceeds the half-modulus. */
  if(B > (Q-1)/2)
    return 1;

  /* It is ok to leak which coefficient violates the bound since
     the probability for each coefficient is independent of secret
     data but we must not leak the sign of the centralized representative. */
  for(i = 0; i < N; ++i) {
    /* Absolute value */
    t = a->coeffs[i] >> 31;
    t = a->coeffs[i] - (t & 2*a->coeffs[i]);

    if(t >= B) {
      DBENCH_STOP(*tsample);
      return 1;
    }
  }

  DBENCH_STOP(*tsample);
  return 0;
}

/*************************************************
* Name:        rej_uniform
*
* Description: Sample uniformly random coefficients in [0, Q-1] by
*              performing rejection sampling on array of random bytes.
*
* Arguments:   - int32_t *a: pointer to output array (allocated)
*              - unsigned int len: number of coefficients to be sampled
*              - const uint8_t *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/
static unsigned int rej_uniform(int32_t *a,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos;
  uint32_t t;
  DBENCH_START();

  ctr = pos = 0;
  while(ctr < len && pos + 3 <= buflen) {
    t  = buf[pos++];
    t |= (uint32_t)buf[pos++] << 8;
    t |= (uint32_t)buf[pos++] << 16;
    t &= 0x7FFFFF;

    if(t < Q)
      a[ctr++] = t;
  }

  DBENCH_STOP(*tsample);
  return ctr;
}

/*************************************************
* Name:        poly_uniform
*
* Description: Sample polynomial with uniformly random coefficients
*              in [0,Q-1] by performing rejection sampling on the
*              output stream of SHAKE128(seed|nonce)
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length SEEDBYTES
*              - uint16_t nonce: 2-byte nonce
**************************************************/
#define POLY_UNIFORM_NBLOCKS ((768 + STREAM128_BLOCKBYTES - 1)/STREAM128_BLOCKBYTES)
void poly_uniform(poly *a,
                  const uint8_t seed[SEEDBYTES],
                  uint16_t nonce)
{
  unsigned int i, ctr, off;
  unsigned int buflen = POLY_UNIFORM_NBLOCKS*STREAM128_BLOCKBYTES;
  uint8_t buf[POLY_UNIFORM_NBLOCKS*STREAM128_BLOCKBYTES + 2];
  stream128_state state;

  stream128_init(&state, seed, nonce);
  stream128_squeezeblocks(buf, POLY_UNIFORM_NBLOCKS, &state);

  ctr = rej_uniform(a->coeffs, N, buf, buflen);

  while(ctr < N) {
    off = buflen % 3;
    for(i = 0; i < off; ++i)
      buf[i] = buf[buflen - off + i];

    stream128_squeezeblocks(buf + off, 1, &state);
    buflen = STREAM128_BLOCKBYTES + off;
    ctr += rej_uniform(a->coeffs + ctr, N - ctr, buf, buflen);
  }
}

/*************************************************
* Name:        rej_eta
*
* Description: Sample uniformly random coefficients in [-ETA, ETA] by
*              performing rejection sampling on array of random bytes.
*
* Arguments:   - int32_t *a: pointer to output array (allocated)
*              - unsigned int len: number of coefficients to be sampled
*              - const uint8_t *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/
static unsigned int rej_eta(int32_t *a,
                            unsigned int len,
                            const uint8_t *buf,
                            unsigned int buflen)
{
  unsigned int ctr, pos;
  uint32_t t0, t1;
  DBENCH_START();

  ctr = pos = 0;
  while(ctr < len && pos < buflen) {
    t0 = buf[pos] & 0x0F;
    t1 = buf[pos++] >> 4;

#if ETA == 2
    if(t0 < 15) {
      t0 = t0 - (205*t0 >> 10)*5;
      a[ctr++] = 2 - t0;
    }
    if(t1 < 15 && ctr < len) {
      t1 = t1 - (205*t1 >> 10)*5;
      a[ctr++] = 2 - t1;
    }
#elif ETA == 4
    if(t0 < 9)
      a[ctr++] = 4 - t0;
    if(t1 < 9 && ctr < len)
      a[ctr++] = 4 - t1;
#endif
  }

  DBENCH_STOP(*tsample);
  return ctr;
}

/*************************************************
* Name:        poly_uniform_eta
*
* Description: Sample polynomial with uniformly random coefficients
*              in [-ETA,ETA] by performing rejection sampling on the
*              output stream from SHAKE256(seed|nonce)
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length CRHBYTES
*              - uint16_t nonce: 2-byte nonce
**************************************************/
#if ETA == 2
#define POLY_UNIFORM_ETA_NBLOCKS ((136 + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
#elif ETA == 4
#define POLY_UNIFORM_ETA_NBLOCKS ((227 + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
#endif
void poly_uniform_eta(poly *a,
                      const uint8_t seed[CRHBYTES],
                      uint16_t nonce)
{
  unsigned int ctr;
  unsigned int buflen = POLY_UNIFORM_ETA_NBLOCKS*STREAM256_BLOCKBYTES;
  uint8_t buf[POLY_UNIFORM_ETA_NBLOCKS*STREAM256_BLOCKBYTES];
  stream256_state state;

  stream256_init(&state, seed, nonce);
  stream256_squeezeblocks(buf, POLY_UNIFORM_ETA_NBLOCKS, &state);

  ctr = rej_eta(a->coeffs, N, buf, buflen);

  while(ctr < N) {
    stream256_squeezeblocks(buf, 1, &state);
    ctr += rej_eta(a->coeffs + ctr, N - ctr, buf, STREAM256_BLOCKBYTES);
  }
}

/*************************************************
* Name:        poly_uniform_gamma1m1
*
* Description: Sample polynomial with uniformly random coefficients
*              in [-(GAMMA1 - 1), GAMMA1] by unpacking output stream
*              of SHAKE256(seed|nonce)
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length CRHBYTES
*              - uint16_t nonce: 16-bit nonce
**************************************************/
#define POLY_UNIFORM_GAMMA1_NBLOCKS ((POLYZ_PACKEDBYTES + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
void poly_uniform_gamma1(poly *a,
                         const uint8_t seed[CRHBYTES],
                         uint16_t nonce)
{
  uint8_t buf[POLY_UNIFORM_GAMMA1_NBLOCKS*STREAM256_BLOCKBYTES];
  stream256_state state;

  stream256_init(&state, seed, nonce);
  stream256_squeezeblocks(buf, POLY_UNIFORM_GAMMA1_NBLOCKS, &state);
  polyz_unpack(a, buf);
}

/*************************************************
* Name:        challenge
*
* Description: Implementation of H. Samples polynomial with TAU nonzero
*              coefficients in {-1,1} using the output stream of
*              SHAKE256(seed).
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const uint8_t mu[]: byte array containing seed of length CTILDEBYTES
**************************************************/
void poly_challenge(poly *c, const uint8_t seed[CTILDEBYTES]) {
  unsigned int i, b, pos;
  uint64_t signs;
  uint8_t buf[SHAKE256_RATE];
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, seed, CTILDEBYTES);
  shake256_finalize(&state);
  shake256_squeezeblocks(buf, 1, &state);

  signs = 0;
  for(i = 0; i < 8; ++i)
    signs |= (uint64_t)buf[i] << 8*i;
  pos = 8;

  for(i = 0; i < N; ++i)
    c->coeffs[i] = 0;
  for(i = N-TAU; i < N; ++i) {
    do {
      if(pos >= SHAKE256_RATE) {
        shake256_squeezeblocks(buf, 1, &state);
        pos = 0;
      }

      b = buf[pos++];
    } while(b > i);

    c->coeffs[i] = c->coeffs[b];
    c->coeffs[b] = 1 - 2*(signs & 1);
    signs >>= 1;
  }
}

/*************************************************
* Name:        polyeta_pack
*
* Description: Bit-pack polynomial with coefficients in [-ETA,ETA].
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYETA_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyeta_pack(uint8_t *r, const poly *a) {
  unsigned int i;
  uint8_t t[8];
  DBENCH_START();

#if ETA == 2
  for(i = 0; i < N/8; ++i) {
    t[0] = ETA - a->coeffs[8*i+0];
    t[1] = ETA - a->coeffs[8*i+1];
    t[2] = ETA - a->coeffs[8*i+2];
    t[3] = ETA - a->coeffs[8*i+3];
    t[4] = ETA - a->coeffs[8*i+4];
    t[5] = ETA - a->coeffs[8*i+5];
    t[6] = ETA - a->coeffs[8*i+6];
    t[7] = ETA - a->coeffs[8*i+7];

    r[3*i+0]  = (t[0] >> 0) | (t[1] << 3) | (t[2] << 6);
    r[3*i+1]  = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7);
    r[3*i+2]  = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5);
  }
#elif ETA == 4
  for(i = 0; i < N/2; ++i) {
    t[0] = ETA - a->coeffs[2*i+0];
    t[1] = ETA - a->coeffs[2*i+1];
    r[i] = t[0] | (t[1] << 4);
  }
#endif

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polyeta_unpack
*
* Description: Unpack polynomial with coefficients in [-ETA,ETA].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyeta_unpack(poly *r, const uint8_t *a) {
  unsigned int i;
  DBENCH_START();

#if ETA == 2
  for(i = 0; i < N/8; ++i) {
    r->coeffs[8*i+0] =  (a[3*i+0] >> 0) & 7;
    r->coeffs[8*i+1] =  (a[3*i+0] >> 3) & 7;
    r->coeffs[8*i+2] = ((a[3*i+0] >> 6) | (a[3*i+1] << 2)) & 7;
    r->coeffs[8*i+3] =  (a[3*i+1] >> 1) & 7;
    r->coeffs[8*i+4] =  (a[3*i+1] >> 4) & 7;
    r->coeffs[8*i+5] = ((a[3*i+1] >> 7) | (a[3*i+2] << 1)) & 7;
    r->coeffs[8*i+6] =  (a[3*i+2] >> 2) & 7;
    r->coeffs[8*i+7] =  (a[3*i+2] >> 5) & 7;

    r->coeffs[8*i+0] = ETA - r->coeffs[8*i+0];
    r->coeffs[8*i+1] = ETA - r->coeffs[8*i+1];
    r->coeffs[8*i+2] = ETA - r->coeffs[8*i+2];
    r->coeffs[8*i+3] = ETA - r->coeffs[8*i+3];
    r->coeffs[8*i+4] = ETA - r->coeffs[8*i+4];
    r->coeffs[8*i+5] = ETA - r->coeffs[8*i+5];
    r->coeffs[8*i+6] = ETA - r->coeffs[8*i+6];
    r->coeffs[8*i+7] = ETA - r->coeffs[8*i+7];
  }
#elif ETA == 4
  for(i = 0; i < N/2; ++i) {
    r->coeffs[2*i+0] = a[i] & 0x0F;
    r->coeffs[2*i+1] = a[i] >> 4;
    r->coeffs[2*i+0] = ETA - r->coeffs[2*i+0];
    r->coeffs[2*i+1] = ETA - r->coeffs[2*i+1];
  }
#endif

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polyt1_pack
*
* Description: Bit-pack polynomial t1 with coefficients fitting in 10 bits.
*              Input coefficients are assumed to be standard representatives.
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYT1_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyt1_pack(uint8_t *r, const poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N/4; ++i) {
    r[5*i+0] = (a->coeffs[4*i+0] >> 0);
    r[5*i+1] = (a->coeffs[4*i+0] >> 8) | (a->coeffs[4*i+1] << 2);
    r[5*i+2] = (a->coeffs[4*i+1] >> 6) | (a->coeffs[4*i+2] << 4);
    r[5*i+3] = (a->coeffs[4*i+2] >> 4) | (a->coeffs[4*i+3] << 6);
    r[5*i+4] = (a->coeffs[4*i+3] >> 2);
  }

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polyt1_unpack
*
* Description: Unpack polynomial t1 with 10-bit coefficients.
*              Output coefficients are standard representatives.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyt1_unpack(poly *r, const uint8_t *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N/4; ++i) {
    r->coeffs[4*i+0] = ((a[5*i+0] >> 0) | ((uint32_t)a[5*i+1] << 8)) & 0x3FF;
    r->coeffs[4*i+1] = ((a[5*i+1] >> 2) | ((uint32_t)a[5*i+2] << 6)) & 0x3FF;
    r->coeffs[4*i+2] = ((a[5*i+2] >> 4) | ((uint32_t)a[5*i+3] << 4)) & 0x3FF;
    r->coeffs[4*i+3] = ((a[5*i+3] >> 6) | ((uint32_t)a[5*i+4] << 2)) & 0x3FF;
  }

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polyt_pack
*
* Description: Bit-pack polynomial t with coefficients fitting in 23 bits.
*              Input coefficients are assumed to be standard representatives.
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYT_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void sclls_polyt_pack(uint8_t *r, const poly *a) {
  unsigned int i;
  uint64_t buffer;
  unsigned int bits;
  unsigned int out_idx;
  DBENCH_START();

  /* Tight LSB-first packing: 256 * 23 bits = 5888 bits = 736 bytes (no padding). */
  buffer = 0;
  bits = 0;
  out_idx = 0;
  for (i = 0; i < N; ++i) {
    uint32_t c = (uint32_t)(a->coeffs[i] & 0x7FFFFF);
    buffer |= ((uint64_t)c) << bits;
    bits += 23;
    while (bits >= 8) {
      r[out_idx++] = (uint8_t)(buffer & 0xFF);
      buffer >>= 8;
      bits -= 8;
    }
  }
  /* N*23 divisible by 8 => bits==0, out_idx==SCLLS_POLYT_PACKEDBYTES */

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polyt_unpack
*
* Description: Unpack polynomial t with 23-bit coefficients.
*              Output coefficients are standard representatives.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void sclls_polyt_unpack(poly *r, const uint8_t *a) {
  unsigned int i;
  uint64_t buffer;
  unsigned int bits;
  unsigned int in_idx;
  DBENCH_START();

  buffer = 0;
  bits = 0;
  in_idx = 0;
  for (i = 0; i < N; ++i) {
    while (bits < 23) {
      buffer |= ((uint64_t)a[in_idx++]) << bits;
      bits += 8;
    }
    r->coeffs[i] = (int32_t)(buffer & 0x7FFFFF);
    buffer >>= 23;
    bits -= 23;
  }

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polyt0_pack
*
* Description: Bit-pack polynomial t0 with coefficients in ]-2^{D-1}, 2^{D-1}].
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYT0_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyt0_pack(uint8_t *r, const poly *a) {
  unsigned int i;
  uint32_t t[8];
  DBENCH_START();

  for(i = 0; i < N/8; ++i) {
    t[0] = (1 << (D-1)) - a->coeffs[8*i+0];
    t[1] = (1 << (D-1)) - a->coeffs[8*i+1];
    t[2] = (1 << (D-1)) - a->coeffs[8*i+2];
    t[3] = (1 << (D-1)) - a->coeffs[8*i+3];
    t[4] = (1 << (D-1)) - a->coeffs[8*i+4];
    t[5] = (1 << (D-1)) - a->coeffs[8*i+5];
    t[6] = (1 << (D-1)) - a->coeffs[8*i+6];
    t[7] = (1 << (D-1)) - a->coeffs[8*i+7];

    r[13*i+ 0]  =  t[0];
    r[13*i+ 1]  =  t[0] >>  8;
    r[13*i+ 1] |=  t[1] <<  5;
    r[13*i+ 2]  =  t[1] >>  3;
    r[13*i+ 3]  =  t[1] >> 11;
    r[13*i+ 3] |=  t[2] <<  2;
    r[13*i+ 4]  =  t[2] >>  6;
    r[13*i+ 4] |=  t[3] <<  7;
    r[13*i+ 5]  =  t[3] >>  1;
    r[13*i+ 6]  =  t[3] >>  9;
    r[13*i+ 6] |=  t[4] <<  4;
    r[13*i+ 7]  =  t[4] >>  4;
    r[13*i+ 8]  =  t[4] >> 12;
    r[13*i+ 8] |=  t[5] <<  1;
    r[13*i+ 9]  =  t[5] >>  7;
    r[13*i+ 9] |=  t[6] <<  6;
    r[13*i+10]  =  t[6] >>  2;
    r[13*i+11]  =  t[6] >> 10;
    r[13*i+11] |=  t[7] <<  3;
    r[13*i+12]  =  t[7] >>  5;
  }

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polyt0_unpack
*
* Description: Unpack polynomial t0 with coefficients in ]-2^{D-1}, 2^{D-1}].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyt0_unpack(poly *r, const uint8_t *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N/8; ++i) {
    r->coeffs[8*i+0]  = a[13*i+0];
    r->coeffs[8*i+0] |= (uint32_t)a[13*i+1] << 8;
    r->coeffs[8*i+0] &= 0x1FFF;

    r->coeffs[8*i+1]  = a[13*i+1] >> 5;
    r->coeffs[8*i+1] |= (uint32_t)a[13*i+2] << 3;
    r->coeffs[8*i+1] |= (uint32_t)a[13*i+3] << 11;
    r->coeffs[8*i+1] &= 0x1FFF;

    r->coeffs[8*i+2]  = a[13*i+3] >> 2;
    r->coeffs[8*i+2] |= (uint32_t)a[13*i+4] << 6;
    r->coeffs[8*i+2] &= 0x1FFF;

    r->coeffs[8*i+3]  = a[13*i+4] >> 7;
    r->coeffs[8*i+3] |= (uint32_t)a[13*i+5] << 1;
    r->coeffs[8*i+3] |= (uint32_t)a[13*i+6] << 9;
    r->coeffs[8*i+3] &= 0x1FFF;

    r->coeffs[8*i+4]  = a[13*i+6] >> 4;
    r->coeffs[8*i+4] |= (uint32_t)a[13*i+7] << 4;
    r->coeffs[8*i+4] |= (uint32_t)a[13*i+8] << 12;
    r->coeffs[8*i+4] &= 0x1FFF;

    r->coeffs[8*i+5]  = a[13*i+8] >> 1;
    r->coeffs[8*i+5] |= (uint32_t)a[13*i+9] << 7;
    r->coeffs[8*i+5] &= 0x1FFF;

    r->coeffs[8*i+6]  = a[13*i+9] >> 6;
    r->coeffs[8*i+6] |= (uint32_t)a[13*i+10] << 2;
    r->coeffs[8*i+6] |= (uint32_t)a[13*i+11] << 10;
    r->coeffs[8*i+6] &= 0x1FFF;

    r->coeffs[8*i+7]  = a[13*i+11] >> 3;
    r->coeffs[8*i+7] |= (uint32_t)a[13*i+12] << 5;
    r->coeffs[8*i+7] &= 0x1FFF;

    r->coeffs[8*i+0] = (1 << (D-1)) - r->coeffs[8*i+0];
    r->coeffs[8*i+1] = (1 << (D-1)) - r->coeffs[8*i+1];
    r->coeffs[8*i+2] = (1 << (D-1)) - r->coeffs[8*i+2];
    r->coeffs[8*i+3] = (1 << (D-1)) - r->coeffs[8*i+3];
    r->coeffs[8*i+4] = (1 << (D-1)) - r->coeffs[8*i+4];
    r->coeffs[8*i+5] = (1 << (D-1)) - r->coeffs[8*i+5];
    r->coeffs[8*i+6] = (1 << (D-1)) - r->coeffs[8*i+6];
    r->coeffs[8*i+7] = (1 << (D-1)) - r->coeffs[8*i+7];
  }

  DBENCH_STOP(*tpack);
}
/*************************************************
* Name:        sclls_polyz0_pack
*
* Description: Bit-pack polynomial z0 with coefficients in
*              [-(GAMMA1 - 1), GAMMA1].
*
*              Encoding:
*              - If GAMMA1 == (1 << 17): 8 coefficients are packed into 19 bytes
*                (each coefficient uses 19 bits, little-endian bit stream).
*              - If GAMMA1 == (1 << 19): 8 coefficients are packed into 21 bytes
*                (each coefficient uses 21 bits, little-endian bit stream).
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYZ_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyz_pack(uint8_t *r, const poly *a) {
  unsigned int i;
  uint32_t t[4];
  DBENCH_START();

#if GAMMA1 == (1 << 17)
  for(i = 0; i < N/4; ++i) {
    t[0] = GAMMA1 - a->coeffs[4*i+0];
    t[1] = GAMMA1 - a->coeffs[4*i+1];
    t[2] = GAMMA1 - a->coeffs[4*i+2];
    t[3] = GAMMA1 - a->coeffs[4*i+3];

    r[9*i+0]  = t[0];
    r[9*i+1]  = t[0] >> 8;
    r[9*i+2]  = t[0] >> 16;
    r[9*i+2] |= t[1] << 2;
    r[9*i+3]  = t[1] >> 6;
    r[9*i+4]  = t[1] >> 14;
    r[9*i+4] |= t[2] << 4;
    r[9*i+5]  = t[2] >> 4;
    r[9*i+6]  = t[2] >> 12;
    r[9*i+6] |= t[3] << 6;
    r[9*i+7]  = t[3] >> 2;
    r[9*i+8]  = t[3] >> 10;
  }
#elif GAMMA1 == (1 << 19)
  for(i = 0; i < N/2; ++i) {
    t[0] = GAMMA1 - a->coeffs[2*i+0];
    t[1] = GAMMA1 - a->coeffs[2*i+1];

    r[5*i+0]  = t[0];
    r[5*i+1]  = t[0] >> 8;
    r[5*i+2]  = t[0] >> 16;
    r[5*i+2] |= t[1] << 4;
    r[5*i+3]  = t[1] >> 4;
    r[5*i+4]  = t[1] >> 12;
  }
#endif

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polyz_unpack
*
* Description: Unpack polynomial z with coefficients
*              in [-(GAMMA1 - 1), GAMMA1].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyz_unpack(poly *r, const uint8_t *a) {
  unsigned int i;
  DBENCH_START();

#if GAMMA1 == (1 << 17)
  for(i = 0; i < N/4; ++i) {
    r->coeffs[4*i+0]  = a[9*i+0];
    r->coeffs[4*i+0] |= (uint32_t)a[9*i+1] << 8;
    r->coeffs[4*i+0] |= (uint32_t)a[9*i+2] << 16;
    r->coeffs[4*i+0] &= 0x3FFFF;

    r->coeffs[4*i+1]  = a[9*i+2] >> 2;
    r->coeffs[4*i+1] |= (uint32_t)a[9*i+3] << 6;
    r->coeffs[4*i+1] |= (uint32_t)a[9*i+4] << 14;
    r->coeffs[4*i+1] &= 0x3FFFF;

    r->coeffs[4*i+2]  = a[9*i+4] >> 4;
    r->coeffs[4*i+2] |= (uint32_t)a[9*i+5] << 4;
    r->coeffs[4*i+2] |= (uint32_t)a[9*i+6] << 12;
    r->coeffs[4*i+2] &= 0x3FFFF;

    r->coeffs[4*i+3]  = a[9*i+6] >> 6;
    r->coeffs[4*i+3] |= (uint32_t)a[9*i+7] << 2;
    r->coeffs[4*i+3] |= (uint32_t)a[9*i+8] << 10;
    r->coeffs[4*i+3] &= 0x3FFFF;

    r->coeffs[4*i+0] = GAMMA1 - r->coeffs[4*i+0];
    r->coeffs[4*i+1] = GAMMA1 - r->coeffs[4*i+1];
    r->coeffs[4*i+2] = GAMMA1 - r->coeffs[4*i+2];
    r->coeffs[4*i+3] = GAMMA1 - r->coeffs[4*i+3];
  }
#elif GAMMA1 == (1 << 19)
  for(i = 0; i < N/2; ++i) {
    r->coeffs[2*i+0]  = a[5*i+0];
    r->coeffs[2*i+0] |= (uint32_t)a[5*i+1] << 8;
    r->coeffs[2*i+0] |= (uint32_t)a[5*i+2] << 16;
    r->coeffs[2*i+0] &= 0xFFFFF;

    r->coeffs[2*i+1]  = a[5*i+2] >> 4;
    r->coeffs[2*i+1] |= (uint32_t)a[5*i+3] << 4;
    r->coeffs[2*i+1] |= (uint32_t)a[5*i+4] << 12;
    /* r->coeffs[2*i+1] &= 0xFFFFF; */ /* No effect, since we're anyway at 20 bits */

    r->coeffs[2*i+0] = GAMMA1 - r->coeffs[2*i+0];
    r->coeffs[2*i+1] = GAMMA1 - r->coeffs[2*i+1];
  }
#endif

  DBENCH_STOP(*tpack);
}
/*************************************************
* Name:        sclls_polyz0_pack
*
* Description: Bit-pack polynomial z0 with coefficients in
*              [-2*(GAMMA1 - 1), 2*(GAMMA1)]. Packing uses a
*              little-endian bit stream.
*
*              - If GAMMA1 == (1<<17): pack 8 coefficients into 19 bytes
*                (each coefficient uses 19 bits).
*              - If GAMMA1 == (1<<19): pack 8 coefficients into 21 bytes
*                (each coefficient uses 21 bits).
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYZ_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void sclls_polyz0_pack(uint8_t *r, const poly *a) {
  unsigned int i;
  uint32_t t[8];
  DBENCH_START();

#if GAMMA1 == (1 << 17)
  for(i = 0; i < N/8; ++i) {
    t[0] = 2*GAMMA1 - a->coeffs[8*i+0];
    t[1] = 2*GAMMA1 - a->coeffs[8*i+1];
    t[2] = 2*GAMMA1 - a->coeffs[8*i+2];
    t[3] = 2*GAMMA1 - a->coeffs[8*i+3];
    t[4] = 2*GAMMA1 - a->coeffs[8*i+4];
    t[5] = 2*GAMMA1 - a->coeffs[8*i+5];
    t[6] = 2*GAMMA1 - a->coeffs[8*i+6];
    t[7] = 2*GAMMA1 - a->coeffs[8*i+7];

    /* Pack 8 coefficients, each 19 bits (little-endian bit stream):
     *   r[19*i + 0..18] stores bits stream[(19*i*8) + 8*j .. + 8*j+7].
     */
    r[19*i+0]  = (uint8_t)(t[0] & 0xFF);
    r[19*i+1]  = (uint8_t)((t[0] >> 8) & 0xFF);
    r[19*i+2]  = (uint8_t)(((t[0] >> 16) & 0x07) | ((t[1] & 0x1Fu) << 3));
    r[19*i+3]  = (uint8_t)((t[1] >> 5) & 0xFF);
    r[19*i+4]  = (uint8_t)(((t[1] >> 13) & 0x3Fu) | ((t[2] & 0x03u) << 6));
    r[19*i+5]  = (uint8_t)((t[2] >> 2) & 0xFF);
    r[19*i+6]  = (uint8_t)((t[2] >> 10) & 0xFF);
    r[19*i+7]  = (uint8_t)(((t[2] >> 18) & 0x01u) | ((t[3] & 0x7Fu) << 1));
    r[19*i+8]  = (uint8_t)((t[3] >> 7) & 0xFF);
    r[19*i+9]  = (uint8_t)(((t[3] >> 15) & 0x0Fu) | ((t[4] & 0x0Fu) << 4));
    r[19*i+10] = (uint8_t)((t[4] >> 4) & 0xFF);
    r[19*i+11] = (uint8_t)(((t[4] >> 12) & 0x7Fu) | ((t[5] & 0x01u) << 7));
    r[19*i+12] = (uint8_t)((t[5] >> 1) & 0xFF);
    r[19*i+13] = (uint8_t)((t[5] >> 9) & 0xFF);
    r[19*i+14] = (uint8_t)(((t[5] >> 17) & 0x03u) | ((t[6] & 0x3Fu) << 2));
    r[19*i+15] = (uint8_t)((t[6] >> 6) & 0xFF);
    r[19*i+16] = (uint8_t)(((t[6] >> 14) & 0x1Fu) | ((t[7] & 0x07u) << 5));
    r[19*i+17] = (uint8_t)((t[7] >> 3) & 0xFF);
    r[19*i+18] = (uint8_t)((t[7] >> 11) & 0xFF);
  }
#elif GAMMA1 == (1 << 19)
  for(i = 0; i < N/8; ++i) {
    t[0] = 2*GAMMA1 - a->coeffs[8*i+0];
    t[1] = 2*GAMMA1 - a->coeffs[8*i+1];
    t[2] = 2*GAMMA1 - a->coeffs[8*i+2];
    t[3] = 2*GAMMA1 - a->coeffs[8*i+3];
    t[4] = 2*GAMMA1 - a->coeffs[8*i+4];
    t[5] = 2*GAMMA1 - a->coeffs[8*i+5];
    t[6] = 2*GAMMA1 - a->coeffs[8*i+6];
    t[7] = 2*GAMMA1 - a->coeffs[8*i+7];

    /* Pack 8 coefficients, each 21 bits (little-endian bit stream).
     * Byte layout (least significant bits first) across t[0..7]:
     *   r[21*i + 0..20] = stream bits [0..167] grouped into 8-bit chunks.
     */
    r[21*i+0]  = (uint8_t)(t[0] & 0xFF);
    r[21*i+1]  = (uint8_t)((t[0] >> 8) & 0xFF);
    r[21*i+2]  = (uint8_t)(((t[0] >> 16) | (t[1] << 5)) & 0xFF);
    r[21*i+3]  = (uint8_t)((t[1] >> 3) & 0xFF);
    r[21*i+4]  = (uint8_t)((t[1] >> 11) & 0xFF);
    r[21*i+5]  = (uint8_t)(((t[1] >> 19) | (t[2] << 2)) & 0xFF);
    r[21*i+6]  = (uint8_t)((t[2] >> 6) & 0xFF);
    r[21*i+7]  = (uint8_t)(((t[2] >> 14) | (t[3] << 7)) & 0xFF);
    r[21*i+8]  = (uint8_t)((t[3] >> 1) & 0xFF);
    r[21*i+9]  = (uint8_t)((t[3] >> 9) & 0xFF);
    r[21*i+10] = (uint8_t)(((t[3] >> 17) | (t[4] << 4)) & 0xFF);
    r[21*i+11] = (uint8_t)((t[4] >> 4) & 0xFF);
    r[21*i+12] = (uint8_t)((t[4] >> 12) & 0xFF);
    r[21*i+13] = (uint8_t)(((t[4] >> 20) | (t[5] << 1)) & 0xFF);
    r[21*i+14] = (uint8_t)((t[5] >> 7) & 0xFF);
    r[21*i+15] = (uint8_t)(((t[5] >> 15) | (t[6] << 6)) & 0xFF);
    r[21*i+16] = (uint8_t)((t[6] >> 2) & 0xFF);
    r[21*i+17] = (uint8_t)((t[6] >> 10) & 0xFF);
    r[21*i+18] = (uint8_t)(((t[6] >> 18) | (t[7] << 3)) & 0xFF);
    r[21*i+19] = (uint8_t)((t[7] >> 5) & 0xFF);
    r[21*i+20] = (uint8_t)((t[7] >> 13) & 0xFF);

  }
#endif

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        sclls_polyz0_unpack
*
* Description: Unpack polynomial z0 with coefficients in
*              [-2*(GAMMA1 - 1), 2*(GAMMA1)] from the little-endian
*              bit stream produced by sclls_polyz0_pack.
*
*              - If GAMMA1 == (1<<17): 8 coefficients are decoded from 19 bytes
*                (each coefficient uses 19 bits).
*              - If GAMMA1 == (1<<19): 8 coefficients are decoded from 21 bytes
*                (each coefficient uses 21 bits).
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void sclls_polyz0_unpack(poly *r, const uint8_t *a) {
  unsigned int i;
  DBENCH_START();

  const uint32_t two_gamma = (uint32_t)(2*GAMMA1);

#if GAMMA1 == (1 << 17)
  /* Inverse of sclls_polyz0_pack for GAMMA1 == 2^17.
   * Pack: 8 coefficients -> 19 bytes, 19-bit little-endian bit stream.
   */
  for(i = 0; i < N/8; ++i) {
    uint32_t b0  = a[19*i+0];
    uint32_t b1  = a[19*i+1];
    uint32_t b2  = a[19*i+2];
    uint32_t b3  = a[19*i+3];
    uint32_t b4  = a[19*i+4];
    uint32_t b5  = a[19*i+5];
    uint32_t b6  = a[19*i+6];
    uint32_t b7  = a[19*i+7];
    uint32_t b8  = a[19*i+8];
    uint32_t b9  = a[19*i+9];
    uint32_t b10 = a[19*i+10];
    uint32_t b11 = a[19*i+11];
    uint32_t b12 = a[19*i+12];
    uint32_t b13 = a[19*i+13];
    uint32_t b14 = a[19*i+14];
    uint32_t b15 = a[19*i+15];
    uint32_t b16 = a[19*i+16];
    uint32_t b17 = a[19*i+17];
    uint32_t b18 = a[19*i+18];

    /* Reconstruct t[0..7] (each 19 bits) */
    uint32_t t0 = b0 | (b1 << 8) | ((b2 & 0x07) << 16);
    uint32_t t1 = (b2 >> 3) | (b3 << 5) | ((b4 & 0x3Fu) << 13);
    uint32_t t2 = ((b4 >> 6) & 0x03u) | (b5 << 2) | (b6 << 10) | ((b7 & 0x01u) << 18);
    uint32_t t3 = ((b7 >> 1) & 0x7Fu) | (b8 << 7) | ((b9 & 0x0Fu) << 15);
    uint32_t t4 = (b9 >> 4) | (b10 << 4) | ((b11 & 0x7Fu) << 12);
    uint32_t t5 = ((b11 >> 7) & 0x01u) | (b12 << 1) | (b13 << 9) | ((b14 & 0x03u) << 17);
    uint32_t t6 = ((b14 >> 2) & 0x3Fu) | (b15 << 6) | ((b16 & 0x1Fu) << 14);
    uint32_t t7 = ((b16 >> 5) & 0x07u) | (b17 << 3) | (b18 << 11);

    r->coeffs[8*i+0] = (int32_t)(two_gamma - t0);
    r->coeffs[8*i+1] = (int32_t)(two_gamma - t1);
    r->coeffs[8*i+2] = (int32_t)(two_gamma - t2);
    r->coeffs[8*i+3] = (int32_t)(two_gamma - t3);
    r->coeffs[8*i+4] = (int32_t)(two_gamma - t4);
    r->coeffs[8*i+5] = (int32_t)(two_gamma - t5);
    r->coeffs[8*i+6] = (int32_t)(two_gamma - t6);
    r->coeffs[8*i+7] = (int32_t)(two_gamma - t7);
  }
#elif GAMMA1 == (1 << 19)
  /* Inverse of sclls_polyz0_pack for GAMMA1 == 2^19.
   * Pack: 8 coefficients -> 21 bytes, 21-bit little-endian bit stream.
   */
  for(i = 0; i < N/8; ++i) {
    uint32_t b0  = a[21*i+0];
    uint32_t b1  = a[21*i+1];
    uint32_t b2  = a[21*i+2];
    uint32_t b3  = a[21*i+3];
    uint32_t b4  = a[21*i+4];
    uint32_t b5  = a[21*i+5];
    uint32_t b6  = a[21*i+6];
    uint32_t b7  = a[21*i+7];
    uint32_t b8  = a[21*i+8];
    uint32_t b9  = a[21*i+9];
    uint32_t b10 = a[21*i+10];
    uint32_t b11 = a[21*i+11];
    uint32_t b12 = a[21*i+12];
    uint32_t b13 = a[21*i+13];
    uint32_t b14 = a[21*i+14];
    uint32_t b15 = a[21*i+15];
    uint32_t b16 = a[21*i+16];
    uint32_t b17 = a[21*i+17];
    uint32_t b18 = a[21*i+18];
    uint32_t b19 = a[21*i+19];
    uint32_t b20 = a[21*i+20];

    /* Reconstruct t[0..7] (each 21 bits) */
    uint32_t t0 = b0 | (b1 << 8) | ((b2 & 0x1Fu) << 16);
    uint32_t t1 = (b2 >> 5) | (b3 << 3) | (b4 << 11) | ((b5 & 0x03u) << 19);
    uint32_t t2 = (b5 >> 2) | (b6 << 6) | ((b7 & 0x7Fu) << 14);
    uint32_t t3 = (b7 >> 7) | (b8 << 1) | (b9 << 9) | ((b10 & 0x0Fu) << 17);
    uint32_t t4 = (b10 >> 4) | (b11 << 4) | (b12 << 12) | ((b13 & 0x01u) << 20);
    uint32_t t5 = ((b13 >> 1) & 0x7Fu) | (b14 << 7) | ((b15 & 0x3Fu) << 15);
    uint32_t t6 = ((b15 >> 6) & 0x03u) | (b16 << 2) | (b17 << 10) | ((b18 & 0x07u) << 18);
    uint32_t t7 = ((b18 >> 3) & 0x1Fu) | (b19 << 5) | (b20 << 13);

    r->coeffs[8*i+0] = (int32_t)(two_gamma - t0);
    r->coeffs[8*i+1] = (int32_t)(two_gamma - t1);
    r->coeffs[8*i+2] = (int32_t)(two_gamma - t2);
    r->coeffs[8*i+3] = (int32_t)(two_gamma - t3);
    r->coeffs[8*i+4] = (int32_t)(two_gamma - t4);
    r->coeffs[8*i+5] = (int32_t)(two_gamma - t5);
    r->coeffs[8*i+6] = (int32_t)(two_gamma - t6);
    r->coeffs[8*i+7] = (int32_t)(two_gamma - t7);
  }
#endif

  DBENCH_STOP(*tpack);
}
/*************************************************
* Name:        clps_polyz0a_pack
*
* Description: Pack polynomial z_0_A with coefficients in
*              (-CLPS_Z0A_BOUND, CLPS_Z0A_BOUND), CLPS_Z0A_BOUND = 3*(GAMMA1-BETA).
*              Each coefficient is stored as coeff + CLPS_Z0A_BOUND using
*              22 bits (little-endian bit stream); 8 coefficients -> 22 bytes.
*
* Arguments:   - uint8_t *r: output, length CLPS_POLYZ0A_PACKEDBYTES
*              - const poly *a: pointer to input polynomial
**************************************************/
void clps_polyz0a_pack(uint8_t *r, const poly *a) {
  unsigned int i;
  uint32_t t[8];
  DBENCH_START();

  for(i = 0; i < N/8; ++i) {
    t[0] = (uint32_t)(a->coeffs[8*i+0] + CLPS_Z0A_BOUND);
    t[1] = (uint32_t)(a->coeffs[8*i+1] + CLPS_Z0A_BOUND);
    t[2] = (uint32_t)(a->coeffs[8*i+2] + CLPS_Z0A_BOUND);
    t[3] = (uint32_t)(a->coeffs[8*i+3] + CLPS_Z0A_BOUND);
    t[4] = (uint32_t)(a->coeffs[8*i+4] + CLPS_Z0A_BOUND);
    t[5] = (uint32_t)(a->coeffs[8*i+5] + CLPS_Z0A_BOUND);
    t[6] = (uint32_t)(a->coeffs[8*i+6] + CLPS_Z0A_BOUND);
    t[7] = (uint32_t)(a->coeffs[8*i+7] + CLPS_Z0A_BOUND);

    /* Pack 8 coefficients, each 22 bits (little-endian bit stream):
     * stream bit b = coefficient (b/22), bit (b%22).
     */
    r[22*i+0]  = (uint8_t)(t[0] & 0xFF);
    r[22*i+1]  = (uint8_t)((t[0] >> 8) & 0xFF);
    r[22*i+2]  = (uint8_t)(((t[0] >> 16) & 0x3Fu) | ((t[1] & 0x03u) << 6));
    r[22*i+3]  = (uint8_t)((t[1] >> 2) & 0xFF);
    r[22*i+4]  = (uint8_t)((t[1] >> 10) & 0xFF);
    r[22*i+5]  = (uint8_t)(((t[1] >> 18) & 0x0Fu) | ((t[2] & 0x0Fu) << 4));
    r[22*i+6]  = (uint8_t)((t[2] >> 4) & 0xFF);
    r[22*i+7]  = (uint8_t)((t[2] >> 12) & 0xFF);
    r[22*i+8]  = (uint8_t)(((t[2] >> 20) & 0x03u) | ((t[3] & 0x3Fu) << 2));
    r[22*i+9]  = (uint8_t)((t[3] >> 6) & 0xFF);
    r[22*i+10] = (uint8_t)((t[3] >> 14) & 0xFF);
    r[22*i+11] = (uint8_t)(t[4] & 0xFF);
    r[22*i+12] = (uint8_t)((t[4] >> 8) & 0xFF);
    r[22*i+13] = (uint8_t)(((t[4] >> 16) & 0x3Fu) | ((t[5] & 0x03u) << 6));
    r[22*i+14] = (uint8_t)((t[5] >> 2) & 0xFF);
    r[22*i+15] = (uint8_t)((t[5] >> 10) & 0xFF);
    r[22*i+16] = (uint8_t)(((t[5] >> 18) & 0x0Fu) | ((t[6] & 0x0Fu) << 4));
    r[22*i+17] = (uint8_t)((t[6] >> 4) & 0xFF);
    r[22*i+18] = (uint8_t)((t[6] >> 12) & 0xFF);
    r[22*i+19] = (uint8_t)(((t[6] >> 20) & 0x03u) | ((t[7] & 0x3Fu) << 2));
    r[22*i+20] = (uint8_t)((t[7] >> 6) & 0xFF);
    r[22*i+21] = (uint8_t)((t[7] >> 14) & 0xFF);
  }

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        clps_polyz0a_unpack
*
* Description: Inverse of clps_polyz0a_pack. Each coefficient is decoded
*              from 22 bits and shifted back by CLPS_Z0A_BOUND.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: packed data, length CLPS_POLYZ0A_PACKEDBYTES
**************************************************/
void clps_polyz0a_unpack(poly *r, const uint8_t *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N/8; ++i) {
    uint32_t b0  = a[22*i+0];
    uint32_t b1  = a[22*i+1];
    uint32_t b2  = a[22*i+2];
    uint32_t b3  = a[22*i+3];
    uint32_t b4  = a[22*i+4];
    uint32_t b5  = a[22*i+5];
    uint32_t b6  = a[22*i+6];
    uint32_t b7  = a[22*i+7];
    uint32_t b8  = a[22*i+8];
    uint32_t b9  = a[22*i+9];
    uint32_t b10 = a[22*i+10];
    uint32_t b11 = a[22*i+11];
    uint32_t b12 = a[22*i+12];
    uint32_t b13 = a[22*i+13];
    uint32_t b14 = a[22*i+14];
    uint32_t b15 = a[22*i+15];
    uint32_t b16 = a[22*i+16];
    uint32_t b17 = a[22*i+17];
    uint32_t b18 = a[22*i+18];
    uint32_t b19 = a[22*i+19];
    uint32_t b20 = a[22*i+20];
    uint32_t b21 = a[22*i+21];

    /* Reconstruct t[0..7] (each 22 bits) */
    uint32_t t0 = b0 | (b1 << 8) | ((b2 & 0x3Fu) << 16);
    uint32_t t1 = (b2 >> 6) | (b3 << 2) | (b4 << 10) | ((b5 & 0x0Fu) << 18);
    uint32_t t2 = (b5 >> 4) | (b6 << 4) | (b7 << 12) | ((b8 & 0x03u) << 20);
    uint32_t t3 = (b8 >> 2) | (b9 << 6) | (b10 << 14);
    uint32_t t4 = b11 | (b12 << 8) | ((b13 & 0x3Fu) << 16);
    uint32_t t5 = (b13 >> 6) | (b14 << 2) | (b15 << 10) | ((b16 & 0x0Fu) << 18);
    uint32_t t6 = (b16 >> 4) | (b17 << 4) | (b18 << 12) | ((b19 & 0x03u) << 20);
    uint32_t t7 = (b19 >> 2) | (b20 << 6) | (b21 << 14);

    r->coeffs[8*i+0] = (int32_t)(t0 - CLPS_Z0A_BOUND);
    r->coeffs[8*i+1] = (int32_t)(t1 - CLPS_Z0A_BOUND);
    r->coeffs[8*i+2] = (int32_t)(t2 - CLPS_Z0A_BOUND);
    r->coeffs[8*i+3] = (int32_t)(t3 - CLPS_Z0A_BOUND);
    r->coeffs[8*i+4] = (int32_t)(t4 - CLPS_Z0A_BOUND);
    r->coeffs[8*i+5] = (int32_t)(t5 - CLPS_Z0A_BOUND);
    r->coeffs[8*i+6] = (int32_t)(t6 - CLPS_Z0A_BOUND);
    r->coeffs[8*i+7] = (int32_t)(t7 - CLPS_Z0A_BOUND);
  }

  DBENCH_STOP(*tpack);
}
/*************************************************
* Name:        sclls_polyhw_pack
*
* Description: Map each allowed coefficient to a 4-bit codeword 0..11, then
*              pack two codewords per byte (N/2 bytes per polynomial).
*              Allowed sets (hw = HighBits(a+b+c) - HighBits(a) - HighBits(b)
*              - HighBits(c), with the modulus wrap counted by m=0,1,2):
*              GAMMA1==(1<<17): {-90..-87,-46..-43,-2..1} -> {0..3,4..7,8..11}
*              GAMMA1==(1<<19): {-34..-31,-18..-15,-2..1} -> {0..3,4..7,8..11}
*
* Arguments:   - uint8_t *r: output, length SCLLS_HWBYTES_PERPOLY (256*4/8)
*              - const poly *a: input polynomial
**************************************************/
void sclls_polyhw_pack(uint8_t *r, const poly *a) {
  unsigned int i;
  unsigned int c0, c1;
  DBENCH_START();
  for(i = 0; i < N/2; ++i) {
#if GAMMA1 == (1 << 17)
    if(a->coeffs[2*i+0] <= -87) c0 = (unsigned int)(a->coeffs[2*i+0] + 90);
    else if(a->coeffs[2*i+0] <= -43) c0 = (unsigned int)(a->coeffs[2*i+0] + 50);
    else c0 = (unsigned int)(a->coeffs[2*i+0] + 10);
    if(a->coeffs[2*i+1] <= -87) c1 = (unsigned int)(a->coeffs[2*i+1] + 90);
    else if(a->coeffs[2*i+1] <= -43) c1 = (unsigned int)(a->coeffs[2*i+1] + 50);
    else c1 = (unsigned int)(a->coeffs[2*i+1] + 10);
#elif GAMMA1 == (1 << 19)
    if(a->coeffs[2*i+0] <= -31) c0 = (unsigned int)(a->coeffs[2*i+0] + 34);
    else if(a->coeffs[2*i+0] <= -15) c0 = (unsigned int)(a->coeffs[2*i+0] + 22);
    else c0 = (unsigned int)(a->coeffs[2*i+0] + 10);
    if(a->coeffs[2*i+1] <= -31) c1 = (unsigned int)(a->coeffs[2*i+1] + 34);
    else if(a->coeffs[2*i+1] <= -15) c1 = (unsigned int)(a->coeffs[2*i+1] + 22);
    else c1 = (unsigned int)(a->coeffs[2*i+1] + 10);
#endif
    r[i] = (uint8_t)((c0 & 0x0F) | ((c1 & 0x0F) << 4));
  }
  DBENCH_STOP(*tpack);
}
/*************************************************
* Name:        sclls_polyhw_unpack
*
* Description: Inverse of sclls_polyhw_pack (4-bit codewords, two per byte).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t *r: packed data, length SCLLS_HWBYTES_PERPOLY
**************************************************/
void sclls_polyhw_unpack(poly *a, const uint8_t *r) {
  unsigned int i;
  uint8_t b, c0, c1;
  DBENCH_START();

  for(i = 0; i < N/2; ++i) {
    b = r[i];
    c0 = b & 0x0F;
    c1 = b >> 4;
#if GAMMA1 == (1 << 17)
    a->coeffs[2*i+0] = (c0 <= 3) ? (int32_t)c0 - 90 : ((c0 <= 7) ? (int32_t)c0 - 50 : (int32_t)c0 - 10);
    a->coeffs[2*i+1] = (c1 <= 3) ? (int32_t)c1 - 90 : ((c1 <= 7) ? (int32_t)c1 - 50 : (int32_t)c1 - 10);
#elif GAMMA1 == (1 << 19)
    a->coeffs[2*i+0] = (c0 <= 3) ? (int32_t)c0 - 34 : ((c0 <= 7) ? (int32_t)c0 - 22 : (int32_t)c0 - 10);
    a->coeffs[2*i+1] = (c1 <= 3) ? (int32_t)c1 - 34 : ((c1 <= 7) ? (int32_t)c1 - 22 : (int32_t)c1 - 10);
#endif
  }
  DBENCH_STOP(*tpack);
}
/*************************************************
* Name:        polyw1_pack
*
* Description: Bit-pack polynomial w1 with coefficients in [0,15] or [0,43].
*              Input coefficients are assumed to be standard representatives.
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYW1_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyw1_pack(uint8_t *r, const poly *a) {
  unsigned int i;
  DBENCH_START();
  /* 每系数6 bit，4系数打包为3字节 (GAMMA2==(Q-1)/88) */
#if GAMMA2 == (Q-1)/88
  for(i = 0; i < N/4; ++i) {
    r[3*i+0]  = a->coeffs[4*i+0];
    r[3*i+0] |= a->coeffs[4*i+1] << 6;
    r[3*i+1]  = a->coeffs[4*i+1] >> 2;
    r[3*i+1] |= a->coeffs[4*i+2] << 4;
    r[3*i+2]  = a->coeffs[4*i+2] >> 4;
    r[3*i+2] |= a->coeffs[4*i+3] << 2;
  }
  /* 每系数4 bit，2系数1字节 (GAMMA2==(Q-1)/32) */
#elif GAMMA2 == (Q-1)/32
  /* 系数为 [0,15] 的模16表示: 双方 w1 的整数可能相差 16(模数),
     必须对低位同样掩码, 否则 bit4 会泄漏进打包字节, 与 unpack 不对称 */
  for(i = 0; i < N/2; ++i)
    r[i] = (a->coeffs[2*i+0] & 0x0F) | ((a->coeffs[2*i+1] & 0x0F) << 4);
#endif

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        sclls_polyw1_unpack
*
* Description: Unpack polynomial w1 from bit-packed bytes.
*              Output coefficients are in [0,15] or [0,43]
*              depending on GAMMA2.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with at least
*                                  POLYW1_PACKEDBYTES bytes
**************************************************/
void sclls_polyw1_unpack(poly *r, const uint8_t *a) {
  unsigned int i;
  DBENCH_START();

#if GAMMA2 == (Q-1)/88
  for(i = 0; i < N/4; ++i) {
    r->coeffs[4*i+0] = a[3*i+0] & 0x3F;
    r->coeffs[4*i+1] = (a[3*i+0] >> 6) | ((a[3*i+1] & 0x0F) << 2);
    r->coeffs[4*i+2] = (a[3*i+1] >> 4) | ((a[3*i+2] & 0x03) << 4);
    r->coeffs[4*i+3] = (a[3*i+2] >> 2) & 0x3F;
  }
#elif GAMMA2 == (Q-1)/32
  for(i = 0; i < N/2; ++i) {
    r->coeffs[2*i+0] = a[i] & 0x0F;
    r->coeffs[2*i+1] = (a[i] >> 4) & 0x0F;
  }
#endif

  DBENCH_STOP(*tpack);
}
