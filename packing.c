#include "params.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"

/*************************************************
* Name:        sclls_pack_params
*
* Description: Bit-pack params = (rho, t).
*
* Arguments:   - uint8_t params[]: output byte array
*              - const uint8_t rho[]: byte array containing rho
*              - const polyveck *t: pointer to vector t
**************************************************/
void sclls_pack_params(uint8_t params[SCLLS_PARAMSBYTES],
                       const uint8_t rho[SEEDBYTES],
                       const polyveck *t)
{
  unsigned int i;

  for (i = 0; i < SEEDBYTES; ++i)
    params[i] = rho[i];
  params += SEEDBYTES;

  for (i = 0; i < K; ++i)
    sclls_polyt_pack(params + i*SCLLS_POLYT_PACKEDBYTES, &t->vec[i]);
}

/*************************************************
* Name:        sclls_unpack_params
*
* Description: Unpack params = (rho, t).
*
* Arguments:   - uint8_t rho[]: output byte array for rho
*              - polyveck *t: pointer to output vector t
*              - const uint8_t params[]: byte array containing bit-packed params
**************************************************/
void sclls_unpack_params(uint8_t rho[SEEDBYTES],
                         polyveck *t,
                         const uint8_t params[SCLLS_PARAMSBYTES])
{
  unsigned int i;

  for (i = 0; i < SEEDBYTES; ++i)
    rho[i] = params[i];
  params += SEEDBYTES;

  for (i = 0; i < K; ++i)
    sclls_polyt_unpack(&t->vec[i], params + i*SCLLS_POLYT_PACKEDBYTES);
}

/*************************************************
* Name:        sclls_pack_msk
*
* Description: Bit-pack secret key msk = (key, tr, s1, s2).
*
* Arguments:   - uint8_t msk[]: output byte array
*              - const uint8_t key[]: byte array containing key
*              - const uint8_t tr[]: byte array containing tr
*              - const polyvecl *s1: pointer to vector s1
*              - const polyveck *s2: pointer to vector s2
**************************************************/
void sclls_pack_msk(uint8_t msk[SCLLS_MSKBYTES],
                    const uint8_t key[SEEDBYTES],
                    const uint8_t tr[TRBYTES],
                    const polyvecl *s1,
                    const polyveck *s2)
{
  unsigned int i;

  for (i = 0; i < SEEDBYTES; ++i)
    msk[i] = key[i];
  msk += SEEDBYTES;

  for (i = 0; i < TRBYTES; ++i)
    msk[i] = tr[i];
  msk += TRBYTES;

  for (i = 0; i < L; ++i)
    polyeta_pack(msk + i*POLYETA_PACKEDBYTES, &s1->vec[i]);
  msk += L*POLYETA_PACKEDBYTES;

  for (i = 0; i < K; ++i)
    polyeta_pack(msk + i*POLYETA_PACKEDBYTES, &s2->vec[i]);
  msk += K*POLYETA_PACKEDBYTES;
}

/*************************************************
* Name:        sclls_unpack_msk
*
* Description: Unpack secret key msk = (key, tr, s1, s2).
*
* Arguments:   - uint8_t key[]: output byte array for key
*              - uint8_t tr[]: output byte array for tr
*              - polyvecl *s1: pointer to output vector s1
*              - polyveck *s2: pointer to output vector s2
*              - const uint8_t msk[]: byte array containing bit-packed msk
**************************************************/
void sclls_unpack_msk(uint8_t key[SEEDBYTES],
                      uint8_t tr[TRBYTES],
                      polyvecl *s1,
                      polyveck *s2,
                      const uint8_t msk[SCLLS_MSKBYTES])
{
  unsigned int i;
  const uint8_t *p = msk;

  for (i = 0; i < SEEDBYTES; ++i)
    key[i] = p[i];
  p += SEEDBYTES;

  for (i = 0; i < TRBYTES; ++i)
    tr[i] = p[i];
  p += TRBYTES;

  for (i = 0; i < L; ++i)
    polyeta_unpack(&s1->vec[i], p + i*POLYETA_PACKEDBYTES);
  p += L*POLYETA_PACKEDBYTES;

  for (i = 0; i < K; ++i)
    polyeta_unpack(&s2->vec[i], p + i*POLYETA_PACKEDBYTES);
}
/*************************************************
* Name:        sclls_pack_sig
*
* Description: Bit-pack SCLLS signature fields. On-wire layout:
*              cid (CTILDEBYTES) || z0 (L * POLYZ_PACKEDBYTES) ||
*              hw (K * SCLLS_HWBYTES).
*
* Arguments:   - uint8_t sig[]: output byte array of length SCLLS_SIGBYTES
*              - const uint8_t cid[]: challenge hash, length CTILDEBYTES
*              - const polyvecl *z0: pointer to z0 vector
*              - const polyveck *hw: pointer to hw vector
**************************************************/
void sclls_pack_did(uint8_t did[SCLLS_DIDBYTES],
                    const uint8_t cid[CTILDEBYTES],
                    const polyvecl *zid)
{
  unsigned int i;
  for (i = 0; i < CTILDEBYTES; ++i)
    did[i] = cid[i];
  did += CTILDEBYTES;
  for (i = 0; i < L; ++i)
    polyz_pack(did + i*POLYZ_PACKEDBYTES, &zid->vec[i]);
}

/*************************************************
* Name:        sclls_unpack_did
*
* Description: Unpack DID from sclls_pack_did (cid then zid).
*
* Arguments:   - uint8_t cid[]: output challenge hash, length CTILDEBYTES
*              - polyvecl *zid: pointer to output vector z
*              - const uint8_t did[]: byte array of length SCLLS_DIDBYTES
**************************************************/
void sclls_unpack_did(uint8_t cid[CTILDEBYTES],
                      polyvecl *zid,
                      const uint8_t did[SCLLS_DIDBYTES])
{
  unsigned int i;
  const uint8_t *p = did;

  for (i = 0; i < CTILDEBYTES; ++i)
    cid[i] = p[i];
  p += CTILDEBYTES;
  for (i = 0; i < L; ++i)
    polyz_unpack(&zid->vec[i], p + i*POLYZ_PACKEDBYTES);
}
/*************************************************
* Name:        sclls_pack_pid
*
* Description: Bit-pack PID = (t, w1buf).
*              Layout: K polynomials t[i] via sclls_polyt_pack (each
*              SCLLS_POLYT_PACKEDBYTES), then K*w1-packed polynomials
*              concatenated as in polyw1_pack (total K*POLYW1_PACKEDBYTES).
*              Inverse of sclls_unpack_pid on the same buffers.
*
* Arguments:   - uint8_t pid[]: output, length SCLLS_PIDBYTES
*              - const polyveck *t: K-vector t (coeffs in [0,Q) per sclls_polyt_pack)
*              - const uint8_t w1buf[]: packed w1 hints, length K*POLYW1_PACKEDBYTES
**************************************************/
void sclls_pack_pid(uint8_t pid[SCLLS_PIDBYTES],
                    const polyveck *t,
                    const uint8_t w1buf[K*POLYW1_PACKEDBYTES])
{
  unsigned int i;
  for (i = 0; i < K; ++i)
    sclls_polyt_pack(pid + i*SCLLS_POLYT_PACKEDBYTES, &t->vec[i]);
  pid += K*SCLLS_POLYT_PACKEDBYTES;
  for (i = 0; i < K*POLYW1_PACKEDBYTES; ++i)
    pid[i] = w1buf[i];
}

/*************************************************
* Name:        sclls_unpack_pid
*
* Description: Unpack PID from sclls_pack_pid (exact inverse: t via
*              sclls_polyt_unpack, w1buf byte-identical).
*
* Arguments:   - uint8_t w1buf[]: output w1 hints, length K*POLYW1_PACKEDBYTES
*              - polyveck *t: output K-vector t
*              - const uint8_t pid[]: input, length SCLLS_PIDBYTES
**************************************************/
void sclls_unpack_pid(uint8_t w1buf[K*POLYW1_PACKEDBYTES],
                      polyveck *t,
                      const uint8_t pid[SCLLS_PIDBYTES])
{
  unsigned int i;
  for (i = 0; i < K; ++i)
    sclls_polyt_unpack(&t->vec[i], pid + i*SCLLS_POLYT_PACKEDBYTES);
  pid += K*SCLLS_POLYT_PACKEDBYTES;
  for (i = 0; i < K*POLYW1_PACKEDBYTES; ++i)
    w1buf[i] = pid[i];
}
/*************************************************
* Name:        sclls_pack_sid
*
* Description: Bit-pack SID = (s1, s2, cid, zid).
*              Layout: L times polyeta_pack(s1[i]); K times polyeta_pack(s2[i]);
*              then raw cid (SCLLS_CIDBYTES), then raw zid (L*POLYZ_PACKEDBYTES).
*              Inverse of sclls_unpack_sid.
*
* Arguments:   - uint8_t sid[]: output, length SCLLS_SIDBYTES
*              - const polyvecl *s1, const polyveck *s2: secret vectors
*              - const uint8_t cid[]: challenge bytes (length SCLLS_CIDBYTES)
*              - const uint8_t zid[]: packed z (length SCLLS_ZIDBYTES)
**************************************************/
void sclls_pack_sid(uint8_t sid[SCLLS_SIDBYTES],
                    const polyvecl *s1,
                    const polyveck *s2,
                    const uint8_t cid[SCLLS_CIDBYTES],                    
                    const uint8_t zid[SCLLS_ZIDBYTES])
{
  unsigned int i;
  for (i = 0; i < L; ++i)
    polyeta_pack(sid + i*POLYETA_PACKEDBYTES, &s1->vec[i]);
  sid += L*POLYETA_PACKEDBYTES;
  for (i = 0; i < K; ++i)
    polyeta_pack(sid + i*POLYETA_PACKEDBYTES, &s2->vec[i]);
  sid += K*POLYETA_PACKEDBYTES;
  for (i = 0; i < SCLLS_CIDBYTES; ++i)
    sid[i] = cid[i];
  sid += SCLLS_CIDBYTES;
  for (i = 0; i < SCLLS_ZIDBYTES; ++i)
    sid[i] = zid[i];
}
/*************************************************
* Name:        sclls_unpack_sid
*
* Description: Unpack SID from sclls_pack_sid (exact inverse: s1/s2 via
*              polyeta_unpack; cid and zid byte-identical).
*
* Arguments:   - polyvecl *s1, polyveck *s2: output vectors
*              - uint8_t cid[], zid[]: output buffers
*              - const uint8_t sid[]: input, length SCLLS_SIDBYTES
**************************************************/
void sclls_unpack_sid(polyvecl *s1,
                      polyveck *s2,
                      uint8_t cid[SCLLS_CIDBYTES],
                      uint8_t zid[SCLLS_ZIDBYTES],
                      const uint8_t sid[SCLLS_SIDBYTES])
{
  unsigned int i;
  const uint8_t *p = sid;
  for (i = 0; i < L; ++i)
    polyeta_unpack(&s1->vec[i], p + i*POLYETA_PACKEDBYTES);
  p += L*POLYETA_PACKEDBYTES;
  for (i = 0; i < K; ++i)
    polyeta_unpack(&s2->vec[i], p + i*POLYETA_PACKEDBYTES);
  p += K*POLYETA_PACKEDBYTES;
  for (i = 0; i < SCLLS_CIDBYTES; ++i)
    cid[i] = p[i];
  p += SCLLS_CIDBYTES;
  for (i = 0; i < SCLLS_ZIDBYTES; ++i)
    zid[i] = p[i];
}
/*************************************************
* Name:        sclls_pack_sig
*
* Description: Bit-pack signature sig = (cid, z0, hw).
*              On-wire: cid[CTILDEBYTES] || z0[L] || hw[K].
*
* Arguments:   - uint8_t sig[]: output byte array of length SCLLS_SIGBYTES
*              - const uint8_t cid[]: challenge hash, length CTILDEBYTES
*              - const polyvecl *z0: vector z0 (packed with polyz_pack)
*              - const polyveck *hw: vector hw (packed with sclls_polyhw_pack)
**************************************************/
void sclls_pack_sig(uint8_t sig[SCLLS_SIGBYTES],
                    const uint8_t cid[CTILDEBYTES],
                    const polyvecl *z0,
                    const polyveck *hw)
{
  unsigned int i;
  for (i = 0; i < CTILDEBYTES; ++i)
    sig[i] = cid[i];
  sig += CTILDEBYTES;
  for (i = 0; i < L; ++i)
    sclls_polyz0_pack(sig + i*SCLLS_POLYZ0_PACKEDBYTES, &z0->vec[i]);
  sig += L*SCLLS_POLYZ0_PACKEDBYTES;
  for (i = 0; i < K; ++i)
    sclls_polyhw_pack(sig + i*SCLLS_HWBYTES_PERPOLY, &hw->vec[i]);
}

/*************************************************
* Name:        sclls_unpack_sig
*
* Description: Unpack signature sig from sclls_pack_sig.
*              On-wire: cid[CTILDEBYTES] || z0[L] || hw[K].
*
* Arguments:   - uint8_t cid[]: output challenge hash, length CTILDEBYTES
*              - polyvecl *z0: output vector z0
*              - polyveck *hw: output vector hw
*              - const uint8_t sig[]: input byte array of length SCLLS_SIGBYTES
**************************************************/
void sclls_unpack_sig(uint8_t cid[CTILDEBYTES],
                      polyvecl *z0,
                      polyveck *hw,
                      const uint8_t sig[SCLLS_SIGBYTES])
{
  unsigned int i;
  const uint8_t *p = sig;

  for (i = 0; i < CTILDEBYTES; ++i)
    cid[i] = p[i];
  p += CTILDEBYTES;

  for (i = 0; i < L; ++i)
    sclls_polyz0_unpack(&z0->vec[i], p + i*SCLLS_POLYZ0_PACKEDBYTES);
  p += L*SCLLS_POLYZ0_PACKEDBYTES;

  for (i = 0; i < K; ++i)
    sclls_polyhw_unpack(&hw->vec[i], p + i*SCLLS_HWBYTES_PERPOLY);
}
/*************************************************
* Name:        pack_pk
*
* Description: Bit-pack public key pk = (rho, t1).
*
* Arguments:   - uint8_t pk[]: output byte array
*              - const uint8_t rho[]: byte array containing rho
*              - const polyveck *t1: pointer to vector t1
**************************************************/
void pack_pk(uint8_t pk[CRYPTO_PUBLICKEYBYTES],
             const uint8_t rho[SEEDBYTES],
             const polyveck *t1)
{
  unsigned int i;

  for(i = 0; i < SEEDBYTES; ++i)
    pk[i] = rho[i];
  pk += SEEDBYTES;

  for(i = 0; i < K; ++i)
    polyt1_pack(pk + i*POLYT1_PACKEDBYTES, &t1->vec[i]);
}

/*************************************************
* Name:        unpack_pk
*
* Description: Unpack public key pk = (rho, t1).
*
* Arguments:   - const uint8_t rho[]: output byte array for rho
*              - const polyveck *t1: pointer to output vector t1
*              - uint8_t pk[]: byte array containing bit-packed pk
**************************************************/
void unpack_pk(uint8_t rho[SEEDBYTES],
               polyveck *t1,
               const uint8_t pk[CRYPTO_PUBLICKEYBYTES])
{
  unsigned int i;

  for(i = 0; i < SEEDBYTES; ++i)
    rho[i] = pk[i];
  pk += SEEDBYTES;

  for(i = 0; i < K; ++i)
    polyt1_unpack(&t1->vec[i], pk + i*POLYT1_PACKEDBYTES);
}

/*************************************************
* Name:        pack_sk
*
* Description: Bit-pack secret key sk = (rho, tr, key, t0, s1, s2).
*
* Arguments:   - uint8_t sk[]: output byte array
*              - const uint8_t rho[]: byte array containing rho
*              - const uint8_t tr[]: byte array containing tr
*              - const uint8_t key[]: byte array containing key
*              - const polyveck *t0: pointer to vector t0
*              - const polyvecl *s1: pointer to vector s1
*              - const polyveck *s2: pointer to vector s2
**************************************************/
void pack_sk(uint8_t sk[CRYPTO_SECRETKEYBYTES],
             const uint8_t rho[SEEDBYTES],
             const uint8_t tr[TRBYTES],
             const uint8_t key[SEEDBYTES],
             const polyveck *t0,
             const polyvecl *s1,
             const polyveck *s2)
{
  unsigned int i;

  for(i = 0; i < SEEDBYTES; ++i)
    sk[i] = rho[i];
  sk += SEEDBYTES;

  for(i = 0; i < SEEDBYTES; ++i)
    sk[i] = key[i];
  sk += SEEDBYTES;

  for(i = 0; i < TRBYTES; ++i)
    sk[i] = tr[i];
  sk += TRBYTES;

  for(i = 0; i < L; ++i)
    polyeta_pack(sk + i*POLYETA_PACKEDBYTES, &s1->vec[i]);
  sk += L*POLYETA_PACKEDBYTES;

  for(i = 0; i < K; ++i)
    polyeta_pack(sk + i*POLYETA_PACKEDBYTES, &s2->vec[i]);
  sk += K*POLYETA_PACKEDBYTES;

  for(i = 0; i < K; ++i)
    polyt0_pack(sk + i*POLYT0_PACKEDBYTES, &t0->vec[i]);
}

/*************************************************
* Name:        unpack_sk
*
* Description: Unpack secret key sk = (rho, tr, key, t0, s1, s2).
*
* Arguments:   - const uint8_t rho[]: output byte array for rho
*              - const uint8_t tr[]: output byte array for tr
*              - const uint8_t key[]: output byte array for key
*              - const polyveck *t0: pointer to output vector t0
*              - const polyvecl *s1: pointer to output vector s1
*              - const polyveck *s2: pointer to output vector s2
*              - uint8_t sk[]: byte array containing bit-packed sk
**************************************************/
void unpack_sk(uint8_t rho[SEEDBYTES],
               uint8_t tr[TRBYTES],
               uint8_t key[SEEDBYTES],
               polyveck *t0,
               polyvecl *s1,
               polyveck *s2,
               const uint8_t sk[CRYPTO_SECRETKEYBYTES])
{
  unsigned int i;

  for(i = 0; i < SEEDBYTES; ++i)
    rho[i] = sk[i];
  sk += SEEDBYTES;

  for(i = 0; i < SEEDBYTES; ++i)
    key[i] = sk[i];
  sk += SEEDBYTES;

  for(i = 0; i < TRBYTES; ++i)
    tr[i] = sk[i];
  sk += TRBYTES;

  for(i=0; i < L; ++i)
    polyeta_unpack(&s1->vec[i], sk + i*POLYETA_PACKEDBYTES);
  sk += L*POLYETA_PACKEDBYTES;

  for(i=0; i < K; ++i)
    polyeta_unpack(&s2->vec[i], sk + i*POLYETA_PACKEDBYTES);
  sk += K*POLYETA_PACKEDBYTES;

  for(i=0; i < K; ++i)
    polyt0_unpack(&t0->vec[i], sk + i*POLYT0_PACKEDBYTES);
}

/*************************************************
* Name:        pack_sig
*
* Description: Bit-pack signature sig = (c, z, h).
*
* Arguments:   - uint8_t sig[]: output byte array
*              - const uint8_t *c: pointer to challenge hash length SEEDBYTES
*              - const polyvecl *z: pointer to vector z
*              - const polyveck *h: pointer to hint vector h
**************************************************/
void pack_sig(uint8_t sig[CRYPTO_BYTES],
              const uint8_t c[CTILDEBYTES],
              const polyvecl *z,
              const polyveck *h)
{
  unsigned int i, j, k;

  for(i=0; i < CTILDEBYTES; ++i)
    sig[i] = c[i];
  sig += CTILDEBYTES;

  for(i = 0; i < L; ++i)
    polyz_pack(sig + i*POLYZ_PACKEDBYTES, &z->vec[i]);
  sig += L*POLYZ_PACKEDBYTES;

  /* Encode h */
  for(i = 0; i < OMEGA + K; ++i)
    sig[i] = 0;

  k = 0;
  for(i = 0; i < K; ++i) {
    for(j = 0; j < N; ++j)
      if(h->vec[i].coeffs[j] != 0)
        sig[k++] = j;

    sig[OMEGA + i] = k;
  }
}

/*************************************************
* Name:        unpack_sig
*
* Description: Unpack signature sig = (c, z, h).
*
* Arguments:   - uint8_t *c: pointer to output challenge hash
*              - polyvecl *z: pointer to output vector z
*              - polyveck *h: pointer to output hint vector h
*              - const uint8_t sig[]: byte array containing
*                bit-packed signature
*
* Returns 1 in case of malformed signature; otherwise 0.
**************************************************/
int unpack_sig(uint8_t c[CTILDEBYTES],
               polyvecl *z,
               polyveck *h,
               const uint8_t sig[CRYPTO_BYTES])
{
  unsigned int i, j, k;

  for(i = 0; i < CTILDEBYTES; ++i)
    c[i] = sig[i];
  sig += CTILDEBYTES;

  for(i = 0; i < L; ++i)
    polyz_unpack(&z->vec[i], sig + i*POLYZ_PACKEDBYTES);
  sig += L*POLYZ_PACKEDBYTES;

  /* Decode h */
  k = 0;
  for(i = 0; i < K; ++i) {
    for(j = 0; j < N; ++j)
      h->vec[i].coeffs[j] = 0;

    if(sig[OMEGA + i] < k || sig[OMEGA + i] > OMEGA)
      return 1;

    for(j = k; j < sig[OMEGA + i]; ++j) {
      /* Coefficients are ordered for strong unforgeability */
      if(j > k && sig[j] <= sig[j-1]) return 1;
      h->vec[i].coeffs[sig[j]] = 1;
    }

    k = sig[OMEGA + i];
  }

  /* Extra indices are zero for strong unforgeability */
  for(j = k; j < OMEGA; ++j)
    if(sig[j])
      return 1;

  return 0;
}
