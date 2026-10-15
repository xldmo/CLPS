#ifndef PACKING_H
#define PACKING_H

#include <stdint.h>
#include "params.h"
#include "polyvec.h"

#define sclls_pack_params DILITHIUM_NAMESPACE(sclls_pack_params)
void sclls_pack_params(uint8_t params[SCLLS_PARAMSBYTES],
                       const uint8_t rho[SEEDBYTES],
                       const polyveck *t);

#define sclls_unpack_params DILITHIUM_NAMESPACE(sclls_unpack_params)
void sclls_unpack_params(uint8_t rho[SEEDBYTES],
                         polyveck *t,
                         const uint8_t params[SCLLS_PARAMSBYTES]);

#define sclls_pack_msk DILITHIUM_NAMESPACE(sclls_pack_msk)
void sclls_pack_msk(uint8_t msk[SCLLS_MSKBYTES],
                    const uint8_t key[SEEDBYTES],
                    const uint8_t tr[TRBYTES],
                    const polyvecl *s1,
                    const polyveck *s2);

#define sclls_unpack_msk DILITHIUM_NAMESPACE(sclls_unpack_msk)
void sclls_unpack_msk(uint8_t key[SEEDBYTES],
                      uint8_t tr[TRBYTES],
                      polyvecl *s1,
                      polyveck *s2,
                      const uint8_t msk[SCLLS_MSKBYTES]);

#define sclls_pack_did DILITHIUM_NAMESPACE(sclls_pack_did)
void sclls_pack_did(uint8_t did[SCLLS_DIDBYTES],
                    const uint8_t cid[CTILDEBYTES],
                    const polyvecl *zid);

#define sclls_unpack_did DILITHIUM_NAMESPACE(sclls_unpack_did)
void sclls_unpack_did(uint8_t cid[CTILDEBYTES],
                      polyvecl *zid,
                      const uint8_t did[SCLLS_DIDBYTES]);

#define sclls_pack_pid DILITHIUM_NAMESPACE(sclls_pack_pid)
void sclls_pack_pid(uint8_t pid[SCLLS_PIDBYTES],
                    const polyveck *t,
                    const uint8_t w1buf[K*POLYW1_PACKEDBYTES]);

#define sclls_unpack_pid DILITHIUM_NAMESPACE(sclls_unpack_pid)
void sclls_unpack_pid(uint8_t w1buf[K*POLYW1_PACKEDBYTES],
                      polyveck *t,
                      const uint8_t pid[SCLLS_PIDBYTES]);

#define sclls_pack_sid DILITHIUM_NAMESPACE(sclls_pack_sid)
void sclls_pack_sid(uint8_t sid[SCLLS_SIDBYTES],
                    const polyvecl *s1,
                    const polyveck *s2,
                    const uint8_t cid[SCLLS_CIDBYTES],
                    const uint8_t zid[SCLLS_ZIDBYTES]);

#define sclls_unpack_sid DILITHIUM_NAMESPACE(sclls_unpack_sid)
void sclls_unpack_sid(polyvecl *s1,
                      polyveck *s2,
                      uint8_t cid[SCLLS_ZIDBYTES],
                      uint8_t zid[SCLLS_CIDBYTES],
                      const uint8_t sid[SCLLS_SIDBYTES]);

#define sclls_pack_sig DILITHIUM_NAMESPACE(sclls_pack_sig)
void sclls_pack_sig(uint8_t sig[SCLLS_SIGBYTES],
                    const uint8_t cid[CTILDEBYTES],
                    const polyvecl *z0,
                    const polyveck *hw);

#define sclls_unpack_sig DILITHIUM_NAMESPACE(sclls_unpack_sig)
void sclls_unpack_sig(uint8_t cid[CTILDEBYTES],
                      polyvecl *z0,
                      polyveck *hw,
                      const uint8_t sig[SCLLS_SIGBYTES]);

#define pack_pk DILITHIUM_NAMESPACE(pack_pk)
void pack_pk(uint8_t pk[CRYPTO_PUBLICKEYBYTES], const uint8_t rho[SEEDBYTES], const polyveck *t1);

#define pack_sk DILITHIUM_NAMESPACE(pack_sk)
void pack_sk(uint8_t sk[CRYPTO_SECRETKEYBYTES],
             const uint8_t rho[SEEDBYTES],
             const uint8_t tr[TRBYTES],
             const uint8_t key[SEEDBYTES],
             const polyveck *t0,
             const polyvecl *s1,
             const polyveck *s2);

#define pack_sig DILITHIUM_NAMESPACE(pack_sig)
void pack_sig(uint8_t sig[CRYPTO_BYTES], const uint8_t c[CTILDEBYTES], const polyvecl *z, const polyveck *h);

#define unpack_pk DILITHIUM_NAMESPACE(unpack_pk)
void unpack_pk(uint8_t rho[SEEDBYTES], polyveck *t1, const uint8_t pk[CRYPTO_PUBLICKEYBYTES]);

#define unpack_sk DILITHIUM_NAMESPACE(unpack_sk)
void unpack_sk(uint8_t rho[SEEDBYTES],
               uint8_t tr[TRBYTES],
               uint8_t key[SEEDBYTES],
               polyveck *t0,
               polyvecl *s1,
               polyveck *s2,
               const uint8_t sk[CRYPTO_SECRETKEYBYTES]);

#define unpack_sig DILITHIUM_NAMESPACE(unpack_sig)
int unpack_sig(uint8_t c[CTILDEBYTES], polyvecl *z, polyveck *h, const uint8_t sig[CRYPTO_BYTES]);

#endif
