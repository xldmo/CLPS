/**
 * @file clps_packing.h
 * @brief CLPS 方案的字节打包/解包函数声明
 *
 * CLPS 中各数据结构的字节布局:
 *   params      = rho[SEEDBYTES] || t_KGC[0..K-1] (每个 SCLLS_POLYT_PACKEDBYTES 字节)
 *   msk         = rho'[CRHBYTES] || s1_KGC[0..L-1] || s2_KGC[0..K-1]
 *   pk_IDi      = t_IDi[0..K-1] (SCLLS_POLYT_PACKEDBYTES) || ~c_IDi[CTILDEBYTES]
 *   partialkey_i= ~c_IDi[CTILDEBYTES] || z_IDi[0..L-1] (polyz_pack)
 *   sk_IDi      = s1_IDi[0..L-1] || s2_IDi[0..K-1] || usk_i[SEEDBYTES]
 *                 || ~c_IDi[CTILDEBYTES] || z_IDi[0..L-1] || tr_IDi[TRBYTES]
 *   proxykey_PA = ID_P[SEEDBYTES] || ID_A[SEEDBYTES] || z_P[0..L-1]
 *                 || ~c_PA[CTILDEBYTES] || z_PA[0..L-1]
 *   sigma_PA    = ~c_A[CTILDEBYTES] || ~c_PA[CTILDEBYTES] || z_PA[0..L-1]
 *                 || z_0_A[0..L-1] (CLPS_POLYZ0A_PACKEDBYTES, 22-bit) || h_w[0..K-1]
 */

#ifndef CLPS_PACKING_H
#define CLPS_PACKING_H

#include <stdint.h>
#include "params.h"
#include "polyvec.h"
#include "clps.h"

/*===============================================================
 * params: (rho, t_KGC)
 *===============================================================*/
#define clps_pack_params CLPS_NAMESPACE(clps_pack_params)
void clps_pack_params(uint8_t params[CLPS_PARAMSBYTES],
                      const uint8_t rho[SEEDBYTES],
                      const polyveck* t);

#define clps_unpack_params CLPS_NAMESPACE(clps_unpack_params)
void clps_unpack_params(uint8_t rho[SEEDBYTES],
                        polyveck* t,
                        const uint8_t params[CLPS_PARAMSBYTES]);

/*===============================================================
 * msk: (rho', s1_KGC, s2_KGC)
 *===============================================================*/
#define clps_pack_msk CLPS_NAMESPACE(clps_pack_msk)
void clps_pack_msk(uint8_t msk[CLPS_MSKBYTES],
                   const uint8_t rhoprime[CRHBYTES],
                   const polyvecl* s1,
                   const polyveck* s2);

#define clps_unpack_msk CLPS_NAMESPACE(clps_unpack_msk)
void clps_unpack_msk(uint8_t rhoprime[CRHBYTES],
                     polyvecl* s1,
                     polyveck* s2,
                     const uint8_t msk[CLPS_MSKBYTES]);

/*===============================================================
 * pk_IDi: (t, ~c_IDi)
 *===============================================================*/
#define clps_pack_pk CLPS_NAMESPACE(clps_pack_pk)
void clps_pack_pk(uint8_t pk[CLPS_PKIDBYTES],
                  const polyveck* t,
                  const uint8_t cid[CTILDEBYTES]);

#define clps_unpack_pk CLPS_NAMESPACE(clps_unpack_pk)
void clps_unpack_pk(polyveck* t,
                    uint8_t cid[CTILDEBYTES],
                    const uint8_t pk[CLPS_PKIDBYTES]);

/*===============================================================
 * partialkey_i: (~c_IDi, z_IDi)
 *===============================================================*/
#define clps_pack_partialkey CLPS_NAMESPACE(clps_pack_partialkey)
void clps_pack_partialkey(uint8_t pk[CLPS_PARTIALKEYBYTES],
                          const uint8_t cid[CTILDEBYTES],
                          const polyvecl* z);

#define clps_unpack_partialkey CLPS_NAMESPACE(clps_unpack_partialkey)
void clps_unpack_partialkey(uint8_t cid[CTILDEBYTES],
                            polyvecl* z,
                            const uint8_t pk[CLPS_PARTIALKEYBYTES]);

/*===============================================================
 * sk_IDi: (s1, s2, usk, ~c_IDi, z_IDi, tr)
 *===============================================================*/
#define clps_pack_sk CLPS_NAMESPACE(clps_pack_sk)
void clps_pack_sk(uint8_t sk[CLPS_SKIDBYTES],
                  const polyvecl* s1,
                  const polyveck* s2,
                  const uint8_t usk[SEEDBYTES],
                  const uint8_t cid[CTILDEBYTES],
                  const polyvecl* z,
                  const uint8_t tr[TRBYTES]);

#define clps_unpack_sk CLPS_NAMESPACE(clps_unpack_sk)
void clps_unpack_sk(polyvecl* s1,
                    polyveck* s2,
                    uint8_t usk[SEEDBYTES],
                    uint8_t cid[CTILDEBYTES],
                    polyvecl* z,
                    uint8_t tr[TRBYTES],
                    const uint8_t sk[CLPS_SKIDBYTES]);

/*===============================================================
 * proxykey_PA: (ID_P, ID_A, z_P, ~c_PA, z_PA)
 *===============================================================*/
#define clps_pack_proxykey CLPS_NAMESPACE(clps_pack_proxykey)
void clps_pack_proxykey(uint8_t t[CLPS_PROXYKEYBYTES],
                        const uint8_t idp[SEEDBYTES],
                        const uint8_t ida[SEEDBYTES],
                        const polyvecl* zp,
                        const uint8_t cidpa[CTILDEBYTES],
                        const polyvecl* zpa);

#define clps_unpack_proxykey CLPS_NAMESPACE(clps_unpack_proxykey)
void clps_unpack_proxykey(uint8_t idp[SEEDBYTES],
                          uint8_t ida[SEEDBYTES],
                          polyvecl* zp,
                          uint8_t cidpa[CTILDEBYTES],
                          polyvecl* zpa,
                          const uint8_t t[CLPS_PROXYKEYBYTES]);

/*===============================================================
 * sigma_PA: (~c_A, ~c_PA, z_PA, z_0_A, h_w)
 *===============================================================*/
#define clps_pack_sigma CLPS_NAMESPACE(clps_pack_sigma)
void clps_pack_sigma(uint8_t sig[CLPS_SIGBYTES],
                     const uint8_t cida[CTILDEBYTES],
                     const uint8_t cidpa[CTILDEBYTES],
                     const polyvecl* zpa,
                     const polyvecl* z0a,
                     const polyveck* hw);

#define clps_unpack_sigma CLPS_NAMESPACE(clps_unpack_sigma)
void clps_unpack_sigma(uint8_t cida[CTILDEBYTES],
                       uint8_t cidpa[CTILDEBYTES],
                       polyvecl* zpa,
                       polyvecl* z0a,
                       polyveck* hw,
                       const uint8_t sig[CLPS_SIGBYTES]);

#endif /* CLPS_PACKING_H */
