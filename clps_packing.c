/**
 * @file clps_packing.c
 * @brief CLPS 方案的字节打包/解包函数实现
 *
 * 数据布局约定:
 *   params       = rho[SEEDBYTES] || t_KGC[0..K-1] 每个 SCLLS_POLYT_PACKEDBYTES (736)
 *   msk          = rho'[CRHBYTES] || s1_KGC[0..L-1] || s2_KGC[0..K-1] (polyeta_pack)
 *   pk_IDi       = t_IDi[0..K-1] (736) || ~c_IDi[CTILDEBYTES]
 *   partialkey_i = ~c_IDi[CTILDEBYTES] || z_IDi[0..L-1] (polyz_pack)
 *   sk_IDi       = s1[L] || s2[K] (polyeta_pack) || usk[SEEDBYTES] || ~c[CTILDEBYTES]
 *                  || z[L] (polyz_pack) || tr[TRBYTES]
 *   proxykey_PA  = ID_P[SEEDBYTES] || ID_A[SEEDBYTES] || z_P[L] (polyz_pack)
 *                  || ~c_PA[CTILDEBYTES] || z_PA[L] (polyz_pack)
 *   sigma_PA     = ~c_A[CTILDEBYTES] || ~c_PA[CTILDEBYTES] || z_PA[L] (polyz_pack)
 *                  || z_0_A[L] (CLPS_POLYZ0A_PACKEDBYTES, 22-bit) || h_w[K] (sclls_polyhw_pack)
 */

#include <stdint.h>
#include "params.h"
#include "clps_packing.h"
#include "polyvec.h"
#include "poly.h"

/*===============================================================
 * params: (rho, t_KGC)
 *===============================================================*/

/**
 * @brief 打包 params = (rho || t_KGC)
 *
 * @param params  输出 buffer (CLPS_PARAMSBYTES 字节)
 * @param rho     输入 rho 种子
 * @param t       输入 KGC 公钥向量 t_KGC
 */
void clps_pack_params(uint8_t params[CLPS_PARAMSBYTES],
                      const uint8_t rho[SEEDBYTES],
                      const polyveck* t)
{
  unsigned int i;

  for (i = 0; i < SEEDBYTES; ++i)
    params[i] = rho[i];
  params += SEEDBYTES;

  for (i = 0; i < K; ++i)
    sclls_polyt_pack(params + i * SCLLS_POLYT_PACKEDBYTES, &t->vec[i]);
}

/**
 * @brief 解包 params 得到 (rho, t_KGC)
 *
 * @param rho    输出 rho 种子
 * @param t      输出 KGC 公钥向量 t_KGC
 * @param params 输入 params buffer
 */
void clps_unpack_params(uint8_t rho[SEEDBYTES],
                        polyveck* t,
                        const uint8_t params[CLPS_PARAMSBYTES])
{
  unsigned int i;

  for (i = 0; i < SEEDBYTES; ++i)
    rho[i] = params[i];
  params += SEEDBYTES;

  for (i = 0; i < K; ++i)
    sclls_polyt_unpack(&t->vec[i], params + i * SCLLS_POLYT_PACKEDBYTES);
}

/*===============================================================
 * msk: (rho', s1_KGC, s2_KGC)
 *===============================================================*/

/**
 * @brief 打包 msk = (rho' || s1_KGC || s2_KGC)
 *
 * @param msk      输出 buffer (CLPS_MSKBYTES 字节)
 * @param rhoprime 输入 rho' 种子 (CRHBYTES)
 * @param s1       输入 KGC 秘密向量 s1_KGC
 * @param s2       输入 KGC 秘密向量 s2_KGC
 */
void clps_pack_msk(uint8_t msk[CLPS_MSKBYTES],
                   const uint8_t rhoprime[CRHBYTES],
                   const polyvecl* s1,
                   const polyveck* s2)
{
  unsigned int i;

  for (i = 0; i < CRHBYTES; ++i)
    msk[i] = rhoprime[i];
  msk += CRHBYTES;

  for (i = 0; i < L; ++i)
    polyeta_pack(msk + i * POLYETA_PACKEDBYTES, &s1->vec[i]);
  msk += L * POLYETA_PACKEDBYTES;

  for (i = 0; i < K; ++i)
    polyeta_pack(msk + i * POLYETA_PACKEDBYTES, &s2->vec[i]);
}

/**
 * @brief 解包 msk 得到 (rho', s1_KGC, s2_KGC)
 *
 * @param rhoprime 输出 rho' 种子
 * @param s1       输出 s1_KGC
 * @param s2       输出 s2_KGC
 * @param msk      输入 msk buffer
 */
void clps_unpack_msk(uint8_t rhoprime[CRHBYTES],
                     polyvecl* s1,
                     polyveck* s2,
                     const uint8_t msk[CLPS_MSKBYTES])
{
  unsigned int i;
  const uint8_t* p = msk;

  for (i = 0; i < CRHBYTES; ++i)
    rhoprime[i] = p[i];
  p += CRHBYTES;

  for (i = 0; i < L; ++i)
    polyeta_unpack(&s1->vec[i], p + i * POLYETA_PACKEDBYTES);
  p += L * POLYETA_PACKEDBYTES;

  for (i = 0; i < K; ++i)
    polyeta_unpack(&s2->vec[i], p + i * POLYETA_PACKEDBYTES);
}

/*===============================================================
 * pk_IDi: (t, ~c_IDi)
 *===============================================================*/

/**
 * @brief 打包 pk_IDi = (t_IDi || ~c_IDi)
 *
 * @param pk  输出 buffer (CLPS_PKIDBYTES 字节)
 * @param t   输入用户公钥向量 t_IDi
 * @param cid 输入挑战哈希 ~c_IDi
 */
void clps_pack_pk(uint8_t pk[CLPS_PKIDBYTES],
                  const polyveck* t,
                  const uint8_t cid[CTILDEBYTES])
{
  unsigned int i;

  for (i = 0; i < K; ++i)
    sclls_polyt_pack(pk + i * SCLLS_POLYT_PACKEDBYTES, &t->vec[i]);
  pk += K * SCLLS_POLYT_PACKEDBYTES;

  for (i = 0; i < CTILDEBYTES; ++i)
    pk[i] = cid[i];
}

/**
 * @brief 解包 pk_IDi 得到 (t, ~c_IDi)
 *
 * @param t   输出用户公钥向量 t_IDi
 * @param cid 输出挑战哈希 ~c_IDi
 * @param pk  输入 pk buffer
 */
void clps_unpack_pk(polyveck* t,
                    uint8_t cid[CTILDEBYTES],
                    const uint8_t pk[CLPS_PKIDBYTES])
{
  unsigned int i;

  for (i = 0; i < K; ++i)
    sclls_polyt_unpack(&t->vec[i], pk + i * SCLLS_POLYT_PACKEDBYTES);
  pk += K * SCLLS_POLYT_PACKEDBYTES;

  for (i = 0; i < CTILDEBYTES; ++i)
    cid[i] = pk[i];
}

/*===============================================================
 * partialkey_i: (~c_IDi, z_IDi)
 *===============================================================*/

/**
 * @brief 打包 partialkey_i = (~c_IDi || z_IDi)
 *
 * @param pk  输出 buffer (CLPS_PARTIALKEYBYTES 字节)
 * @param cid 输入挑战哈希 ~c_IDi
 * @param z   输入 z_IDi 向量 (系数范围 [-(GAMMA1-1), GAMMA1])
 */
void clps_pack_partialkey(uint8_t pk[CLPS_PARTIALKEYBYTES],
                          const uint8_t cid[CTILDEBYTES],
                          const polyvecl* z)
{
  unsigned int i;

  for (i = 0; i < CTILDEBYTES; ++i)
    pk[i] = cid[i];
  pk += CTILDEBYTES;

  for (i = 0; i < L; ++i)
    polyz_pack(pk + i * POLYZ_PACKEDBYTES, &z->vec[i]);
}

/**
 * @brief 解包 partialkey_i 得到 (~c_IDi, z_IDi)
 *
 * @param cid 输出挑战哈希 ~c_IDi
 * @param z   输出 z_IDi 向量
 * @param pk  输入 partialkey buffer
 */
void clps_unpack_partialkey(uint8_t cid[CTILDEBYTES],
                            polyvecl* z,
                            const uint8_t pk[CLPS_PARTIALKEYBYTES])
{
  unsigned int i;
  const uint8_t* p = pk;

  for (i = 0; i < CTILDEBYTES; ++i)
    cid[i] = p[i];
  p += CTILDEBYTES;

  for (i = 0; i < L; ++i)
    polyz_unpack(&z->vec[i], p + i * POLYZ_PACKEDBYTES);
}

/*===============================================================
 * sk_IDi: (s1, s2, usk, ~c_IDi, z_IDi, tr)
 *===============================================================*/

/**
 * @brief 打包 sk_IDi
 *
 * 布局: s1[L] || s2[K] (polyeta_pack) || usk[SEEDBYTES] || ~c[CTILDEBYTES]
 *       || z[L] (polyz_pack) || tr[TRBYTES]
 *
 * @param sk  输出 buffer (CLPS_SKIDBYTES 字节)
 * @param s1  输入 s1_IDi (系数范围 [-ETA, ETA])
 * @param s2  输入 s2_IDi (系数范围 [-ETA, ETA])
 * @param usk 输入用户秘密种子 usk_i
 * @param cid 输入挑战哈希 ~c_IDi
 * @param z   输入 z_IDi 向量
 * @param tr  输入追踪值 tr_IDi
 */
void clps_pack_sk(uint8_t sk[CLPS_SKIDBYTES],
                  const polyvecl* s1,
                  const polyveck* s2,
                  const uint8_t usk[SEEDBYTES],
                  const uint8_t cid[CTILDEBYTES],
                  const polyvecl* z,
                  const uint8_t tr[TRBYTES])
{
  unsigned int i;

  for (i = 0; i < L; ++i)
    polyeta_pack(sk + i * POLYETA_PACKEDBYTES, &s1->vec[i]);
  sk += L * POLYETA_PACKEDBYTES;

  for (i = 0; i < K; ++i)
    polyeta_pack(sk + i * POLYETA_PACKEDBYTES, &s2->vec[i]);
  sk += K * POLYETA_PACKEDBYTES;

  for (i = 0; i < SEEDBYTES; ++i)
    sk[i] = usk[i];
  sk += SEEDBYTES;

  for (i = 0; i < CTILDEBYTES; ++i)
    sk[i] = cid[i];
  sk += CTILDEBYTES;

  for (i = 0; i < L; ++i)
    polyz_pack(sk + i * POLYZ_PACKEDBYTES, &z->vec[i]);
  sk += L * POLYZ_PACKEDBYTES;

  for (i = 0; i < TRBYTES; ++i)
    sk[i] = tr[i];
}

/**
 * @brief 解包 sk_IDi 得到 (s1, s2, usk, ~c, z, tr)
 *
 * @param s1  输出 s1_IDi
 * @param s2  输出 s2_IDi
 * @param usk 输出 usk_i
 * @param cid 输出 ~c_IDi
 * @param z   输出 z_IDi
 * @param tr  输出 tr_IDi
 * @param sk  输入 sk buffer
 */
void clps_unpack_sk(polyvecl* s1,
                    polyveck* s2,
                    uint8_t usk[SEEDBYTES],
                    uint8_t cid[CTILDEBYTES],
                    polyvecl* z,
                    uint8_t tr[TRBYTES],
                    const uint8_t sk[CLPS_SKIDBYTES])
{
  unsigned int i;
  const uint8_t* p = sk;

  for (i = 0; i < L; ++i)
    polyeta_unpack(&s1->vec[i], p + i * POLYETA_PACKEDBYTES);
  p += L * POLYETA_PACKEDBYTES;

  for (i = 0; i < K; ++i)
    polyeta_unpack(&s2->vec[i], p + i * POLYETA_PACKEDBYTES);
  p += K * POLYETA_PACKEDBYTES;

  for (i = 0; i < SEEDBYTES; ++i)
    usk[i] = p[i];
  p += SEEDBYTES;

  for (i = 0; i < CTILDEBYTES; ++i)
    cid[i] = p[i];
  p += CTILDEBYTES;

  for (i = 0; i < L; ++i)
    polyz_unpack(&z->vec[i], p + i * POLYZ_PACKEDBYTES);
  p += L * POLYZ_PACKEDBYTES;

  for (i = 0; i < TRBYTES; ++i)
    tr[i] = p[i];
}

/*===============================================================
 * proxykey_PA: (ID_P, ID_A, z_P, ~c_PA, z_PA)
 *===============================================================*/

/**
 * @brief 打包代理密钥 proxykey_PA
 *
 * 布局: ID_P[SEEDBYTES] || ID_A[SEEDBYTES] || z_P[L] (polyz_pack)
 *       || ~c_PA[CTILDEBYTES] || z_PA[L] (polyz_pack)
 *
 * @param t      输出 buffer (CLPS_PROXYKEYBYTES 字节)
 * @param idp    输入授权者身份 ID_P
 * @param ida    输入代理者身份 ID_A
 * @param zp     输入 z_P 向量 (来自主用户 PartialKeyGen)
 * @param cidpa  输入挑战哈希 ~c_PA
 * @param zpa    输入 z_PA 向量
 */
void clps_pack_proxykey(uint8_t t[CLPS_PROXYKEYBYTES],
                        const uint8_t idp[SEEDBYTES],
                        const uint8_t ida[SEEDBYTES],
                        const polyvecl* zp,
                        const uint8_t cidpa[CTILDEBYTES],
                        const polyvecl* zpa)
{
  unsigned int i;

  for (i = 0; i < SEEDBYTES; ++i)
    t[i] = idp[i];
  t += SEEDBYTES;

  for (i = 0; i < SEEDBYTES; ++i)
    t[i] = ida[i];
  t += SEEDBYTES;

  for (i = 0; i < L; ++i)
    polyz_pack(t + i * POLYZ_PACKEDBYTES, &zp->vec[i]);
  t += L * POLYZ_PACKEDBYTES;

  for (i = 0; i < CTILDEBYTES; ++i)
    t[i] = cidpa[i];
  t += CTILDEBYTES;

  for (i = 0; i < L; ++i)
    polyz_pack(t + i * POLYZ_PACKEDBYTES, &zpa->vec[i]);
}

/**
 * @brief 解包代理密钥 proxykey_PA
 *
 * @param idp    输出 ID_P
 * @param ida    输出 ID_A
 * @param zp     输出 z_P
 * @param cidpa  输出 ~c_PA
 * @param zpa    输出 z_PA
 * @param t      输入 proxykey buffer
 */
void clps_unpack_proxykey(uint8_t idp[SEEDBYTES],
                          uint8_t ida[SEEDBYTES],
                          polyvecl* zp,
                          uint8_t cidpa[CTILDEBYTES],
                          polyvecl* zpa,
                          const uint8_t t[CLPS_PROXYKEYBYTES])
{
  unsigned int i;
  const uint8_t* p = t;

  for (i = 0; i < SEEDBYTES; ++i)
    idp[i] = p[i];
  p += SEEDBYTES;

  for (i = 0; i < SEEDBYTES; ++i)
    ida[i] = p[i];
  p += SEEDBYTES;

  for (i = 0; i < L; ++i)
    polyz_unpack(&zp->vec[i], p + i * POLYZ_PACKEDBYTES);
  p += L * POLYZ_PACKEDBYTES;

  for (i = 0; i < CTILDEBYTES; ++i)
    cidpa[i] = p[i];
  p += CTILDEBYTES;

  for (i = 0; i < L; ++i)
    polyz_unpack(&zpa->vec[i], p + i * POLYZ_PACKEDBYTES);
}

/*===============================================================
 * sigma_PA: (~c_A, ~c_PA, z_PA, z_0_A, h_w)
 *===============================================================*/

/**
 * @brief 打包代理签名 sigma_PA
 *
 * 布局: ~c_A[CTILDEBYTES] || ~c_PA[CTILDEBYTES] || z_PA[L] (polyz_pack)
 *       || z_0_A[L] (clps_polyz0a_pack, 22-bit) || h_w[K] (sclls_polyhw_pack)
 *
 * @param sig   输出 buffer (CLPS_SIGBYTES 字节)
 * @param cida  输入 ~c_A
 * @param cidpa 输入 ~c_PA
 * @param zpa   输入 z_PA 向量
 * @param z0a   输入 z_0_A 向量 (系数范围 (-3*(GAMMA1-BETA), 3*(GAMMA1-BETA)))
 * @param hw    输入 h_w 向量 (提示)
 */
void clps_pack_sigma(uint8_t sig[CLPS_SIGBYTES],
                     const uint8_t cida[CTILDEBYTES],
                     const uint8_t cidpa[CTILDEBYTES],
                     const polyvecl* zpa,
                     const polyvecl* z0a,
                     const polyveck* hw)
{
  unsigned int i;

  for (i = 0; i < CTILDEBYTES; ++i)
    sig[i] = cida[i];
  sig += CTILDEBYTES;

  for (i = 0; i < CTILDEBYTES; ++i)
    sig[i] = cidpa[i];
  sig += CTILDEBYTES;

  for (i = 0; i < L; ++i)
    polyz_pack(sig + i * POLYZ_PACKEDBYTES, &zpa->vec[i]);
  sig += L * POLYZ_PACKEDBYTES;

  for (i = 0; i < L; ++i)
    clps_polyz0a_pack(sig + i * CLPS_POLYZ0A_PACKEDBYTES, &z0a->vec[i]);
  sig += L * CLPS_POLYZ0A_PACKEDBYTES;

  for (i = 0; i < K; ++i)
    sclls_polyhw_pack(sig + i * SCLLS_HWBYTES_PERPOLY, &hw->vec[i]);
}

/**
 * @brief 解包代理签名 sigma_PA
 *
 * @param cida  输出 ~c_A
 * @param cidpa 输出 ~c_PA
 * @param zpa   输出 z_PA 向量
 * @param z0a   输出 z_0_A 向量
 * @param hw    输出 h_w 向量
 * @param sig   输入 sigma buffer
 */
void clps_unpack_sigma(uint8_t cida[CTILDEBYTES],
                       uint8_t cidpa[CTILDEBYTES],
                       polyvecl* zpa,
                       polyvecl* z0a,
                       polyveck* hw,
                       const uint8_t sig[CLPS_SIGBYTES])
{
  unsigned int i;
  const uint8_t* p = sig;

  for (i = 0; i < CTILDEBYTES; ++i)
    cida[i] = p[i];
  p += CTILDEBYTES;

  for (i = 0; i < CTILDEBYTES; ++i)
    cidpa[i] = p[i];
  p += CTILDEBYTES;

  for (i = 0; i < L; ++i)
    polyz_unpack(&zpa->vec[i], p + i * POLYZ_PACKEDBYTES);
  p += L * POLYZ_PACKEDBYTES;

  for (i = 0; i < L; ++i)
    clps_polyz0a_unpack(&z0a->vec[i], p + i * CLPS_POLYZ0A_PACKEDBYTES);
  p += L * CLPS_POLYZ0A_PACKEDBYTES;

  for (i = 0; i < K; ++i)
    sclls_polyhw_unpack(&hw->vec[i], p + i * SCLLS_HWBYTES_PERPOLY);
}
