/**
 * @file clps.c
 * @brief CLPS (Certificateless Proxy Signature) 方案核心实现
 *
 * 基于 ML-DSA (Dilithium) 格密码框架，实现无证书代理签名机制。
 * 包含以下算法:
 *   1. Setup                     - 系统参数生成 (params, msk)
 *   2. PartialKeyGen             - 部分私钥生成 (KGC 签发)
 *   3. FullKeyGen                - 完整用户密钥对生成
 *   4. Principal-Authorize-Agent - 授权代理凭证生成 (proxykey)
 *   5. Agent-Sign                - 代理签名
 *   6. Agent-Verify              - 代理签名验证
 *
 * 参考文献:
 *   - CLPS_scheme_终版.tex (CLPS 伪代码, 唯一实现规格)
 *   - Dilithium / ML-DSA
 *
 * 域状态约定:
 *   - 标准域: 多项式系数在 [-Q, Q) 或 [0, Q) 范围内
 *   - NTT 域: 系数经过 ntt() 转换，用于加速乘法
 *   - 约减操作:
 *     poly_reduce / polyveck_reduce: 模 Q 约减
 *     poly_caddq / polyveck_caddq:   规范化到 [0, Q) 标准代表元
 *   - 矩阵乘法流程: ntt -> pointwise_montgomery -> reduce -> invntt_tomont -> caddq
 *
 * 关键约定 (与 tex 算法正文一致):
 *   - z_0_A 范数界: ||z_0_A||_inf < 3*(GAMMA1-BETA) = CLPS_Z0A_BOUND
 *   - z_0_A 打包: 每系数 22 bit (CLPS_POLYZ0A_PACKEDBYTES = 704 字节/多项式)
 *   - h_w: Makeh_w 三输入版, 值集 {-1,0,1,-17,-16,-15} 复用 sclls_polyhw_pack
 */

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "params.h"
#include "clps.h"
#include "clps_packing.h"
#include "polyvec.h"
#include "poly.h"
#include "symmetric.h"
#include "fips202.h"

/* ---- 基准测试插桩: 全局拒绝采样计数器 (bench.c 使用) ----
 * clps_rej_partial_keygen : clps_partial_key_gen 的拒绝采样循环迭代次数
 * clps_rej_authorize      : clps_principal_authorize_agent 的拒绝采样循环迭代次数
 * clps_rej_agent_sign     : clps_agent_sign 的拒绝采样循环迭代次数
 */
unsigned long long clps_rej_partial_keygen = 0;
unsigned long long clps_rej_authorize      = 0;
unsigned long long clps_rej_agent_sign     = 0;

/*===============================================================
 * 内部辅助: tr = H1(rho || t)
 *===============================================================*/
static void clps_h1_tr(uint8_t tr[TRBYTES],
                       const uint8_t rho[SEEDBYTES],
                       const polyveck *t)
{
  uint8_t pbuf[CLPS_PARAMSBYTES];
  polyveck tq;
  keccak_state state;

  tq = *t;
  polyveck_caddq(&tq);
  clps_pack_params(pbuf, rho, &tq);
  shake256_init(&state);
  shake256_absorb(&state, pbuf, CLPS_PARAMSBYTES);
  shake256_finalize(&state);
  shake256_squeeze(tr, TRBYTES, &state);
}

/*===============================================================
 * 内部辅助: tr = H1(params) 其中 params = (rho || t_KGC 打包)
 *===============================================================*/
static void clps_h1_tr_params(uint8_t tr[TRBYTES],
                              const uint8_t params[CLPS_PARAMSBYTES])
{
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, params, CLPS_PARAMSBYTES);
  shake256_finalize(&state);
  shake256_squeeze(tr, TRBYTES, &state);
}


/*===============================================================
 * clps_setup —— Setup: 系统参数生成
 *
 * 伪代码对应 (Algorithm Setup):
 *   rho  <- {0,1}^256
 *   rho' <- {0,1}^512
 *   A    := ExpandA(rho)
 *   (s1_KGC, s2_KGC) := ExpandS(rho')
 *   t_KGC := A*s1_KGC + s2_KGC
 *   tr    := H1(rho || t_KGC)
 *   params := (rho, t_KGC)
 *   msk    := (rho', s1_KGC, s2_KGC)
 *
 * 代码作用:
 *   从随机种子 seed 经单次 shake256 输出 96 字节派生:
 *     databuf[0..31]  = rho  (矩阵 A 的生成种子)
 *     databuf[32..95] = rho' (KGC 私钥的派生种子)
 *
 * 输出:
 *   params = rho || t_KGC (CLPS_PARAMSBYTES 字节)
 *   msk    = rho' || s1_KGC || s2_KGC (CLPS_MSKBYTES 字节)
 *===============================================================*/
int clps_setup(uint8_t* params, uint8_t* msk, const uint8_t* seed)
{
  uint8_t databuf[SEEDBYTES + CRHBYTES];  /* rho(32) || rho'(64) */
  uint8_t seed_in[SEEDBYTES + 2];
  polyvecl mat[K], s1, s1hat;
  polyveck s2, t;

  /* 从 seed 派生 rho 与 rho' (域名标签 K, L 确保独立) */
  memcpy(seed_in, seed, SEEDBYTES);
  seed_in[SEEDBYTES + 0] = K;
  seed_in[SEEDBYTES + 1] = L;
  shake256(databuf, SEEDBYTES + CRHBYTES, seed_in, SEEDBYTES + 2);

  /* A = ExpandA(rho) */
  polyvec_matrix_expand(mat, databuf);

  /* (s1_KGC, s2_KGC) = ExpandS(rho') */
  polyvecl_uniform_eta(&s1, databuf + SEEDBYTES, 0);
  polyveck_uniform_eta(&s2, databuf + SEEDBYTES, L);

  /* t_KGC = A*s1_KGC + s2_KGC */
  s1hat = s1;
  polyvecl_ntt(&s1hat);
  polyvec_matrix_pointwise_montgomery(&t, mat, &s1hat);
  polyveck_reduce(&t);
  polyveck_invntt_tomont(&t);
  polyveck_add(&t, &t, &s2);
  polyveck_caddq(&t);

  /* params := (rho || t_KGC), msk := (rho' || s1_KGC || s2_KGC) */
  clps_pack_params(params, databuf, &t);
  clps_pack_msk(msk, databuf + SEEDBYTES, &s1, &s2);

  return 0;
}


/*===============================================================
 * clps_partial_key_gen —— PartialKeyGen: 部分私钥生成 (KGC 签发)
 *
 * 伪代码对应 (Algorithm PartialKeyGen):
 *   mu_i = H2(tr || ID_i)
 *   rho_i = H2(rho' || mu_i)
 *   循环:
 *     y_i = ExpandMask(rho_i, kappa)
 *     w_i = A*y_i
 *     w1_i = HighBits(w_i, 2*gamma2)
 *     ~c_i = H3(mu_i || Packw1(w1_i))
 *     c_i = SampleInBall(~c_i)
 *     z_i = y_i + c_i*s1_KGC
 *     r0_i = LowBits(w_i - c_i*s2_KGC, 2*gamma2)
 *     (拒绝采样: ||z_i||_inf < gamma1-beta 且 ||r0_i||_inf < gamma2-beta)
 *   partialkey_i := (~c_i, z_i)
 *
 * 注: tex 中 z_i := y_P + c_i*s1_KGC 的 "y_P" 为笔误, 应为 y_i。
 *===============================================================*/
int clps_partial_key_gen(uint8_t* partialkey,
                         const uint8_t* params,
                         const uint8_t* msk,
                         const uint8_t* id)
{
  uint8_t tr[TRBYTES];
  uint8_t rho[SEEDBYTES];
  uint8_t rhoprime[CRHBYTES];
  uint8_t mu[CRHBYTES], rhoi[CRHBYTES];
  uint8_t w1buf[K * POLYW1_PACKEDBYTES];
  uint8_t cid[CTILDEBYTES];
  uint16_t nonce = 0;
  unsigned int rej_count = 0;
  keccak_state state;
  polyvecl mat[K], s1, y, z;
  polyveck s2, t, w, w1, w0;
  poly cp;
  polyvecl yhat;

  /* 解包 params -> rho (t_KGC 在本算法中不使用), 解包 msk -> rho', s1_KGC, s2_KGC */
  clps_unpack_params(rho, &t, params);
  clps_unpack_msk(rhoprime, &s1, &s2, msk);

  /* tr = H1(rho || t_KGC), 直接对打包后的 params 字节计算 */
  clps_h1_tr_params(tr, params);

  /* mu_i = H2(tr || ID_i) */
  shake256_init(&state);
  shake256_absorb(&state, tr, TRBYTES);
  shake256_absorb(&state, id, SEEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

  /* rho_i = H2(rho' || mu_i) */
  shake256_init(&state);
  shake256_absorb(&state, rhoprime, CRHBYTES);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_finalize(&state);
  shake256_squeeze(rhoi, CRHBYTES, &state);

  polyvec_matrix_expand(mat, rho);
  polyvecl_ntt(&s1);
  polyveck_ntt(&s2);

rej:
  rej_count++;
  /* 采样遮掩向量 y_i */
  polyvecl_uniform_gamma1(&y, rhoi, nonce);
  nonce = nonce + L;

  /* w_i = A*y_i */
  yhat = y;
  polyvecl_ntt(&yhat);
  polyvec_matrix_pointwise_montgomery(&w, mat, &yhat);
  polyveck_reduce(&w);
  polyveck_invntt_tomont(&w);

  /* 分解 w_i -> (w1_i, w0_i) */
  polyveck_caddq(&w);
  polyveck_decompose(&w1, &w0, &w);
  polyveck_pack_w1(w1buf, &w1);

  /* ~c_i = H3(mu_i || w_packi) */
  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, w1buf, K * POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(cid, CTILDEBYTES, &state);

  /* c_i = SampleInBall(~c_i) */
  poly_challenge(&cp, cid);
  poly_ntt(&cp);

  /* z_i = y_i + c_i*s1_KGC */
  {
    polyvecl cs1;
    polyvecl_pointwise_poly_montgomery(&cs1, &cp, &s1);
    polyvecl_invntt_tomont(&cs1);
    polyvecl_add(&z, &y, &cs1);
    polyvecl_reduce(&z);
  }

  /* 拒绝采样检查 #1: ||z_i||_inf < GAMMA1 - BETA */
  if (polyvecl_chknorm(&z, GAMMA1 - BETA))
    goto rej;

  /* r0_i = LowBits(w_i - c_i*s2_KGC) */
  {
    polyveck cs2;
    polyveck_pointwise_poly_montgomery(&cs2, &cp, &s2);
    polyveck_reduce(&cs2);
    polyveck_invntt_tomont(&cs2);
    polyveck_sub(&w0, &w0, &cs2);
    polyveck_reduce(&w0);
  }

  /* 拒绝采样检查 #2: ||r0_i||_inf < GAMMA2 - BETA */
  if (polyveck_chknorm(&w0, GAMMA2 - BETA))
    goto rej;

  /* 打包 partialkey_i = (~c_i || z_i) */
  clps_pack_partialkey(partialkey, cid, &z);

#ifndef CLPS_NO_REJ_PRINTF
  printf("[PartialKeyGen] sampling iterations (rej loop): %u\n", rej_count);
#endif
  clps_rej_partial_keygen += rej_count;

  return 0;
}


/*===============================================================
 * clps_full_key_gen —— FullKeyGen: 完整用户密钥对生成
 *
 * 伪代码对应 (Algorithm FullKeyGen):
 *   c_i = SampleInBall(~c_i)
 *   w'1_i = HighBits(A*z_i - c_i*t_KGC, 2*gamma2)
 *   mu_i = H2(tr || ID_i)
 *   ~c'_i = H3(mu_i || Packw1(w'1_i))
 *   若 ||z_i||_inf < gamma1-beta 且 ~c'_i = ~c_i:
 *     xi_i <- {0,1}^256
 *     theta_i = H2(xi_i || usk_i)
 *     (s1_i, s2_i) = ExpandS(theta_i)
 *     t_i = A*s1_i + s2_i
 *     tr_i = H1(rho || t_i)
 *     sk_i := (s1_i, s2_i, usk_i, ~c_i, z_i, tr_i)
 *     pk_i := (t_i, ~c_i)
 *   否则返回失败
 *
 * 注: xi_i 从输入 seed 经带域名标签的 KDF 派生。
 *===============================================================*/
int clps_full_key_gen(uint8_t* pk,
                      uint8_t* sk,
                      const uint8_t* params,
                      const uint8_t* partialkey,
                      const uint8_t* id,
                      const uint8_t* usk,
                      const uint8_t* seed)
{
  uint8_t tr[TRBYTES], mu[CRHBYTES];
  uint8_t rho[SEEDBYTES];
  uint8_t cid[CTILDEBYTES], cid2[CTILDEBYTES];
  uint8_t w1buf[K * POLYW1_PACKEDBYTES];
  uint8_t xi[SEEDBYTES];
  uint8_t theta[CRHBYTES];
  uint8_t tr_i[TRBYTES];
  unsigned int i;
  keccak_state state;
  polyvecl mat[K], s1, z, s1hat;
  polyveck s2, tkgc, t, az, az1, az0;
  poly cp;

  /* 解包 params -> rho, t_KGC; 解包 partialkey -> ~c_i, z_i */
  clps_unpack_params(rho, &tkgc, params);
  clps_unpack_partialkey(cid, &z, partialkey);

  /* 先检查 ||z_i||_inf < GAMMA1 - BETA (z 仍为标准域) */
  if (polyvecl_chknorm(&z, GAMMA1 - BETA))
    return -1;

  /* tr = H1(rho || t_KGC) */
  clps_h1_tr_params(tr, params);

  /* mu_i = H2(tr || ID_i) */
  shake256_init(&state);
  shake256_absorb(&state, tr, TRBYTES);
  shake256_absorb(&state, id, SEEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

  polyvec_matrix_expand(mat, rho);

  /* c_i = SampleInBall(~c_i) */
  poly_challenge(&cp, cid);
  poly_ntt(&cp);

  /* w'1_i = HighBits(A*z_i - c_i*t_KGC, 2*gamma2) */
  polyvecl_ntt(&z);
  polyvec_matrix_pointwise_montgomery(&az, mat, &z);
  polyveck_reduce(&az);
  polyveck_invntt_tomont(&az);
  {
    polyveck ct;
    polyveck_ntt(&tkgc);
    polyveck_pointwise_poly_montgomery(&ct, &cp, &tkgc);
    polyveck_reduce(&ct);
    polyveck_invntt_tomont(&ct);
    polyveck_sub(&az, &az, &ct);
    polyveck_reduce(&az);
  }
  polyveck_caddq(&az);
  polyveck_decompose(&az1, &az0, &az);
  polyveck_pack_w1(w1buf, &az1);

  /* ~c'_i = H3(mu_i || w'packi) */
  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, w1buf, K * POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(cid2, CTILDEBYTES, &state);

  /* 校验 ~c'_i = ~c_i */
  for (i = 0; i < CTILDEBYTES; ++i)
    if (cid[i] != cid2[i])
      return -1;

  /* xi_i <- {0,1}^256 (从 seed 派生, 域名标签 0x5A) */
  {
    uint8_t seed_in[SEEDBYTES + 1];
    memcpy(seed_in, seed, SEEDBYTES);
    seed_in[SEEDBYTES] = 0x5A;
    shake256(xi, SEEDBYTES, seed_in, SEEDBYTES + 1);
  }

  /* theta_i = H2(xi_i || usk_i) */
  shake256_init(&state);
  shake256_absorb(&state, xi, SEEDBYTES);
  shake256_absorb(&state, usk, SEEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(theta, CRHBYTES, &state);

  /* (s1_i, s2_i) = ExpandS(theta_i) */
  polyvecl_uniform_eta(&s1, theta, 0);
  polyveck_uniform_eta(&s2, theta, L);

  /* t_i = A*s1_i + s2_i */
  s1hat = s1;
  polyvecl_ntt(&s1hat);
  polyvec_matrix_pointwise_montgomery(&t, mat, &s1hat);
  polyveck_reduce(&t);
  polyveck_invntt_tomont(&t);
  polyveck_add(&t, &t, &s2);
  polyveck_caddq(&t);

  /* tr_i = H1(rho || t_i) */
  clps_h1_tr(tr_i, rho, &t);

  /* sk_i := (s1_i, s2_i, usk_i, ~c_i, z_i, tr_i) */
  clps_pack_sk(sk, &s1, &s2, usk, cid, &z, tr_i);
  /* pk_i := (t_i, ~c_i) */
  clps_pack_pk(pk, &t, cid);

  return 0;
}


/*===============================================================
 * clps_principal_authorize_agent —— Principal-Authorize-Agent:
 *                                    授权代理凭证生成
 *
 * 伪代码对应 (Algorithm Principal-Authotize-Agent):
 *   c_P = SampleInBall(~c_P)
 *   tr_A = H1(rho || t_A)
 *   mu_P,A = H2(tr || ID_P || ID_A || tr_A || tr_P)
 *   nu = H2(usk_P || mu_P,A)
 *   循环:
 *     y_P,A = ExpandMask(nu, kappa)
 *     w_P,A = A*y_P,A
 *     ~c_P,A = H3(mu_P,A || Packw1(HighBits(w_P,A)))
 *     c_P,A = SampleInBall(~c_P,A)
 *     z_P,A = y_P,A + c_P,A*s1_P
 *     r0_P,A = LowBits(w_P,A - c_P,A*s2_P, 2*gamma2)
 *   proxykey_P,A := (ID_P, ID_A, z_P, ~c_P,A, z_P,A)
 *
 * 注: tex 中 w'1_P := HighBits(A*z_P - c_P*t_KGC) 已计算但后续未使用,
 *     本实现省略该无效计算。
 *===============================================================*/
int clps_principal_authorize_agent(uint8_t* proxykey,
                                   const uint8_t* params,
                                   const uint8_t* skp,
                                   const uint8_t* idp,
                                   const uint8_t* ida,
                                   const uint8_t* pka)
{
  uint8_t tr[TRBYTES], tra[TRBYTES], trp[TRBYTES];
  uint8_t rho[SEEDBYTES];
  uint8_t mu[CRHBYTES], nu[CRHBYTES];
  uint8_t cidp[CTILDEBYTES], cidpa[CTILDEBYTES];
  uint8_t cid_unused[CTILDEBYTES];
  uint8_t uskp[SEEDBYTES];
  uint8_t w1buf[K * POLYW1_PACKEDBYTES];
  uint16_t nonce = 0;
  unsigned int rej_count = 0;
  keccak_state state;
  polyvecl mat[K], s1p, ypa, zpa, zp;
  polyveck s2p, tkgc, ta, w, w1, w0;
  poly cp;
  polyvecl yhat;

  /* 解包授权者私钥 sk_P -> (s1_P, s2_P, usk_P, ~c_P, z_P, tr_P) */
  clps_unpack_sk(&s1p, &s2p, uskp, cidp, &zp, trp, skp);

  /* 解包代理者公钥 pk_A -> t_A (~c_A 不使用) */
  clps_unpack_pk(&ta, cid_unused, pka);

  /* 解包 params -> rho, t_KGC */
  clps_unpack_params(rho, &tkgc, params);

  polyvec_matrix_expand(mat, rho);

  /* tr = H1(rho || t_KGC) */
  clps_h1_tr_params(tr, params);

  /* tr_A = H1(rho || t_A) */
  clps_h1_tr(tra, rho, &ta);

  /* mu_P,A = H2(tr || ID_P || ID_A || tr_A || tr_P) */
  shake256_init(&state);
  shake256_absorb(&state, tr, TRBYTES);
  shake256_absorb(&state, idp, SEEDBYTES);
  shake256_absorb(&state, ida, SEEDBYTES);
  shake256_absorb(&state, tra, TRBYTES);
  shake256_absorb(&state, trp, TRBYTES);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

  /* nu = H2(usk_P || mu_P,A) */
  shake256_init(&state);
  shake256_absorb(&state, uskp, SEEDBYTES);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_finalize(&state);
  shake256_squeeze(nu, CRHBYTES, &state);

  polyvecl_ntt(&s1p);
  polyveck_ntt(&s2p);

rej:
  rej_count++;
  /* 采样遮掩向量 y_P,A */
  polyvecl_uniform_gamma1(&ypa, nu, nonce);
  nonce = nonce + L;

  /* w_P,A = A*y_P,A */
  yhat = ypa;
  polyvecl_ntt(&yhat);
  polyvec_matrix_pointwise_montgomery(&w, mat, &yhat);
  polyveck_reduce(&w);
  polyveck_invntt_tomont(&w);

  /* 分解 w_P,A -> (w1, w0) */
  polyveck_caddq(&w);
  polyveck_decompose(&w1, &w0, &w);
  polyveck_pack_w1(w1buf, &w1);

  /* ~c_P,A = H3(mu_P,A || w_packP,A) */
  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, w1buf, K * POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(cidpa, CTILDEBYTES, &state);

  /* c_P,A = SampleInBall(~c_P,A) */
  poly_challenge(&cp, cidpa);
  poly_ntt(&cp);

  /* z_P,A = y_P,A + c_P,A*s1_P */
  {
    polyvecl cs1;
    polyvecl_pointwise_poly_montgomery(&cs1, &cp, &s1p);
    polyvecl_invntt_tomont(&cs1);
    polyvecl_add(&zpa, &ypa, &cs1);
    polyvecl_reduce(&zpa);
  }

  /* 拒绝采样检查 #1: ||z_P,A||_inf < GAMMA1 - BETA */
  if (polyvecl_chknorm(&zpa, GAMMA1 - BETA))
    goto rej;

  /* r0_P,A = LowBits(w_P,A - c_P,A*s2_P) */
  {
    polyveck cs2;
    polyveck_pointwise_poly_montgomery(&cs2, &cp, &s2p);
    polyveck_reduce(&cs2);
    polyveck_invntt_tomont(&cs2);
    polyveck_sub(&w0, &w0, &cs2);
    polyveck_reduce(&w0);
  }

  /* 拒绝采样检查 #2: ||r0_P,A||_inf < GAMMA2 - BETA */
  if (polyveck_chknorm(&w0, GAMMA2 - BETA))
    goto rej;

  /* 打包 proxykey := (ID_P || ID_A || z_P || ~c_P,A || z_P,A) */
  clps_pack_proxykey(proxykey, idp, ida, &zp, cidpa, &zpa);

#ifndef CLPS_NO_REJ_PRINTF
  printf("[Principal-Authorize-Agent] sampling iterations (rej loop): %u\n", rej_count);
#endif
  clps_rej_authorize += rej_count;

  return 0;
}


/*===============================================================
 * clps_agent_sign —— Agent-Sign: 代理签名
 *
 * 伪代码对应 (Algorithm Agent-Sign):
 *   c_P = SampleInBall(~c_P)
 *   c_P,A = SampleInBall(~c_P,A)
 *   mu_A = H2(tr_A || M)
 *   nu' = H2(usk_A || mu_A || proxykey_P,A)
 *   w'1_P = HighBits(A*z_P - c_P*t_KGC, 2*gamma2)
 *   循环:
 *     y_A = ExpandMask(nu', kappa)
 *     w_A = A*y_A
 *     w1_A = HighBits(w_A, 2*gamma2) + w'1_P
 *     ~c_A = H3(mu_A || Packw1(w1_A))
 *     c_A = SampleInBall(~c_A)
 *     z_A = y_A + c_A*s1_A
 *     r0_A = LowBits(w_A - c_A*s2_A, 2*gamma2)
 *   z_0_A := z_A + z_P,A + z_P
 *   h_w := Makeh_w(A*z_P - c_P*t_KGC, A*z_P,A - c_P,A*t_P,
 *                  w_A - c_A*s2_A, 2*gamma2)
 *   sigma_P,A := (~c_A, ~c_P,A, z_P,A, z_0_A, h_w)
 *===============================================================*/
int clps_agent_sign(uint8_t* sig,
                    const uint8_t* m,
                    const uint8_t mlen,
                    const uint8_t* params,
                    const uint8_t* proxykey,
                    const uint8_t* pkp,
                    const uint8_t* sk)
{
  uint8_t mu[CRHBYTES], nu2[CRHBYTES];
  uint8_t tra[TRBYTES];
  uint8_t rho[SEEDBYTES];
  uint8_t cidp[CTILDEBYTES], cidpa[CTILDEBYTES], cida[CTILDEBYTES];
  uint8_t cidsk[CTILDEBYTES];
  uint8_t usk[SEEDBYTES];
  uint8_t idp_in[SEEDBYTES], ida_in[SEEDBYTES];
  uint8_t w1buf[K * POLYW1_PACKEDBYTES];
  uint8_t w1idpbuf[K * POLYW1_PACKEDBYTES];
  uint16_t nonce = 0;
  unsigned int rej_count = 0;
  keccak_state state;
  polyvecl mat[K], s1, y, zpa, zp, zsk, za, z0a;
  polyveck s2, tkgc, tp;
  polyveck w_full, w1a, wa1, w0a;
  polyveck az_p, b, w_minus_cs2, hw;
  poly cp, cpa, cidp_poly;
  polyvecl yhat;

  /* 解包授权者公钥 pk_P -> (t_P, ~c_P) */
  clps_unpack_pk(&tp, cidp, pkp);

  /* 解包代理凭证 proxykey -> (ID_P, ID_A, z_P, ~c_P,A, z_P,A)
   * 原始 proxykey 字节保留用于 nu' = H2(usk_A || mu_A || proxykey) */
  clps_unpack_proxykey(idp_in, ida_in, &zp, cidpa, &zpa, proxykey);

  /* 解包代理者私钥 sk_A -> (s1_A, s2_A, usk_A, ~c_A, z_A, tr_A)
   * 其中 ~c_A, z_A 部分在签名时不需要 */
  clps_unpack_sk(&s1, &s2, usk, cidsk, &zsk, tra, sk);

  /* 解包 params -> rho, t_KGC */
  clps_unpack_params(rho, &tkgc, params);

  polyvec_matrix_expand(mat, rho);

  /* c_P = SampleInBall(~c_P), c_P,A = SampleInBall(~c_P,A) */
  poly_challenge(&cidp_poly, cidp);
  poly_ntt(&cidp_poly);
  poly_challenge(&cpa, cidpa);
  poly_ntt(&cpa);

  /* mu_A = H2(tr_A || M) */
  shake256_init(&state);
  shake256_absorb(&state, tra, TRBYTES);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

  /* nu' = H2(usk_A || mu_A || proxykey_P,A) */
  shake256_init(&state);
  shake256_absorb(&state, usk, SEEDBYTES);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, proxykey, CLPS_PROXYKEYBYTES);
  shake256_finalize(&state);
  shake256_squeeze(nu2, CRHBYTES, &state);

  /* 预计算 a = A*z_P - c_P*t_KGC (reduce 后), 并分解得到 w'1_P */
  {
    polyvecl zp_hat;
    polyveck ct;
    polyveck tmp1, tmp0;

    zp_hat = zp;
    polyvecl_ntt(&zp_hat);
    polyvec_matrix_pointwise_montgomery(&az_p, mat, &zp_hat);
    polyveck_reduce(&az_p);
    polyveck_invntt_tomont(&az_p);

    polyveck_ntt(&tkgc);
    polyveck_pointwise_poly_montgomery(&ct, &cidp_poly, &tkgc);
    polyveck_reduce(&ct);
    polyveck_invntt_tomont(&ct);
    polyveck_sub(&az_p, &az_p, &ct);
    polyveck_reduce(&az_p);

    polyveck_caddq(&az_p);
    polyveck_decompose(&tmp1, &tmp0, &az_p);
    polyveck_pack_w1(w1idpbuf, &tmp1);
  }

  /* 预计算 b = A*z_P,A - c_P,A*t_P (reduce 后, Makeh_w 第二输入) */
  {
    polyvecl zpa_hat;
    polyveck ct;

    zpa_hat = zpa;
    polyvecl_ntt(&zpa_hat);
    polyvec_matrix_pointwise_montgomery(&b, mat, &zpa_hat);
    polyveck_reduce(&b);
    polyveck_invntt_tomont(&b);

    polyveck_ntt(&tp);
    polyveck_pointwise_poly_montgomery(&ct, &cpa, &tp);
    polyveck_reduce(&ct);
    polyveck_invntt_tomont(&ct);
    polyveck_sub(&b, &b, &ct);
    polyveck_reduce(&b);
  }

  /* s1_A, s2_A 转换到 NTT 域 */
  polyvecl_ntt(&s1);
  polyveck_ntt(&s2);

rej:
  rej_count++;
  /* 采样遮掩向量 y_A */
  polyvecl_uniform_gamma1(&y, nu2, nonce);
  nonce = nonce + L;

  /* w_A = A*y_A (保留完整 w_A 用于 w_A - c_A*s2_A) */
  yhat = y;
  polyvecl_ntt(&yhat);
  polyvec_matrix_pointwise_montgomery(&w_full, mat, &yhat);
  polyveck_reduce(&w_full);
  polyveck_invntt_tomont(&w_full);

  /* w1_A_solo = HighBits(w_A) */
  {
    polyveck w_tmp, w_copy;
    w_copy = w_full;
    polyveck_caddq(&w_copy);
    polyveck_decompose(&wa1, &w_tmp, &w_copy);
  }

  /* w1_A = HighBits(w_A) + w'1_P, 打包 w1buf */
  sclls_polyveck_unpack_w1(&w1a, w1idpbuf);
  polyveck_add(&w1a, &w1a, &wa1);
  polyveck_pack_w1(w1buf, &w1a);

  /* ~c_A = H3(mu_A || w1_A) */
  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, w1buf, K * POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(cida, CTILDEBYTES, &state);

  /* c_A = SampleInBall(~c_A) */
  poly_challenge(&cp, cida);
  poly_ntt(&cp);

  /* z_A = y_A + c_A*s1_A */
  {
    polyvecl cs1;
    polyvecl_pointwise_poly_montgomery(&cs1, &cp, &s1);
    polyvecl_invntt_tomont(&cs1);
    polyvecl_add(&za, &y, &cs1);
    polyvecl_reduce(&za);
  }

  /* 拒绝采样检查 #1: ||z_A||_inf < GAMMA1 - BETA */
  if (polyvecl_chknorm(&za, GAMMA1 - BETA))
    goto rej;

  /* r0_A = LowBits(w_A - c_A*s2_A) */
  {
    polyveck cs2tmp;
    polyveck_pointwise_poly_montgomery(&cs2tmp, &cp, &s2);
    polyveck_reduce(&cs2tmp);
    polyveck_invntt_tomont(&cs2tmp);
    polyveck_sub(&w_minus_cs2, &w_full, &cs2tmp);
    polyveck_reduce(&w_minus_cs2);
    polyveck_caddq(&w_minus_cs2);
    polyveck_decompose(&wa1, &w0a, &w_minus_cs2);
  }

  /* 拒绝采样检查 #2: ||r0_A||_inf < GAMMA2 - BETA */
  if (polyveck_chknorm(&w0a, GAMMA2 - BETA))
    goto rej;

  /* z_0_A = z_A + z_P,A + z_P */
  polyvecl_add(&z0a, &za, &zpa);
  polyvecl_add(&z0a, &z0a, &zp);

  /* h_w = Makeh_w(a, b, c), 其中 c = w_A - c_A*s2_A (已 caddq) */
  clps_polyveck_make_hw3(&hw, &az_p, &b, &w_minus_cs2);

  /* 打包签名 sigma := (~c_A || ~c_P,A || z_P,A || z_0_A || h_w) */
  clps_pack_sigma(sig, cida, cidpa, &zpa, &z0a, &hw);

#ifndef CLPS_NO_REJ_PRINTF
  printf("[Agent-Sign] sampling iterations (rej loop): %u\n", rej_count);
#endif
  clps_rej_agent_sign += rej_count;

  return 0;
}


/*===============================================================
 * clps_agent_verify —— Agent-Verify: 代理签名验证
 *
 * 伪代码对应 (Algorithm Agent-Verify):
 *   c_P = SampleInBall(~c_P)
 *   c_P,A = SampleInBall(~c_P,A)
 *   c_A = SampleInBall(~c_A)
 *   tr_A = H1(rho || t_A)
 *   mu_A = H2(tr_A || M)
 *   tr_P = H1(rho || t_P)
 *   mu_P,A = H2(tr || ID_P || ID_A || tr_A || tr_P)
 *   w'1_P,A = HighBits(A*z_P,A - c_P,A*t_P, 2*gamma2)
 *   w'1_A = HighBits(A*z_0_A - c_P*t_KGC - c_P,A*t_P - c_A*t_A, 2*gamma2)
 *           - w'1_P,A - h_w
 *   ~c'_A = H3(mu_A || Packw1(w'1_A))
 *   ~c'_P,A = H3(mu_P,A || Packw1(w'1_P,A))
 *   若 ||z_0_A||_inf < 3*(gamma1-beta) 且 ~c_A=~c'_A 且 ~c_P,A=~c'_P,A: True
 *
 * 注: tr_P 未显式存储, 由 pk_P 重新计算 tr_P = H1(rho || t_P)。
 *===============================================================*/
int clps_agent_verify(const uint8_t* m,
                      const uint8_t mlen,
                      const uint8_t* sig,
                      const uint8_t* params,
                      const uint8_t* pkp,
                      const uint8_t* pka,
                      const uint8_t* ida,
                      const uint8_t* idp)
{
  uint8_t mu[CRHBYTES], mupa[CRHBYTES];
  uint8_t tr[TRBYTES], tra[TRBYTES], trp[TRBYTES];
  uint8_t rho[SEEDBYTES];
  uint8_t cidp[CTILDEBYTES], cidpa[CTILDEBYTES], cida[CTILDEBYTES];
  uint8_t cida2[CTILDEBYTES], cidpa2[CTILDEBYTES];
  uint8_t cid_unused[CTILDEBYTES];
  uint8_t w1buf[K * POLYW1_PACKEDBYTES];
  uint8_t w1pabuf[K * POLYW1_PACKEDBYTES];
  unsigned int i;
  keccak_state state;
  polyvecl mat[K], zpa, z0a;
  polyveck tp, ta, tkgc, az, az1, az0, hw, w1pa;
  poly cp, cpa, cidp_poly;

  /* 解包签名 sigma -> (~c_A, ~c_P,A, z_P,A, z_0_A, h_w) */
  clps_unpack_sigma(cida, cidpa, &zpa, &z0a, &hw, sig);

  /* 范数检查: ||z_0_A||_inf < 3*(GAMMA1-BETA) (z_0_A 仍为标准域) */
  if (polyvecl_chknorm(&z0a, CLPS_Z0A_BOUND))
    return -1;

  /* 解包公钥 pk_P -> (t_P, ~c_P), pk_A -> (t_A, ~c_A) */
  clps_unpack_pk(&tp, cidp, pkp);
  clps_unpack_pk(&ta, cid_unused, pka);

  /* 解包 params -> rho, t_KGC */
  clps_unpack_params(rho, &tkgc, params);

  polyvec_matrix_expand(mat, rho);

  /* tr = H1(rho || t_KGC), tr_P = H1(rho || t_P), tr_A = H1(rho || t_A) */
  clps_h1_tr_params(tr, params);
  clps_h1_tr(trp, rho, &tp);
  clps_h1_tr(tra, rho, &ta);

  /* mu_A = H2(tr_A || M) */
  shake256_init(&state);
  shake256_absorb(&state, tra, TRBYTES);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

  /* mu_P,A = H2(tr || ID_P || ID_A || tr_A || tr_P) */
  shake256_init(&state);
  shake256_absorb(&state, tr, TRBYTES);
  shake256_absorb(&state, idp, SEEDBYTES);
  shake256_absorb(&state, ida, SEEDBYTES);
  shake256_absorb(&state, tra, TRBYTES);
  shake256_absorb(&state, trp, TRBYTES);
  shake256_finalize(&state);
  shake256_squeeze(mupa, CRHBYTES, &state);

  /* 三个挑战多项式: c_P, c_P,A, c_A */
  poly_challenge(&cidp_poly, cidp);
  poly_ntt(&cidp_poly);
  poly_challenge(&cpa, cidpa);
  poly_ntt(&cpa);
  poly_challenge(&cp, cida);
  poly_ntt(&cp);

  /* w'1_P,A = HighBits(A*z_P,A - c_P,A*t_P, 2*gamma2) */
  {
    polyvecl zpa_hat;
    polyveck ct;
    zpa_hat = zpa;
    polyvecl_ntt(&zpa_hat);
    polyvec_matrix_pointwise_montgomery(&az, mat, &zpa_hat);
    polyveck_reduce(&az);
    polyveck_invntt_tomont(&az);

    polyveck_ntt(&tp);
    polyveck_pointwise_poly_montgomery(&ct, &cpa, &tp);
    polyveck_reduce(&ct);
    polyveck_invntt_tomont(&ct);
    polyveck_sub(&az, &az, &ct);
    polyveck_reduce(&az);
  }
  polyveck_caddq(&az);
  polyveck_decompose(&az1, &az0, &az);
  polyveck_pack_w1(w1pabuf, &az1);

  /* ~c'_P,A = H3(mu_P,A || w'packP,A) */
  shake256_init(&state);
  shake256_absorb(&state, mupa, CRHBYTES);
  shake256_absorb(&state, w1pabuf, K * POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(cidpa2, CTILDEBYTES, &state);

  /* w'1_A = HighBits(A*z_0_A - c_P*t_KGC - c_P,A*t_P - c_A*t_A) - w'1_P,A - h_w */
  polyvecl_ntt(&z0a);
  polyvec_matrix_pointwise_montgomery(&az, mat, &z0a);
  polyveck_reduce(&az);
  polyveck_invntt_tomont(&az);

  /* 减去 c_P*t_KGC */
  {
    polyveck ct;
    polyveck_ntt(&tkgc);
    polyveck_pointwise_poly_montgomery(&ct, &cidp_poly, &tkgc);
    polyveck_reduce(&ct);
    polyveck_invntt_tomont(&ct);
    polyveck_sub(&az, &az, &ct);
    polyveck_reduce(&az);
  }

  /* 减去 c_P,A*t_P (t_P 已在 w'1_P,A 中转换到 NTT 域) */
  {
    polyveck ct;
    polyveck_pointwise_poly_montgomery(&ct, &cpa, &tp);
    polyveck_reduce(&ct);
    polyveck_invntt_tomont(&ct);
    polyveck_sub(&az, &az, &ct);
    polyveck_reduce(&az);
  }

  /* 减去 c_A*t_A */
  {
    polyveck ct;
    polyveck_ntt(&ta);
    polyveck_pointwise_poly_montgomery(&ct, &cp, &ta);
    polyveck_reduce(&ct);
    polyveck_invntt_tomont(&ct);
    polyveck_sub(&az, &az, &ct);
    polyveck_reduce(&az);
  }

  /* 分解并减去 w'1_P,A 与 h_w */
  polyveck_caddq(&az);
  polyveck_decompose(&az1, &az0, &az);
  sclls_polyveck_unpack_w1(&w1pa, w1pabuf);
  polyveck_sub(&az1, &az1, &w1pa);
  polyveck_sub(&az1, &az1, &hw);
  polyveck_pack_w1(w1buf, &az1);

  /* ~c'_A = H3(mu_A || w'packA) */
  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, w1buf, K * POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(cida2, CTILDEBYTES, &state);

  /* 比较 ~c_A = ~c'_A 且 ~c_P,A = ~c'_P,A */
  for (i = 0; i < CTILDEBYTES; ++i)
    if (cida[i] != cida2[i])
      return -1;
  for (i = 0; i < CTILDEBYTES; ++i)
    if (cidpa[i] != cidpa2[i])
      return -1;

  return 0;
}

/* 临时调试函数: 打印两串的差异 */
int clps_dbg_verify(const uint8_t* m,
                    const uint8_t mlen,
                    const uint8_t* sig,
                    const uint8_t* params,
                    const uint8_t* pkp,
                    const uint8_t* pka,
                    const uint8_t* ida,
                    const uint8_t* idp)
{
  uint8_t mu[CRHBYTES], mupa[CRHBYTES];
  uint8_t tr[TRBYTES], tra[TRBYTES], trp[TRBYTES];
  uint8_t rho[SEEDBYTES];
  uint8_t cidp[CTILDEBYTES], cidpa[CTILDEBYTES], cida[CTILDEBYTES];
  uint8_t cida2[CTILDEBYTES], cidpa2[CTILDEBYTES];
  uint8_t cid_unused[CTILDEBYTES];
  uint8_t w1buf[K * POLYW1_PACKEDBYTES];
  uint8_t w1pabuf[K * POLYW1_PACKEDBYTES];
  unsigned int i, j;
  int dif_ca = 0, dif_cpa = 0;
  keccak_state state;
  polyvecl mat[K], zpa, z0a;
  polyveck tp, ta, tkgc, az, az1, az0, hw, w1pa;
  poly cp, cpa, cidp_poly;
  int32_t min_a1 = 9999, max_a1 = -9999;
  int32_t min_hw = 9999, max_hw = -9999;

  clps_unpack_sigma(cida, cidpa, &zpa, &z0a, &hw, sig);

  for (i = 0; i < K; ++i)
    for (j = 0; j < N; ++j) {
      if (hw.vec[i].coeffs[j] < min_hw) min_hw = hw.vec[i].coeffs[j];
      if (hw.vec[i].coeffs[j] > max_hw) max_hw = hw.vec[i].coeffs[j];
    }
  printf("[DBG] z0a norm check bound=%d, hw min=%d max=%d\n", CLPS_Z0A_BOUND, min_hw, max_hw);

  clps_unpack_pk(&tp, cidp, pkp);
  clps_unpack_pk(&ta, cid_unused, pka);
  clps_unpack_params(rho, &tkgc, params);
  polyvec_matrix_expand(mat, rho);

  clps_h1_tr_params(tr, params);
  clps_h1_tr(trp, rho, &tp);
  clps_h1_tr(tra, rho, &ta);

  shake256_init(&state);
  shake256_absorb(&state, tra, TRBYTES);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

  shake256_init(&state);
  shake256_absorb(&state, tr, TRBYTES);
  shake256_absorb(&state, idp, SEEDBYTES);
  shake256_absorb(&state, ida, SEEDBYTES);
  shake256_absorb(&state, tra, TRBYTES);
  shake256_absorb(&state, trp, TRBYTES);
  shake256_finalize(&state);
  shake256_squeeze(mupa, CRHBYTES, &state);

  poly_challenge(&cidp_poly, cidp);
  poly_ntt(&cidp_poly);
  poly_challenge(&cpa, cidpa);
  poly_ntt(&cpa);
  poly_challenge(&cp, cida);
  poly_ntt(&cp);

  {
    polyvecl zpa_hat;
    polyveck ct;
    zpa_hat = zpa;
    polyvecl_ntt(&zpa_hat);
    polyvec_matrix_pointwise_montgomery(&az, mat, &zpa_hat);
    polyveck_reduce(&az);
    polyveck_invntt_tomont(&az);
    polyveck_ntt(&tp);
    polyveck_pointwise_poly_montgomery(&ct, &cpa, &tp);
    polyveck_reduce(&ct);
    polyveck_invntt_tomont(&ct);
    polyveck_sub(&az, &az, &ct);
    polyveck_reduce(&az);
  }
  polyveck_caddq(&az);
  polyveck_decompose(&az1, &az0, &az);
  polyveck_pack_w1(w1pabuf, &az1);

  shake256_init(&state);
  shake256_absorb(&state, mupa, CRHBYTES);
  shake256_absorb(&state, w1pabuf, K * POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(cidpa2, CTILDEBYTES, &state);

  polyvecl_ntt(&z0a);
  polyvec_matrix_pointwise_montgomery(&az, mat, &z0a);
  polyveck_reduce(&az);
  polyveck_invntt_tomont(&az);

  {
    polyveck ct;
    polyveck_ntt(&tkgc);
    polyveck_pointwise_poly_montgomery(&ct, &cidp_poly, &tkgc);
    polyveck_reduce(&ct);
    polyveck_invntt_tomont(&ct);
    polyveck_sub(&az, &az, &ct);
    polyveck_reduce(&az);
  }
  {
    polyveck ct;
    polyveck_pointwise_poly_montgomery(&ct, &cpa, &tp);
    polyveck_reduce(&ct);
    polyveck_invntt_tomont(&ct);
    polyveck_sub(&az, &az, &ct);
    polyveck_reduce(&az);
  }
  {
    polyveck ct;
    polyveck_ntt(&ta);
    polyveck_pointwise_poly_montgomery(&ct, &cp, &ta);
    polyveck_reduce(&ct);
    polyveck_invntt_tomont(&ct);
    polyveck_sub(&az, &az, &ct);
    polyveck_reduce(&az);
  }

  polyveck_caddq(&az);
  polyveck_decompose(&az1, &az0, &az);
  sclls_polyveck_unpack_w1(&w1pa, w1pabuf);
  polyveck_sub(&az1, &az1, &w1pa);
  polyveck_sub(&az1, &az1, &hw);
  for (i = 0; i < K; ++i)
    for (j = 0; j < N; ++j) {
      if (az1.vec[i].coeffs[j] < min_a1) min_a1 = az1.vec[i].coeffs[j];
      if (az1.vec[i].coeffs[j] > max_a1) max_a1 = az1.vec[i].coeffs[j];
    }
  printf("[DBG] w'1_A coeff min=%d max=%d (expected in [0,43])\n", min_a1, max_a1);
  polyveck_pack_w1(w1buf, &az1);

  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, w1buf, K * POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(cida2, CTILDEBYTES, &state);

  for (i = 0; i < CTILDEBYTES; ++i)
    if (cida[i] != cida2[i]) dif_ca = 1;
  for (i = 0; i < CTILDEBYTES; ++i)
    if (cidpa[i] != cidpa2[i]) dif_cpa = 1;

  printf("[DBG] cida match=%d, cidpa match=%d\n", !dif_ca, !dif_cpa);
  if (dif_ca) {
    printf("[DBG] cida : ");
    for (i = 0; i < CTILDEBYTES; ++i) printf("%02x", cida[i]);
    printf("\n[DBG] cida2: ");
    for (i = 0; i < CTILDEBYTES; ++i) printf("%02x", cida2[i]);
    printf("\n");
  }
  if (dif_cpa) {
    printf("[DBG] cidpa : ");
    for (i = 0; i < CTILDEBYTES; ++i) printf("%02x", cidpa[i]);
    printf("\n[DBG] cidpa2: ");
    for (i = 0; i < CTILDEBYTES; ++i) printf("%02x", cidpa2[i]);
    printf("\n");
  }
  return 0;
}
