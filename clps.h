/**
 * @file clps.h
 * @brief CLPS (Certificateless Proxy Signature) 方案头文件
 *
 * 本方案基于 ML-DSA (Dilithium) 格密码框架，实现无证书代理签名机制。
 * 包含以下核心算法:
 *   1. Setup                    - 系统参数生成 (params, msk)
 *   2. PartialKeyGen            - 部分私钥生成 (KGC 签发)
 *   3. FullKeyGen               - 完整用户密钥对生成
 *   4. Principal-Authorize-Agent - 授权代理凭证生成 (proxykey)
 *   5. Agent-Sign               - 代理签名
 *   6. Agent-Verify             - 代理签名验证
 *
 * 参考文献:
 *   - CLPS_scheme_终版.tex (CLPS 伪代码, 唯一实现规格)
 *   - Dilithium / ML-DSA (基础格密码框架)
 */

#ifndef CLPS_H
#define CLPS_H

#include <stddef.h>
#include <stdint.h>
#include "params.h"
#include "polyvec.h"
#include "poly.h"

/* 命名空间宏 (与 config.h 保持一致) */
#ifndef CLPS_NAMESPACE
#define CLPS_NAMESPACE(s) clps_##s
#endif

/*===============================================================
 * CLPS 方案函数声明
 *===============================================================*/

/**
 * @brief Setup: 系统参数生成
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
 * @param params    输出公共参数 (CLPS_PARAMSBYTES 字节, rho||t_KGC)
 * @param msk       输出主私钥 (CLPS_MSKBYTES 字节, rho'||s1_KGC||s2_KGC)
 * @param seed      输入随机种子 (SEEDBYTES 字节, 派生 rho 与 rho')
 * @return          成功返回 0
 */
int clps_setup(uint8_t* params, uint8_t* msk, const uint8_t* seed);

/**
 * @brief PartialKeyGen: 部分私钥生成 (KGC 签发)
 *
 * 伪代码对应 (Algorithm PartialKeyGen):
 *   mu_i = H2(tr || ID_i)
 *   rho_i = H2(rho' || mu_i)
 *   (循环) y_i = ExpandMask(rho_i, kappa)
 *          w_i = A*y_i
 *          w1_i = HighBits(w_i, 2*gamma2)
 *          ~c_i = H3(mu_i || Packw1(w1_i))
 *          c_i = SampleInBall(~c_i)
 *          z_i = y_i + c_i*s1_KGC
 *          r0_i = LowBits(w_i - c_i*s2_KGC, 2*gamma2)
 *   partialkey_i := (~c_i, z_i)
 *
 * @param partialkey  输出部分私钥 (CLPS_PARTIALKEYBYTES 字节, ~c_i||z_i)
 * @param params      输入公共参数
 * @param msk         输入主私钥
 * @param id          输入用户身份信息 (SEEDBYTES 字节)
 * @return            成功返回 0
 */
int clps_partial_key_gen(uint8_t* partialkey,
                         const uint8_t* params,
                         const uint8_t* msk,
                         const uint8_t* id);

/**
 * @brief FullKeyGen: 完整用户密钥对生成
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
 * @param pk          输出公钥 (CLPS_PKIDBYTES 字节, t_i||~c_i)
 * @param sk          输出私钥 (CLPS_SKIDBYTES 字节)
 * @param params      输入公共参数
 * @param partialkey  输入部分私钥 (来自 KGC)
 * @param id          输入用户身份信息 (SEEDBYTES 字节)
 * @param usk         输入用户秘密值 (SEEDBYTES 字节)
 * @param seed        输入随机种子 (SEEDBYTES 字节, 派生 xi_i)
 * @return            成功返回 0, 部分私钥无效返回 -1
 */
int clps_full_key_gen(uint8_t* pk,
                      uint8_t* sk,
                      const uint8_t* params,
                      const uint8_t* partialkey,
                      const uint8_t* id,
                      const uint8_t* usk,
                      const uint8_t* seed);

/**
 * @brief Principal-Authorize-Agent: 授权代理凭证生成
 *
 * 伪代码对应 (Algorithm Principal-Authotize-Agent):
 *   c_P = SampleInBall(~c_P)
 *   tr_A = H1(rho || t_A)
 *   mu_P,A = H2(tr || ID_P || ID_A || tr_A || tr_P)
 *   nu = H2(usk_P || mu_P,A)
 *   (循环) y_P,A = ExpandMask(nu, kappa)
 *          w_P,A = A*y_P,A
 *          ~c_P,A = H3(mu_P,A || Packw1(HighBits(w_P,A)))
 *          z_P,A = y_P,A + c_P,A*s1_P
 *          r0_P,A = LowBits(w_P,A - c_P,A*s2_P, 2*gamma2)
 *   proxykey_P,A := (ID_P, ID_A, z_P, ~c_P,A, z_P,A)
 *
 * @param proxykey  输出代理凭证 (CLPS_PROXYKEYBYTES 字节)
 * @param params    输入公共参数
 * @param skp       输入授权者私钥 sk_P
 * @param idp       输入授权者身份 ID_P (SEEDBYTES 字节)
 * @param ida       输入代理者身份 ID_A (SEEDBYTES 字节)
 * @param pka       输入代理者公钥 pk_A
 * @return          成功返回 0
 */
int clps_principal_authorize_agent(uint8_t* proxykey,
                                   const uint8_t* params,
                                   const uint8_t* skp,
                                   const uint8_t* idp,
                                   const uint8_t* ida,
                                   const uint8_t* pka);

/**
 * @brief Agent-Sign: 代理签名
 *
 * 伪代码对应 (Algorithm Agent-Sign):
 *   c_P = SampleInBall(~c_P)
 *   c_P,A = SampleInBall(~c_P,A)
 *   mu_A = H2(tr_A || M)
 *   nu' = H2(usk_A || mu_A || proxykey_P,A)
 *   w'1_P = HighBits(A*z_P - c_P*t_KGC, 2*gamma2)
 *   (循环) y_A = ExpandMask(nu', kappa)
 *          w_A = A*y_A
 *          w1_A = HighBits(w_A, 2*gamma2) + w'1_P
 *          ~c_A = H3(mu_A || Packw1(w1_A))
 *          z_A = y_A + c_A*s1_A
 *          r0_A = LowBits(w_A - c_A*s2_A, 2*gamma2)
 *   z_0_A := z_A + z_P,A + z_P
 *   h_w := Makeh_w(A*z_P - c_P*t_KGC, A*z_P,A - c_P,A*t_P,
 *                  w_A - c_A*s2_A, 2*gamma2)
 *   sigma_P,A := (~c_A, ~c_P,A, z_P,A, z_0_A, h_w)
 *
 * @param sig       输出代理签名 (CLPS_SIGBYTES 字节)
 * @param m         输入消息
 * @param mlen      消息长度
 * @param params    输入公共参数
 * @param proxykey  输入代理凭证 (来自 Principal-Authorize-Agent)
 * @param pkp       输入授权者公钥 pk_P
 * @param sk        输入代理者私钥 sk_A
 * @return          成功返回 0
 */
int clps_agent_sign(uint8_t* sig,
                    const uint8_t* m,
                    const uint8_t mlen,
                    const uint8_t* params,
                    const uint8_t* proxykey,
                    const uint8_t* pkp,
                    const uint8_t* sk);

/**
 * @brief Agent-Verify: 代理签名验证
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
 * @param m         输入消息
 * @param mlen      消息长度
 * @param sig       输入代理签名
 * @param params    输入公共参数
 * @param pkp       输入授权者公钥 pk_P
 * @param pka       输入代理者公钥 pk_A
 * @param ida       输入代理者身份 ID_A
 * @param idp       输入授权者身份 ID_P
 * @return          验证成功返回 0, 失败返回 -1
 */
int clps_agent_verify(const uint8_t* m,
                      const uint8_t mlen,
                      const uint8_t* sig,
                      const uint8_t* params,
                      const uint8_t* pkp,
                      const uint8_t* pka,
                      const uint8_t* ida,
                      const uint8_t* idp);

/* 临时调试: 打印验证各阶段诊断信息 (仅测试用) */
int clps_dbg_verify(const uint8_t* m,
                    const uint8_t mlen,
                    const uint8_t* sig,
                    const uint8_t* params,
                    const uint8_t* pkp,
                    const uint8_t* pka,
                    const uint8_t* ida,
                    const uint8_t* idp);

#endif /* CLPS_H */
