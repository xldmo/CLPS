/**
 * @file bench.c
 * @brief CLPS 方案基准测试驱动 (拒绝采样筛选 + 平均耗时统计)
 *
 * 用法:
 *   bench scan <start> <count> <rej_thresh> [verbose]
 *       从下标 start 开始连续检测 count 个随机数, 统计每个随机数驱动
 *       完整 CLPS 流程时四次拒绝采样调用的迭代次数, 输出:
 *         F,<idx>,<failmask>                          算法失败的随机数
 *         Q/V,<idx>,<rej_partial_P>,<rej_partial_A>,<rej_authorize>,<rej_sign>,<total>,<fail>
 *             各子算法调用迭代次数均<=thresh
 *         S,<scanned>,<fail>,<qualified>              汇总
 *
 *   bench time <idxlist_file> <reps>
 *       读取合格随机数下标列表(每行一个), 每个下标运行 reps 次并计时, 输出:
 *         T,<idx>,<avg_total_us>,<avg_rej_total>
 *         A,<runs>,<avg_total_us>,<avg_setup_us>,<avg_partial_pair_us>,
 *           <avg_fullkey_pair_us>,<avg_authorize_us>,<avg_agentsign_us>,
 *           <avg_verify_us>,<avg_rej_partial_P>,<avg_rej_partial_A>,
 *           <avg_rej_authorize>,<avg_rej_sign>,<avg_rej_total>,<failed_runs>
 *
 * 随机数定义: 32 字节主随机数 R (由下标经 shake256 确定性派生);
 *   R 直接作为 Setup 的种子, 其余随机量由 R 经域标签派生,
 *   使得同一个随机数在两个方案/三个安全等级下都表示 "同一个随机数"。
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif

#include "params.h"
#include "clps.h"
#include "clps_packing.h"
#include "fips202.h"

/* 由 clps.c 导出的拒绝采样计数器 */
extern unsigned long long clps_rej_partial_keygen;
extern unsigned long long clps_rej_authorize;
extern unsigned long long clps_rej_agent_sign;

/*---------------------------------------------------------------
 * 计时
 *---------------------------------------------------------------*/
static int64_t now_ns(void)
{
#ifdef _WIN32
    LARGE_INTEGER freq, counter;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&counter);
    return (int64_t)(counter.QuadPart * 1000000 / freq.QuadPart);
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (int64_t)ts.tv_sec * 1000000000LL + ts.tv_nsec;
#endif
}

/*---------------------------------------------------------------
 * 随机数派生
 *---------------------------------------------------------------*/
/* 由下标确定性派生 32 字节主随机数 R (两个方案共用同一规则) */
static void master_from_index(uint8_t master[SEEDBYTES], unsigned long long idx)
{
    uint8_t in[8 + 8];
    unsigned int i;
    for (i = 0; i < 8; ++i)
        in[i] = (uint8_t)((idx >> (8 * i)) & 0xFF);
    memcpy(in + 8, "masterRN", 8);
    shake256(master, SEEDBYTES, in, sizeof(in));
}

/* 由主随机数 R 经域标签派生子随机量 */
static void derive(uint8_t out[SEEDBYTES], const uint8_t master[SEEDBYTES], const char *tag)
{
    keccak_state st;
    shake256_init(&st);
    shake256_absorb(&st, master, SEEDBYTES);
    shake256_absorb(&st, (const uint8_t *)tag, (size_t)strlen(tag));
    shake256_finalize(&st);
    shake256_squeeze(out, SEEDBYTES, &st);
}

/*---------------------------------------------------------------
 * 单次运行统计
 *---------------------------------------------------------------*/
typedef struct {
    unsigned int rej[4];        /* PartialKeyGen(P/A) / Authorize / AgentSign */
    unsigned long long total;   /* 四次拒绝采样循环之和 */
    int fail;                   /* bit0: FullKeyGen(P) 失败, bit1: FullKeyGen(A) 失败,
                                 * bit2: AgentVerify 失败, bit3: dbg_verify 失败 */
    int64_t us[8];              /* Setup / PartialKeyGen(P) / PartialKeyGen(A) /
                                 * FullKeyGen(P) / FullKeyGen(A) / Authorize /
                                 * AgentSign / AgentVerify 耗时(us) */
    int64_t us_total;           /* 全流程耗时(us) */
} stat_t;

/* 固定输入 (与 main.c 保持一致) */
static const uint8_t g_msg[SEEDBYTES] = {
    0xD5, 0xC2, 0x93, 0xCE, 0x6B, 0xDE, 0xD2, 0x1E,
    0x3A, 0xAD, 0x35, 0x61, 0x00, 0xCA, 0x8D, 0x6E,
    0x42, 0xE3, 0x94, 0x50, 0xCB, 0x1C, 0x43, 0x5B,
    0x7D, 0xD8, 0x02, 0xEA, 0x43, 0x74, 0x65, 0x32
};

static void run_once(const uint8_t master[SEEDBYTES], stat_t *st)
{
    uint8_t params[CLPS_PARAMSBYTES];
    uint8_t msk[CLPS_MSKBYTES];
    uint8_t partialkey_p[CLPS_PARTIALKEYBYTES];
    uint8_t partialkey_a[CLPS_PARTIALKEYBYTES];
    uint8_t pkp[CLPS_PKIDBYTES], skp[CLPS_SKIDBYTES];
    uint8_t pka[CLPS_PKIDBYTES], ska[CLPS_SKIDBYTES];
    uint8_t seed_p[SEEDBYTES], seed_a[SEEDBYTES];
    uint8_t uskp[SEEDBYTES], uska[SEEDBYTES];
    uint8_t proxykey[CLPS_PROXYKEYBYTES];
    uint8_t sig[CLPS_SIGBYTES];
    uint8_t idp[SEEDBYTES], ida[SEEDBYTES];
    int r;
    int64_t t0, t1, t2;

    memset(st, 0, sizeof(*st));

    memset(idp, 0x01, SEEDBYTES);
    memset(ida, 0x02, SEEDBYTES);
    derive(seed_p, master, "CLPS.seed_P");
    derive(seed_a, master, "CLPS.seed_A");
    derive(uskp,   master, "CLPS.usk_P");
    derive(uska,   master, "CLPS.usk_A");

    t0 = now_ns();

    /* Setup */
    t1 = now_ns();
    clps_setup(params, msk, master);
    t2 = now_ns();
    st->us[0] = t2 - t1;

    /* PartialKeyGen (Principal) */
    clps_rej_partial_keygen = 0;
    t1 = now_ns();
    r = clps_partial_key_gen(partialkey_p, params, msk, idp);
    t2 = now_ns();
    st->us[1] = t2 - t1;
    st->rej[0] = (unsigned int)clps_rej_partial_keygen;
    if (r != 0)
        st->fail |= 0x08;

    /* PartialKeyGen (Agent) */
    clps_rej_partial_keygen = 0;
    t1 = now_ns();
    r = clps_partial_key_gen(partialkey_a, params, msk, ida);
    t2 = now_ns();
    st->us[2] = t2 - t1;
    st->rej[1] = (unsigned int)clps_rej_partial_keygen;
    if (r != 0)
        st->fail |= 0x10;

    /* FullKeyGen (Principal) */
    t1 = now_ns();
    r = clps_full_key_gen(pkp, skp, params, partialkey_p, idp, uskp, seed_p);
    t2 = now_ns();
    st->us[3] = t2 - t1;
    if (r != 0)
        st->fail |= 0x01;

    /* FullKeyGen (Agent) */
    t1 = now_ns();
    r = clps_full_key_gen(pka, ska, params, partialkey_a, ida, uska, seed_a);
    t2 = now_ns();
    st->us[4] = t2 - t1;
    if (r != 0)
        st->fail |= 0x02;

    /* 授权凭证 */
    clps_rej_authorize = 0;
    t1 = now_ns();
    r = clps_principal_authorize_agent(proxykey, params, skp, idp, ida, pka);
    t2 = now_ns();
    st->us[5] = t2 - t1;
    st->rej[2] = (unsigned int)clps_rej_authorize;
    if (r != 0)
        st->fail |= 0x20;

    /* 代理签名 */
    clps_rej_agent_sign = 0;
    t1 = now_ns();
    r = clps_agent_sign(sig, g_msg, SEEDBYTES, params, proxykey, pkp, ska);
    t2 = now_ns();
    st->us[6] = t2 - t1;
    st->rej[3] = (unsigned int)clps_rej_agent_sign;
    if (r != 0)
        st->fail |= 0x40;

    /* 验证 */
    t1 = now_ns();
    r = clps_agent_verify(g_msg, SEEDBYTES, sig, params, pkp, pka, ida, idp);
    t2 = now_ns();
    st->us[7] = t2 - t1;
    if (r != 0)
        st->fail |= 0x04;

    st->us_total = now_ns() - t0;

    st->total  = (unsigned long long)st->rej[0] + st->rej[1] +
                 st->rej[2] + st->rej[3];
}


int main(int argc, char **argv) {
    if (argc != 2) return 2;
    FILE *fp = fopen(argv[1], "r");
    if (!fp) return 2;
    unsigned long long idx;
    unsigned count=0, failures=0;
    uint8_t master[SEEDBYTES];
    stat_t *records = calloc(1000, sizeof(stat_t));
    if (!records) return 2;
    fprintf(stderr, "mode=%d K=%d L=%d params=%d msk=%d partial=%d pk=%d sk=%d token=%d signature=%d\n",
        DILITHIUM_MODE,K,L,CLPS_PARAMSBYTES,CLPS_MSKBYTES,CLPS_PARTIALKEYBYTES,
        CLPS_PKIDBYTES,CLPS_SKIDBYTES,CLPS_PROXYKEYBYTES,CLPS_SIGBYTES);
    puts("level,index,rep,total_ns,setup_ns,partial_principal_ns,partial_agent_ns,fullkey_principal_ns,fullkey_agent_ns,authorize_ns,agentsign_ns,verify_ns,rej_partial_principal,rej_partial_agent,rej_authorize,rej_agentsign,failmask");
    while (fscanf(fp, "%llu", &idx)==1) {
        master_from_index(master, idx);
        stat_t warm;
        for (int k=0;k<10;k++) { run_once(master,&warm); if(warm.fail) return 3; }
        for (int k=0;k<1000;k++) run_once(master,&records[k]);
        for (int k=0;k<1000;k++) {
            stat_t *s=&records[k];
            printf("%d,%llu,%d,%lld",DILITHIUM_MODE,idx,k+1,(long long)s->us_total);
            for(int j=0;j<8;j++) printf(",%lld",(long long)s->us[j]);
            for(int j=0;j<4;j++) printf(",%u",s->rej[j]);
            printf(",%d\n",s->fail);
            failures += s->fail != 0;
        }
        count++; fflush(stdout);
        fprintf(stderr,"completed_seeds=%u failures=%u\n",count,failures);
    }
    fclose(fp); free(records);
    return count==100 && failures==0 ? 0 : 4;
}
