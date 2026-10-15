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
static int64_t now_us(void)
{
#ifdef _WIN32
    LARGE_INTEGER freq, counter;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&counter);
    return (int64_t)(counter.QuadPart * 1000000 / freq.QuadPart);
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (int64_t)ts.tv_sec * 1000000 + ts.tv_nsec / 1000;
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

    t0 = now_us();

    /* Setup */
    t1 = now_us();
    clps_setup(params, msk, master);
    t2 = now_us();
    st->us[0] = t2 - t1;

    /* PartialKeyGen (Principal) */
    clps_rej_partial_keygen = 0;
    t1 = now_us();
    r = clps_partial_key_gen(partialkey_p, params, msk, idp);
    t2 = now_us();
    st->us[1] = t2 - t1;
    st->rej[0] = (unsigned int)clps_rej_partial_keygen;
    if (r != 0)
        st->fail |= 0x08;

    /* PartialKeyGen (Agent) */
    clps_rej_partial_keygen = 0;
    t1 = now_us();
    r = clps_partial_key_gen(partialkey_a, params, msk, ida);
    t2 = now_us();
    st->us[2] = t2 - t1;
    st->rej[1] = (unsigned int)clps_rej_partial_keygen;
    if (r != 0)
        st->fail |= 0x10;

    /* FullKeyGen (Principal) */
    t1 = now_us();
    r = clps_full_key_gen(pkp, skp, params, partialkey_p, idp, uskp, seed_p);
    t2 = now_us();
    st->us[3] = t2 - t1;
    if (r != 0)
        st->fail |= 0x01;

    /* FullKeyGen (Agent) */
    t1 = now_us();
    r = clps_full_key_gen(pka, ska, params, partialkey_a, ida, uska, seed_a);
    t2 = now_us();
    st->us[4] = t2 - t1;
    if (r != 0)
        st->fail |= 0x02;

    /* 授权凭证 */
    clps_rej_authorize = 0;
    t1 = now_us();
    r = clps_principal_authorize_agent(proxykey, params, skp, idp, ida, pka);
    t2 = now_us();
    st->us[5] = t2 - t1;
    st->rej[2] = (unsigned int)clps_rej_authorize;
    if (r != 0)
        st->fail |= 0x20;

    /* 代理签名 */
    clps_rej_agent_sign = 0;
    t1 = now_us();
    r = clps_agent_sign(sig, g_msg, SEEDBYTES, params, proxykey, pkp, ska);
    t2 = now_us();
    st->us[6] = t2 - t1;
    st->rej[3] = (unsigned int)clps_rej_agent_sign;
    if (r != 0)
        st->fail |= 0x40;

    /* 验证 */
    t1 = now_us();
    r = clps_agent_verify(g_msg, SEEDBYTES, sig, params, pkp, pka, ida, idp);
    t2 = now_us();
    st->us[7] = t2 - t1;
    if (r != 0)
        st->fail |= 0x04;

    st->us_total = now_us() - t0;

    st->total  = (unsigned long long)st->rej[0] + st->rej[1] +
                 st->rej[2] + st->rej[3];
}

/*---------------------------------------------------------------
 * scan 模式
 *---------------------------------------------------------------*/
static int do_scan(unsigned long long start, unsigned long long count,
                   unsigned long long thresh, int verbose)
{
    unsigned long long i;
    unsigned long long scanned = 0, nfail = 0, nqual = 0;
    uint8_t master[SEEDBYTES];
    stat_t st;

    printf("# mode=scan start=%llu count=%llu thresh=%llu level=%d\n",
           start, count, thresh, DILITHIUM_MODE);

    for (i = 0; i < count; ++i) {
        unsigned long long idx = start + i;
        master_from_index(master, idx);
        run_once(master, &st);
        scanned++;
        if (st.fail) {
            nfail++;
            printf("F,%llu,%d\n", idx, st.fail);
            fflush(stdout);
        }
        /* 合格判定: 验证成功且每次子算法调用的迭代次数都 <= thresh */
        if (st.fail == 0 && st.rej[0] <= thresh && st.rej[1] <= thresh &&
            st.rej[2] <= thresh && st.rej[3] <= thresh) {
            nqual++;
            printf("Q,%llu,%u,%u,%u,%u,%llu,%d\n",
                   idx, st.rej[0], st.rej[1], st.rej[2], st.rej[3],
                   st.total, st.fail);
            fflush(stdout);
        } else if (verbose) {
            printf("V,%llu,%u,%u,%u,%u,%llu,%d\n",
                   idx, st.rej[0], st.rej[1], st.rej[2], st.rej[3],
                   st.total, st.fail);
        }
    }

    printf("S,%llu,%llu,%llu\n", scanned, nfail, nqual);
    return 0;
}

/*---------------------------------------------------------------
 * time 模式
 *---------------------------------------------------------------*/
static int do_time(const char *listfile, int reps)
{
    FILE *fp = fopen(listfile, "r");
    unsigned long long idx;
    uint8_t master[SEEDBYTES];
    stat_t st;
    double sum_total = 0.0, sum_us[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    double sum_rej[4] = {0, 0, 0, 0};
    unsigned long long runs = 0, nfail = 0;
    double sum_rejcnt = 0.0;
    char line[256];

    if (!fp) {
        fprintf(stderr, "cannot open %s\n", listfile);
        return 1;
    }

    printf("# mode=time list=%s reps=%d level=%d\n", listfile, reps, DILITHIUM_MODE);

    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r')
            continue;
        idx = strtoull(line, NULL, 10);
        master_from_index(master, idx);

        {
            double avg_total = 0.0;
            double avg_rej = 0.0;
            int k;
            int64_t acc = 0;
            unsigned long long acc_rej = 0;
            for (k = 0; k < reps; ++k) {
                run_once(master, &st);
                acc += st.us_total;
                acc_rej += st.total;
                sum_total += (double)st.us_total;
                for (int j = 0; j < 8; ++j)
                    sum_us[j] += (double)st.us[j];
                for (int j = 0; j < 4; ++j)
                    sum_rej[j] += (double)st.rej[j];
                sum_rejcnt += (double)st.total;
                if (st.fail)
                    nfail++;
                runs++;
            }
            avg_total = (double)acc / reps;
            avg_rej = (double)acc_rej / reps;
            printf("T,%llu,%.2f,%.2f\n", idx, avg_total, avg_rej);
            fflush(stdout);
        }
    }
    fclose(fp);

    if (runs == 0) {
        printf("A,0,0,0,0,0\n");
        return 0;
    }

    /* A: runs,total_us,setup_us,partial_P_us,partial_A_us,fullkey_P_us,
     *    fullkey_A_us,authorize_us,agentsign_us,verify_us,
     *    rej_P,rej_A,rej_authorize,rej_sign,
     *    rej_total,failed_runs */
    printf("A,%llu,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%llu\n",
           runs,
           sum_total / runs,
           sum_us[0] / runs,
           sum_us[1] / runs,
           sum_us[2] / runs,
           sum_us[3] / runs,
           sum_us[4] / runs,
           sum_us[5] / runs,
           sum_us[6] / runs,
           sum_us[7] / runs,
           sum_rej[0] / runs,
           sum_rej[1] / runs,
           sum_rej[2] / runs,
           sum_rej[3] / runs,
           sum_rejcnt / runs,
           nfail);
    return 0;
}

/*---------------------------------------------------------------
 * dump 模式: 把合格随机数索引转换成实际 32 字节随机数值(hex)导出
 *   bench dump <idxlist_file> > out.txt
 *   每行: <idx> <64位hex随机数>
 *---------------------------------------------------------------*/
static int do_dump(const char *listfile)
{
    FILE *fp = fopen(listfile, "r");
    uint8_t master[SEEDBYTES];
    unsigned long long idx;
    char line[256];
    unsigned int i;

    if (!fp) {
        fprintf(stderr, "cannot open %s\n", listfile);
        return 1;
    }

    printf("# dump list=%s level=%d\n", listfile, DILITHIUM_MODE);
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r')
            continue;
        idx = strtoull(line, NULL, 10);
        master_from_index(master, idx);
        printf("%llu,", idx);
        for (i = 0; i < SEEDBYTES; ++i)
            printf("%02X", master[i]);
        printf("\n");
    }
    fclose(fp);
    return 0;
}

/*---------------------------------------------------------------
 * scanlist 模式: 仅输出合格随机数索引(每行一个数字), 无其他行, 减少 I/O
 *   bench scanlist <start> <count> <thresh>
 *---------------------------------------------------------------*/
static int do_scanlist(unsigned long long start, unsigned long long count,
                       unsigned long long thresh)
{
    unsigned long long i;
    unsigned long long scanned = 0, nfail = 0, nqual = 0;
    uint8_t master[SEEDBYTES];
    stat_t st;

    for (i = 0; i < count; ++i) {
        unsigned long long idx = start + i;
        master_from_index(master, idx);
        run_once(master, &st);
        scanned++;
        if (st.fail)
            nfail++;
        if (st.fail == 0 && st.rej[0] <= thresh && st.rej[1] <= thresh &&
            st.rej[2] <= thresh && st.rej[3] <= thresh) {
            nqual++;
            printf("%llu\n", idx);
        }
    }
    fprintf(stderr, "S,%llu,%llu,%llu\n", scanned, nfail, nqual);
    return 0;
}

/* 对已有候选索引做二次筛选: bench scanfile <idxlist_file> <thresh> */
static int do_scanfile(const char *listfile, unsigned long long thresh)
{
    FILE *fp = fopen(listfile, "r");
    unsigned long long idx;
    unsigned long long scanned = 0, nfail = 0, nqual = 0;
    uint8_t master[SEEDBYTES];
    stat_t st;
    char line[256];

    if (!fp) {
        fprintf(stderr, "cannot open %s\n", listfile);
        return 1;
    }
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r')
            continue;
        idx = strtoull(line, NULL, 10);
        master_from_index(master, idx);
        run_once(master, &st);
        scanned++;
        if (st.fail)
            nfail++;
        if (st.fail == 0 && st.rej[0] <= thresh && st.rej[1] <= thresh &&
            st.rej[2] <= thresh && st.rej[3] <= thresh) {
            nqual++;
            printf("%llu\n", idx);
        }
    }
    fclose(fp);
    fprintf(stderr, "S,%llu,%llu,%llu\n", scanned, nfail, nqual);
    return 0;
}

/* 输出已有候选索引的详细拒绝采样统计 */
static int do_auditfile(const char *listfile, unsigned long long thresh)
{
    FILE *fp = fopen(listfile, "r");
    unsigned long long idx;
    uint8_t master[SEEDBYTES];
    stat_t st;
    char line[256];

    if (!fp) {
        fprintf(stderr, "cannot open %s\n", listfile);
        return 1;
    }
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r')
            continue;
        idx = strtoull(line, NULL, 10);
        master_from_index(master, idx);
        run_once(master, &st);
        if (st.fail == 0 && st.rej[0] <= thresh && st.rej[1] <= thresh &&
            st.rej[2] <= thresh && st.rej[3] <= thresh)
            printf("Q,%llu,%u,%u,%u,%u,%llu,%d\n", idx, st.rej[0], st.rej[1],
                   st.rej[2], st.rej[3], st.total, st.fail);
        else
            printf("F,%llu,%u,%u,%u,%u,%llu,%d\n", idx, st.rej[0], st.rej[1],
                   st.rej[2], st.rej[3], st.total, st.fail);
    }
    fclose(fp);
    return 0;
}

/*---------------------------------------------------------------
 * main
 *---------------------------------------------------------------*/
int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr,
                "usage:\n"
                "  bench scan <start> <count> <thresh> [verbose]\n"
                "  bench scanfile <idxlist_file> <thresh>\n"
                "  bench auditfile <idxlist_file> <thresh>\n"
                "  bench time <idxlist_file> <reps>\n");
        return 2;
    }

    if (strcmp(argv[1], "scan") == 0) {
        unsigned long long start = 0, count = 0, thresh = 5;
        int verbose = 0;
        if (argc >= 5) {
            start  = strtoull(argv[2], NULL, 10);
            count  = strtoull(argv[3], NULL, 10);
            thresh = strtoull(argv[4], NULL, 10);
        }
        if (argc >= 6)
            verbose = atoi(argv[5]);
        return do_scan(start, count, thresh, verbose);
    }

    if (strcmp(argv[1], "time") == 0) {
        int reps = (argc >= 4) ? atoi(argv[3]) : 10;
        return do_time(argv[2], reps);
    }

    if (strcmp(argv[1], "dump") == 0)
        return do_dump(argv[2]);

    if (strcmp(argv[1], "scanlist") == 0) {
        unsigned long long start = strtoull(argv[2], NULL, 10);
        unsigned long long count = strtoull(argv[3], NULL, 10);
        unsigned long long thresh = strtoull(argv[4], NULL, 10);
        return do_scanlist(start, count, thresh);
    }

    if (strcmp(argv[1], "scanfile") == 0)
        return do_scanfile(argv[2], strtoull(argv[3], NULL, 10));

    if (strcmp(argv[1], "auditfile") == 0)
        return do_auditfile(argv[2], strtoull(argv[3], NULL, 10));

    fprintf(stderr, "unknown mode: %s\n", argv[1]);
    return 2;
}
