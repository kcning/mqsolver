/* mqsolver.c: implementation of structure MQSolver
 */

extern "C" {

#include "mqsolver.h"
#include "mqfix.h"

/* according to number of equations and number of variables, dynamically
 * allocate array to store the equations. Note that this function must be
 * called after header of the challenge has be parsed. I.e. after
 * parse_cha_header returns false.
 */
static void
mqs_alloc_sys(MQSolver* mqs) {
    mqs->orig_sys = SMALLOC(bool*, mqs->eq_num);

    uint64_t i;
    for(i = 0; i < mqs->eq_num; i++) {
        mqs->orig_sys[i] = SMALLOC(bool, mqs->xvar_num);
    }
}

/* for parsing challenge file */
const char* CHA_GF_LINE = "Galois Field";
const char* CHA_VAR_LINE = "Number of variables";
const char* CHA_EQ_LINE = "Number of polynomials";
const char* CHA_SEED_LINE = "Seed";
//const char* CHA_ORD_LINE = "Order";
const char* CHA_EQ_START = "*********";
const size_t MAX_PRE_LEN = 128;

/* testing if pre is a prefix of the string */
static inline bool
check_prefix(const char* pre, const char* str) {
    return !strncmp(pre, str, strnlen(pre, MAX_PRE_LEN));
}

/* parse the header of challenge file, return true is still in header.
 * return false otherwise.
 */
static bool
parse_cha_header(MQSolver* mqs, const char* str) {
    bool verbose = mqs->opts->verbose;
    if(check_prefix(CHA_EQ_START, str)) {
        if(verbose) {
            PRINTF_STAMP("\t\treading equations...\n");
        }
        return false;
    }
    
    if(check_prefix(CHA_VAR_LINE, str)) {
        if(1 != sscanf(str, "%*s %*s %*s %*s : %" PRIu64, &(mqs->var_num))) {
            EXIT_WITH_MSG("[!] cannot parse number of unknowns: %s\n", str);
        }

        if(verbose) {
            PRINTF_STAMP("\t\tnumber of variables: %" PRIu64 "\n", mqs->var_num);
        }

    } else if(check_prefix(CHA_EQ_LINE, str)) {
        if(1 != sscanf(str, "%*s %*s %*s %*s : %" PRIu64, &(mqs->eq_num))) {
            EXIT_WITH_MSG("[!] cannot parse number of equations: %s\n", str);
        }

        if(verbose) {
            PRINTF_STAMP("\t\tnumber of equations: %" PRIu64 "\n", mqs->eq_num);
        }

    } else if(check_prefix(CHA_SEED_LINE, str)) {
        if(1 != sscanf(str, "%*s : %d", &(mqs->seed))) {
            EXIT_WITH_MSG("[!] unable to seed: %s\n", str);
        }

        if(verbose) {
            PRINTF_STAMP("\t\tseed: %d\n", mqs->seed);
        }

    } else if(check_prefix(CHA_GF_LINE, str)) {
        int prime = 0;
        if( (1 != sscanf(str, "%*s %*s : GF(%d)", &prime)) || prime != 2) {
            EXIT_WITH_MSG("[!] unable to process GF(%d)\n", prime);
        }

        if(verbose) {
            PRINTF_STAMP("\t\tfield: GF(%d)\n", prime);
        }
    }

    return true;
}

/* reduce the system to a simpler form. In particular, reduce all
 * x_i^2 to x_i
 */
static void
reduce_sys(MQSolver* mqs) {
    uint64_t eq_idx, var_idx, i, sqr_term_idx;
    for(eq_idx = 0; eq_idx < mqs->eq_num; ++eq_idx) {
        for(var_idx = 0; var_idx < mqs->var_num; ++var_idx) {
            for(i = 0, sqr_term_idx = 0; i < var_idx; i++) {
                sqr_term_idx += i + 2;
            }
            // add the coefficient of x_{var_idx}^2 to x_{var_idx}
            mqs->orig_sys[eq_idx][mqs->xvar_num-2-(mqs->var_num-1-var_idx)] ^=
                mqs->orig_sys[eq_idx][sqr_term_idx];
            mqs->orig_sys[eq_idx][sqr_term_idx] = 0;
        }
    }
}

/* parse the system of challenge file. Note this will destroy the string */
static void
parse_cha_eqs(MQSolver* mqs, char* str, const uint64_t eq_idx) {
    char* ptr = NULL;
    
    uint64_t i = 0;
    ptr = strtok(str, " ;");
    while(NULL != ptr) {
        mqs->orig_sys[eq_idx][i++] = atoi(ptr);
        ptr = strtok(NULL, " ;\n");
    }

    assert(i == mqs->xvar_num);
}

/* function: mqs_readcha
 * usage: read MQ challenge file
 * arguments: a pointer to structure MQSolver
 * return: void
 */
static void
mqs_readcha(MQSolver* mqs) {
    if(NULL == mqs->opts->cha_file) {
        EXIT_WITH_MSG("[!] no challenge file is given\n");
    }
    PRINTF_STAMP("[+] reading challenge from: %s\n", mqs->opts->cha_file);

    FILE* fp = NULL;
    if(NULL == (fp = fopen(mqs->opts->cha_file, "r"))) {
        EXIT_WITH_MSG("[!] cannot read from challenge file: %s\n",
            mqs->opts->cha_file);
    }

    // NOTE: expand the buffer if needed
    const size_t buf_size = 0x1 << 20; // 1MB per line
    char* buf = SMALLOC(char, buf_size);
    bool before_eq = true;
    uint64_t eq_idx = 0;
    while(NULL != fgets(buf, buf_size, fp)) {
        if(true == before_eq) {
            before_eq = parse_cha_header(mqs, buf);

            if(false == before_eq) {
                mqs->xvar_num = sum_binom(mqs->var_num, 2) + mqs->var_num;
                mqs_alloc_sys(mqs);
            }
        } else {
            parse_cha_eqs(mqs, buf, eq_idx++);
        }
    }

    reduce_sys(mqs);

    assert(eq_idx == mqs->eq_num);
    fclose(fp);
    free(buf);
}

/* function: print_sys_info
 * usage: read and print system info
 * arguments: void
 * return: void
 */
static void
print_sys_info() {
    PRINTF_STAMP("[+] system info:\n");
    PRINTF_STAMP("\tCPU frequency: %f\n", get_cpu_freq());

    size_t buf_size = 1024;
    char buf[buf_size];
    if(0 != get_proc_status(buf,  buf_size)) {
        PRINTF_ERR_STAMP("[!] error reading proc status\n");
    } else {
        PRINTF_STAMP("\tproc status:\n");

        char* ptr = NULL;
        ptr = strtok(buf, "\n");
        while(ptr != NULL) {
            PRINTF_STAMP("\t\t%s\n", ptr);
            ptr = strtok(NULL, "\n");
        }
    }
}

/* function: mqs_init
 * usage: init the MQSolver structure
 *      1) parse options
 *      2) read system info
 *      3) init random number generator
 *      4) read the challenge file passed via Options
 * arguments:
 *      1) mqs: a pointer to structure MQSolver
 *      2) argc
 *      3) argv
 * return: void
 */
void
mqs_init(MQSolver* mqs, int argc, char* argv[]) {
    PRINTF_STAMP("[+] record initial timestamp...\n");  /* set __init_stamp */
    PRINTF_DEBUG("[+] debug mode on!\n");

    Options* opts = SMALLOC(Options, 1);
    options_init(opts);
    options_parse(opts, argc, argv);

    if(opts->verbose) print_sys_info();

    /* init random number generator */
    srand(opts->seed);
    PRINTF_STAMP("[+] random seed: 0x%x\n", opts->seed);

    mqs->opts = opts;
    mqs->orig_sys = NULL;
    mqs->var_num = 0;
    mqs->eq_num = 0;
    mqs->xvar_num = 0;
    mqs->seed = 0;

    mqs_readcha(mqs);
}

/* function: mqs_free
 * usage: free the MQSolver structure and its fields
 * arguments: an pointer to structure MQSolver
 * return: void
 */
void
mqs_free(MQSolver* mqs){
    if(NULL != mqs->opts) {
        options_free(mqs->opts);
        SFREE(mqs->opts);
    }

    if(NULL != mqs->orig_sys) {
        uint64_t i;
        for(i = 0; i < mqs->eq_num; i++) {
            SFREE(mqs->orig_sys[i]);
        }
        SFREE(mqs->orig_sys);
    }
}

/* ===========================================================================
 *                              algorithm implementation
 * ===========================================================================
 */

/* compute the partial derivatives of each equation in the system w.r.
 * each variable
 */
static void
find_partial_derivs(MQSolver* mqs, bool* derivs, uint64_t eq_num,
                   uint64_t term_num, uint64_t var_num) {
    uint64_t eq_idx, var_idx;
    for(eq_idx = 0; eq_idx < eq_num; eq_idx++) {
        for(var_idx = 0; var_idx < var_num; var_idx++) {
            diff_eq(mqs->orig_sys[eq_idx], term_num, var_num, var_idx,
                    derivs + eq_idx * var_num * (var_num+1) +
                        var_idx * (var_num+1));
        }
    }
}

/* function: fast_exhaustive
 * usage: use bitsliced fast exhaustive search to solve MQ system with less
 *      than 64 variables and 64 equations.
 * arguments:
 *      1) mqs: a pointer to MQSolver, which contains the system to solve
 *      2) solution: a container for the solution. The size of the container
 *              is defined by mqs->var_num
 * return: true if a solution is found, otherwise false
 */
static bool
fast_exhaustive(MQSolver* mqs, bool* solution) {
    const uint64_t var_num = mqs->var_num;
    const uint64_t eq_num = mqs->eq_num;

    if(var_num > 64 || eq_num > 64) {
        EXIT_WITH_MSG("[!] too many variables and/or equations, "
                      "use crossbred instead\n");
    }

    uint64_t eq_idx, var_idx, i, term;
    bool derivs[eq_num][var_num][var_num+1];
    find_partial_derivs(mqs, (bool*)derivs, eq_num,
                        mqs->xvar_num, var_num);

    // the second order partial dervatives of the system
    // pdiff2[j][i]'s eq_idx-th bit, starting from LSB, holds the partial
    // dervatives of equation eq_idx against first x_i then x_j,
    // which is a constant.
    uint64_t pdiff2[var_num][var_num];
    memset(pdiff2, 0x0, sizeof(uint64_t) * var_num * var_num);
    for(var_idx = 0; var_idx < var_num; ++var_idx) {
        for(i = 0; i < var_num; ++i) {
            for(eq_idx = 0; eq_idx < eq_num; ++eq_idx) {
                term = derivs[eq_idx][var_idx][i];
                assert(0x1UL == term || 0x0UL == term);
                pdiff2[i][var_idx] |= term << eq_idx;
            }
        }
    }

    uint64_t pdiff_eval[var_num];
    memset(pdiff_eval, 0x0, sizeof(uint64_t) * var_num);
    for(var_idx = 0; var_idx < var_num; ++var_idx) {
        for(eq_idx = 0; eq_idx < eq_num; ++eq_idx) {
            if(0 == var_idx) {
                term = derivs[eq_idx][0][var_num];
            } else {
                term = derivs[eq_idx][var_idx][var_num] ^
                    derivs[eq_idx][var_idx][var_idx-1];
            }
            assert(0x1UL == term || 0x0UL == term);
            pdiff_eval[var_idx] |= term << eq_idx;
        }
    }

    PRINTF_STAMP("\t\tbrute forcing...\n");

    // func_eval's eq_idx-th bit (from LSB) holds the result of evaluating
    // equation eq_idx at the current solution
    uint64_t func_eval = 0x0UL;
    // evaluate all functions at initial solution, zero vector
    for(eq_idx = 0; eq_idx < eq_num; ++eq_idx) {
        term = mqs->orig_sys[eq_idx][mqs->xvar_num-1];
        assert(0x1UL == term || 0x0UL == term);
        func_eval |= term << eq_idx;
    }

    uint64_t count = 0;
    const uint64_t bound = (0x1UL << var_num) - 1;
    while(func_eval && count < bound) {
        ++count;
        uint64_t fp_idx = __builtin_ctzll(count);

        if( count & (count-1) ) {
            uint64_t pre_fp_idx = __builtin_ctzll(count ^ (0x1UL << fp_idx));
            pdiff_eval[fp_idx] ^= pdiff2[fp_idx][pre_fp_idx];
        }

        func_eval ^= pdiff_eval[fp_idx];
    }

    // fill in the solution
    for(var_idx = 0; var_idx < var_num; ++var_idx) {
        solution[var_idx] = ((count ^ count >> 1) >> var_idx) & 0x1UL;
    }
    return !func_eval;
}

/* A subroutine of crossbred: simplify an MQ system by reducing square terms */
static inline void
simp_sys(bool* const dst, MQSolver* mqs) {

    const uint64_t term_num = mqs->xvar_num;
    const uint64_t eq_len = mqs->xvar_num - mqs->var_num;

    for(uint64_t i = 0; i < mqs->eq_num; ++i) {

        // constant and linear terms
        for(uint64_t j = 0; j < mqs->var_num+1; ++j) {
            dst[i * eq_len + eq_len-1-j] = mqs->orig_sys[i][term_num-1-j];
        }

        // degree 2 monomials
        uint64_t sqr_term_idx = term_num-1-mqs->var_num-1;
        uint64_t deg2_idx = eq_len-1-mqs->var_num-1;

        for(uint64_t j = 0; j < mqs->var_num; ++j) {
            // square term: x^2 => x
            dst[i * eq_len + eq_len-1-1-j] ^= mqs->orig_sys[i][sqr_term_idx--];

            for(uint64_t k = 0; k < mqs->var_num-1-j; ++k) {
                dst[i * eq_len + deg2_idx--] = mqs->orig_sys[i][sqr_term_idx--];
            }
        }
    }
}

/* check if a solution works for a linear system. The linear system must be
 * free of square terms.
 */
static inline bool
verify_sol(bool* solution, bool* sys, const uint64_t eq_num,
           const uint64_t var_num, const uint64_t term_num,
           const uint64_t start) {
    for(uint64_t i = start; i < eq_num; ++i) { // for each equation
        bool res = sys[i * term_num + term_num-1]; // constant term

        // linear terms
        bool* ptr = sys + i * term_num + term_num-1-1;
        for(uint64_t j = 0; j < var_num; ++j) {
            res ^= solution[var_num-1-j] & *ptr;
            --ptr;
        }

        // quadratic terms
        for(uint64_t j = 0; j < var_num; ++j) {
            if(false == solution[var_num-1-j]) {
                // the var is zero
                ptr -= var_num-1-j;
                continue;
            }

            for(uint64_t k = 1; k < var_num-j; ++k) {
                res ^= *ptr & solution[var_num-1-j-k];
                --ptr;
            }
        }

        if(res) { // the equation is evaluated to 1
            return false;
        }
    }

    return true;
}

/* A subroutine of crossbred, copy the first keq_num eqs in the original MQ
 * system into a column-wise storage.
 */
static void
split_sys(uint32_t* const __restrict__ dst,
          const bool* const __restrict__ sys, const uint64_t eq_num,
          const uint64_t term_num, const uint64_t keq_num) {
    assert(32 >= keq_num);
    assert(eq_num >= keq_num);
    memset(dst, 0x0, sizeof(uint32_t) * term_num);

    for(uint64_t i = 0; i < keq_num; ++i) {
        for(uint64_t j = 0; j < term_num; ++j) {
            if(sys[i * term_num + j]) {
                dst[j] |= 0x1 << i;
            }
        }
    }
}

/* compute the sparsity of the extracted sub-system. */
static double
sparsity(RMac* sys, const uint64_t sub_eq_num, const uint64_t ltnum) {

    const uint64_t start = sys->term_num - ltnum;
    double sparsity = 0.0;

    for(uint64_t i = 0; i < sub_eq_num; ++i) {
        for(uint64_t j = 0; j < ltnum; ++j) {
            if( !rmac_at(sys,  sys->eq_num - sub_eq_num + i, start + j) ) {
                sparsity += 1;
            }
        }
    }

    return sparsity / (sub_eq_num * ltnum);
}

/* A subroutine of crossbred, verify the sub-system row candidates */
static uint64_t
verify_subsys(uint64_t* const valid_eq_indices, RMac* sys,
              const uint64_t pvt_num, const uint64_t keq_num) {
    uint64_t c = 0;

    for(uint64_t i = pvt_num; i < pvt_num + keq_num; ++i) {

        if( !rmac_zrow(sys, i) ) {
            valid_eq_indices[c++] = i;
        }
    }

    return c;
}

/* A subroutine of crossbred, transpose but still keep non-lienar terms, which
 * are guaranteed to be zero in the sub-system. Only at most 32 eqs in the
 * sub-system are kept. Keeping the non-linear terms is helpful for evaluating
 * the sub-system and the partial derivs of it.
 */
static void
pick_subsys(uint32_t* subsys, RMac* sys, const uint64_t mac_term_num,
            const uint64_t keq_num, const uint64_t veq_num,
            const uint64_t* const __restrict__ valid_eq_indices,
            const uint64_t* const __restrict__ l_indices,
            const uint64_t lt_num) {
    // reset subsys
    memset(subsys, 0x0, sizeof(uint32_t) * mac_term_num);
    const uint64_t start_idx = sys->term_num - lt_num;

    // Floyd's random sampling algorithm
    bool* is_used = (bool*) calloc(veq_num, sizeof(bool)); // flags
    if(NULL == is_used) {
        EXIT_WITH_MSG("[!] insufficient memory\n");
    }

    uint64_t cnum = 0; // number of chosen indices
    for(uint64_t i = veq_num - keq_num; i < veq_num && cnum < keq_num; ++i) {
        uint64_t r = rand() % (i + 1);
        if(is_used[r]) {
            r = i;
        }
        
        assert(!is_used[r]);
        is_used[r] = true;
        ++cnum;

        const uint64_t eq_idx = valid_eq_indices[r];
        for(uint64_t j = 0; j < lt_num; ++j) {
            uint32_t b = rmac_at(sys, eq_idx, start_idx + j);
            subsys[l_indices[j]] |= b << i; // start from LSB
        }
    }

    assert(cnum == keq_num);
    free(is_used);
}

/* extract solution from the uint32_t's returned by GPU thread */
static inline void
fill_sol(bool* const sol, const uint32_t knum, const uint32_t tfnum,
         const uint32_t kfnum, const uint32_t sfnum,
         const uint32_t lsys_sol, const uint32_t tf_sol,
         const uint32_t kf_sol, const uint32_t sfval) {

    const uint32_t gc = tf_sol ^ (tf_sol >> 1);
    const uint32_t front_num = knum + tfnum;
    const uint32_t vnum = front_num + kfnum;

    for(uint64_t i = 0; i < knum; ++i) {
        sol[i] = (lsys_sol >> i) & 0x1U;
    }
    
    for(uint64_t i = knum; i < front_num; ++i) {
        sol[i] = (gc >> (front_num-1-i)) & 0x1U;
    }

    for(uint64_t i = front_num; i < vnum; ++i) {
        sol[i] = (kf_sol >> (vnum-1-i)) & 0x1U;
    }

    for(uint64_t i = 0; i < sfnum; ++i) {
        sol[vnum + i] = (sfval >> (sfnum-1-i)) & 0x1U;
    }
}

/* dump the Macaulay matrix into a file in gnuplot format */
static void
dump_mac(Macaulay* mac, const char* fname, const uint64_t eq_num) {
    FILE* f = fopen(fname, "w");

    if(NULL == f) {
        PRINTF_ERR("[!] could not open dump file: %s\n", fname);
        return;
    }

    for(uint64_t i = 0; i < eq_num; ++i) {
        for(uint64_t j = 0; j < mac->term_num; ++j) {
            if(true == mac_at(mac, i, j)) {
                fprintf(f, "%" PRIu64 " %" PRIu64 "\n", j, eq_num-1-i);
            }
        }
    }
    fclose(f);
}

/* dump the reduced Macaulay matrix into a file in gnuplot format */
static void
dump_rmac(RMac* sys, const char* fname) {
    FILE* f = fopen(fname, "w");

    if(NULL == f) {
        PRINTF_ERR("[!] could not open dump file: %s\n", fname);
        return;
    }

    for(uint64_t i = 0; i < sys->eq_num; ++i) {
        for(uint64_t j = 0; j < sys->term_num; ++j) {
            if(true == rmac_at(sys, i, j)) {
                fprintf(f, "%" PRIu64 " %" PRIu64 "\n", j, sys->eq_num-1-i);
            }
        }
    }
    fclose(f);
}

/* subroutine of crossbred, fix variables in the reduced MQ system */
static inline void
fix_mq(bool* dst, const bool* const src, MQFix* mqf, const uint64_t tvnum,
       const uint64_t eq_num, const uint64_t idx, const uint64_t dst_eq_len,
       const uint64_t src_eq_len) {

    const uint64_t vnum = tvnum - mqf->fnum;
    const uint64_t deg2_tnum = binom(vnum, 2);
    const uint64_t src_deg2_tnum = binom(tvnum, 2);

    for(uint64_t i = 0; i < eq_num; ++i) {
        // copy degree 2 monomials that stay the same
        memcpy(dst + i * dst_eq_len, src + i * src_eq_len, sizeof(bool) * deg2_tnum);

        // copy linear terms that stay the same
        memcpy(dst + i * dst_eq_len + deg2_tnum,
               src + i * src_eq_len + src_deg2_tnum, sizeof(bool) * vnum);

        // copy constant term
        dst[i * dst_eq_len + dst_eq_len-1] = src[i * src_eq_len + src_eq_len-1];

        const bool* ptr = src + i * src_eq_len + deg2_tnum;
        // reduce degree 2 monomial
        for(uint64_t j = 0; j < mqf->fnum; ++j) {
            if( !mqfix_lvar(mqf, idx, j) ) {
                ptr += vnum + j;
                continue;
            }

            // reduce into linear term
            for(uint64_t k = 0; k < vnum; ++k) {
                dst[i * dst_eq_len + deg2_tnum + k] ^= ptr[k];
            }

            // reduce into constant
            for(uint64_t k = vnum; k < vnum + j; ++k) {
                if( mqfix_lvar(mqf, idx, k-vnum) ) {
                    dst[i * dst_eq_len + dst_eq_len-1] ^= ptr[k];
                }
            }

            ptr += vnum + j;
        }
        
        ptr += vnum;
        // reduce linear terms into constant term
        for(uint64_t j = 0; j < mqf->fnum; ++j) {
            if( !mqfix_lvar(mqf, idx, j) ) {
                continue;
            }

            dst[i * dst_eq_len + dst_eq_len-1] ^= ptr[j];
        }
    }
}

/* arguments for bf_subsys */
typedef struct {
    uint64_t mac_term_num;
    uint64_t mq_tnum;
    uint64_t var_num;
    uint64_t thread_num;
    uint64_t sfnum;
    uint64_t deg;
    uint64_t kvar_num;
    uint64_t idx;

    uint32_t mempw;
    uint32_t warp_num;
    uint32_t mailbox_size;
    uint32_t max_scnum;
    uint32_t cmem_size;
    uint32_t fsub_tnum;
    uint32_t svnum;
    uint32_t tfnum;
    uint32_t kfnum;
    int dev_id;

    dim3 gdim;
    dim3 bdim;

    uint32_t* host_32eqs;
    uint32_t* subsys;
    bool* solution;
    bool* sys;
    bool* solvable;
    MQFix* mqf;
    MQSolver* mqs;
    
    pthread_mutex_t* bf_done_lock;
    pthread_cond_t* bf_done_cond;

    bool verbose;
    bool* bf_done;
    
} bf_subsys_arg;

/* subroutine of crossbred, fix variables in the sub-system with all possible
 * combination and bruteforce it with GPU. This function is created to pipeline
 * bruteforcing on GPU and Reduced Macaulay matrix computation on CPU and is not
 * thread-safe! There should be no more than one thread executing this function.
 */
static void
bf_subsys(void* dummy) {
    bf_subsys_arg* arg = (bf_subsys_arg*) dummy;

    // reset the device for this host thread
    CUDA_CHECK(cudaSetDevice(arg->dev_id));

    uint32_t* dev_subsys = NULL;
    CUDA_CHECK(cudaMalloc(&dev_subsys, sizeof(uint32_t) * arg->mac_term_num));

    uint32_t* dev_32eqs = NULL;
    CUDA_CHECK(cudaMalloc(&dev_32eqs, sizeof(uint32_t) * arg->mq_tnum));

    bool* tmp = NULL;
    CUDA_CHECK(cudaMalloc(&tmp, sizeof(bool) * arg->var_num * arg->thread_num));

    // copy the filter (32 eqs) to GPU
    CUDA_CHECK(cudaMemcpy(dev_32eqs, arg->host_32eqs,
                          sizeof(uint32_t) * arg->mq_tnum, cudaMemcpyHostToDevice));

    // allocate storage for lsys and graycode enumeration
    uint32_t* gc_data = NULL;
    CUDA_CHECK(cudaMalloc(&gc_data, arg->mempw * arg->warp_num));

    // allocate storage for solution candidates (mailbox)
    uint32_t* gc_mailbox = NULL;
    CUDA_CHECK(cudaMalloc(&gc_mailbox, arg->mailbox_size));

    uint32_t* host_mailbox = NULL;
    CUDA_CHECK(cudaMallocHost(&host_mailbox, arg->mailbox_size));

    // prepare container for sub-system after fixing variables
    uint32_t* fix_subsys = SMALLOC(uint32_t, arg->fsub_tnum);

    // prepare container for last order parial derivs
    uint32_t* cmem = (uint32_t*) malloc(arg->cmem_size);
    if(NULL == cmem) {
        EXIT_WITH_MSG("[!] insufficient memory\n");
    }

    /* ====================================================================
     *                       for each fixed sub-system
     * ====================================================================
     */
    assert( false == *(arg->solvable) );
    const uint64_t iter_num = 0x1UL << arg->sfnum;
    for(uint64_t sfval = 0; !*(arg->solvable) && sfval < iter_num; ++sfval) {

        // reset
        CUDA_CHECK(cudaMemset(gc_data, 0x0, arg->mempw * arg->warp_num));
        CUDA_CHECK(cudaMemset(gc_mailbox, 0x0, arg->mailbox_size));

        if(arg->sfnum) {
            PRINTF_STAMP("[+] fixing %" PRIu64 " variables in the sub-system\n",
                         arg->sfnum);
        }
        fix_sub(fix_subsys, arg->subsys, arg->deg, arg->var_num,
                arg->sfnum, sfval);

        // copy the sub-system to GPU
        CUDA_CHECK(cudaMemcpy(dev_subsys, fix_subsys,
                              sizeof(uint32_t) * arg->fsub_tnum,
                              cudaMemcpyHostToDevice));

        PRINTF_STAMP("[+] initializing graycode data structure...\n");
        gc_init<<<arg->gdim, arg->bdim>>>(gc_data, dev_subsys, arg->deg, arg->svnum,
                                          arg->tfnum, arg->kfnum, arg->kvar_num);

        CUDA_CHECK(cudaPeekAtLastError());
        
        // init constant memory on the device
        gc_cmeminit(cmem, fix_subsys, arg->tfnum, arg->kfnum, arg->kvar_num, arg->deg);
        if(arg->verbose) {
            PRINTF_STAMP("\t\tcopying last order paritial derivatives (%"
                         PRIu32" bytes) to GPU...\n", arg->cmem_size);
        }
        CUDA_CHECK(cudaMemcpyToSymbol(ldderivs, cmem, arg->cmem_size));

        CUDA_CHECK(cudaFuncSetCacheConfig(gc_bfsearch, cudaFuncCachePreferL1));
        CUDA_CHECK(cudaFuncSetSharedMemConfig(gc_bfsearch,
                                              cudaSharedMemBankSizeFourByte));

        PRINTF_STAMP("[+] bruteforcing with GPU...\n");
        gc_bfsearch<<<arg->gdim, arg->bdim>>>(gc_data, arg->deg, arg->tfnum,
                                              arg->kfnum, arg->sfnum, sfval,
                                              gc_mailbox, arg->max_scnum,
                                              dev_32eqs, tmp);

        CUDA_CHECK(cudaPeekAtLastError());

        // TODO: sync on stream
        CUDA_CHECK(cudaMemcpy(host_mailbox, gc_mailbox, arg->mailbox_size,
                              cudaMemcpyDeviceToHost));

        // block the thread until GPU is done
        CUDA_CHECK(cudaDeviceSynchronize());

        PRINTF_STAMP("[+] checking solution candidates...\n");
        uint32_t sol_num = *gc_mbcount(host_mailbox);
        if(arg->verbose) {
            PRINTF_STAMP("\t\tnumber of solution candidates: %" PRIu32 "\n", sol_num);
            PRINTF_STAMP("\t\tnumber of filtered false positives: %"
                         PRIu32"\n", *gc_mbsnum(host_mailbox) - sol_num);
        }

#ifdef GC_DEPC_LSYS
       PRINTF_STAMP("\t\tnumber of dependent linear systems: %" PRIu32 "\n",
                    *gc_mbdepc(host_mailbox));
#endif

        if(arg->max_scnum <= sol_num) {
            PRINTF_ERR_STAMP("[!] the number of solution candidates is more"
                             " than expected\n");
            PRINTF_ERR_STAMP("\t\texpected maximum: %" PRIu32 "\n", arg->max_scnum);
            PRINTF_ERR_STAMP("\t\tsome candidates were therefore dropped\n");
            sol_num = arg->max_scnum;
        }

        for(uint32_t i = 0; !*(arg->solvable) && i < sol_num; ++i) {
            uint32_t* const slot = gc_mbslot(host_mailbox, i);
            fill_sol(arg->solution, arg->kvar_num, arg->tfnum, arg->kfnum,
                     arg->sfnum, *gc_mblsys(slot), *gc_mbtsol(slot),
                     *gc_mbksol(slot), sfval);

            *(arg->solvable) = verify_sol(arg->solution, arg->sys, arg->mqs->eq_num,
                                          arg->var_num, arg->mq_tnum, 32);

            if(*(arg->solvable) && arg->mqf) {
                bool prep_vals[arg->mqf->fnum];
                mqfix_ltob(arg->mqf, arg->idx, prep_vals);

                memcpy(arg->solution + arg->var_num, prep_vals,
                       sizeof(bool) * arg->mqf->fnum);
            }
        }
    } // for each fixed sub-system

    // signal the main thread
    pthread_mutex_lock(arg->bf_done_lock);
    *(arg->bf_done) = true;
    pthread_cond_signal(arg->bf_done_cond);
    pthread_mutex_unlock(arg->bf_done_lock);

    CUDA_CHECK(cudaFreeHost(host_mailbox));
    CUDA_CHECK(cudaFree(dev_subsys));
    CUDA_CHECK(cudaFree(gc_data));
    CUDA_CHECK(cudaFree(gc_mailbox));
    CUDA_CHECK(cudaFree(dev_32eqs));
    CUDA_CHECK(cudaFree(tmp));
    free(fix_subsys);
    free(cmem);
    free(arg->host_32eqs);
    free(arg->subsys);
    free(arg->sys);
    free(arg);
}

const char* dsubsys = "subsys";
const char* dheq = "host_32eqs";
const char* dsys = "sys";
const char* dlim = "==";

/* subroutine of crossbred. Dump the sub-system and relevant info into a file */
static void
dump_subsys(FILE* f, const uint32_t* const subsys, const uint64_t mac_term_num,
            MQFix* mqf, const uint64_t lidx, const uint32_t* const host_32eqs,
            const uint64_t mq_tnum, const bool* const sys, const uint64_t eq_num) {

    // TODO: dump more info to avoid passing parameters again
    fprintf(f, "%s: %" PRIu64 "\n", dsubsys, mac_term_num);
    for(uint64_t i = 0; i < mac_term_num; ++i) {
        uint32_t term = subsys[i];
        fprintf(f, "%" PRIu32 "\n", term);
    }
    
    fprintf(f, "%s:%" PRIu64 "\n", dheq, mq_tnum);
    for(uint64_t i = 0; i < mq_tnum; ++i) {
        uint32_t term = host_32eqs[i];
        fprintf(f, "%" PRIu32 "\n", term);
    }

    fprintf(f, "%s:%" PRIu64 "\n", dsys, eq_num);
    for(uint64_t i = 0; i < eq_num; ++i) {
        for(uint64_t j = 0; j < mq_tnum; ++j) {
            fprintf(f, "%d ", sys[i * mq_tnum + j]);
        }
        fprintf(f, "\n");
    }

    fprintf(f, "%s\n", dlim); // deliminator
}

/* subroutine of crossbred. Parse the output file from dump_subsys to recover the
 * subsystems and relevant info.
 */
static void
parse_subsys(FILE* f, uint32_t* const subsys, uint32_t* const host_32eqs,
             bool* const sys) {
    // NOTE: expand the buffer if needed
    const size_t buf_size = 0x1 << 20; // 1MB per line
    char* buf = SMALLOC(char, buf_size);

    // sub
    uint64_t mac_term_num = 0;
    if( 1 != fscanf(f, "subsys:%" PRIu64 "\n", &mac_term_num) ) {
        EXIT_WITH_MSG("[!] wrong header line for sub-system\n");
    }

    for(uint64_t i = 0; i < mac_term_num; ++i) {
        if ( 1 != fscanf(f, "%" PRIu32 "\n", subsys + i) ) {
            EXIT_WITH_MSG("[!] invalid line for sub-system\n");
        }
    }

    // 32 eqs
    uint64_t mq_tnum = 0;
    if( 1 != fscanf(f, "host_32eqs:%" PRIu64 "\n", &mq_tnum) ) {
        EXIT_WITH_MSG("[!] wrong header line for host equations\n");
    }

    for(uint64_t i = 0; i < mq_tnum; ++i) {
        if ( 1 != fscanf(f, "%" PRIu32 "\n", host_32eqs + i) ) {
            EXIT_WITH_MSG("[!] invalid line for host equations\n");
        }
    }

    // sys
    uint64_t eq_num = 0;
    if( 1 != fscanf(f, "sys:%" PRIu64 "\n", &eq_num) ) {
        EXIT_WITH_MSG("[!] wrong header line for MQ system\n");
    }

    char* digit = NULL;
    for(uint64_t i = 0; i < eq_num; ++i) {
        if(NULL == fgets(buf, buf_size, f)) {
            EXIT_WITH_MSG("[!] failed to read from file\n");
        }

        for(uint64_t j = 0; j < mq_tnum; ++j) {
            digit = strtok(buf, " \n");
            if ( NULL == digit ) {
                EXIT_WITH_MSG("[!] invalid line for system\n");
            }

            sys[i * mq_tnum + j] = atoi(digit);
        }
    }

    if( NULL == fgets(buf, buf_size, f) || strncmp(buf, dlim, strlen(dlim)) ) {
        EXIT_WITH_MSG("[!] missing deliminator\n");
    }

    free(buf);
}

/* function: crossbred
 * usage: use an adaption of Joux and Vitse's crossbred algorithm to solve MQ system
 *
 * arguments:
 *      1) mqs: a pointer to MQSolver, which contains the system to solve
 *      2) solution: a container for the solution. The size of the container
 *              is defined by mqs->var_num
 * return: true if a solution is found, otherwise false
 */
static bool
crossbred(MQSolver* mqs, bool* solution) {
    PRINTF_STAMP("[+] sanity check...\n");

    if(mqs->var_num <= 32 && mqs->eq_num <= 64) {
        PRINTF_ERR_STAMP("[!] less than 32 variables... use fast exhaustive search\n");
        return fast_exhaustive(mqs, solution);
    }

    const uint64_t kvar_num = KNUM;
    //const uint64_t kvar_num = mqs->opts->crossbred_kvar_num;
    if(2 > kvar_num) {
        // see graycode_update in graycode.cu
        // the result is still correct, but it incurs invalid read/write access
        PRINTF_ERR_STAMP("[!] keep at least 2 variables\n");
        return false;
    }

    MQFix* mqf = NULL;
    if(mqs->opts->mq_fix_file) {
        mqf = mqfix_create(mqs->opts->mq_fix_file);
    }

    const bool verbose = mqs->opts->verbose;
    const uint64_t pvnum = mqf ? mqf->fnum : 0;
    const uint64_t sfnum = mqs->opts->crossbred_sub_fnum;
    const uint64_t var_num = mqs->var_num - pvnum;
    const uint64_t mq_tnum = sum_binom(var_num, 2);
    const uint64_t deg = mqs->opts->crossbred_mac_degree;
    const double mq_term_ratio = mqs->opts->crossbred_mq_ratio;
    const uint32_t kfnum = mqs->opts->crossbred_kf_num;
    const uint64_t fvar_num = var_num - kvar_num;
    const uint64_t nlterm_num = mac_nl_num(var_num, deg, fvar_num);
    const uint64_t lterm_num = mac_col_num(var_num, deg) - nlterm_num;

    uint64_t max_indep_rows = mac_row_num(mqs->eq_num, var_num, deg);
    if(4 == deg) {
        max_indep_rows -= binom2(mqs->eq_num) + mqs->eq_num;
    }

    if(pvnum && verbose) {
        PRINTF_STAMP("\t\tnumber of preprocessed variables: %" PRIu64 "\n",
                     pvnum);
        PRINTF_STAMP("\t\tnumber of preprocess configurations: %" PRIu64 "\n",
                     mqf->lnum);
    }

    if(verbose) {
        PRINTF_STAMP("\t\tspecified maximal MQ system density: %f\n", mq_term_ratio);
        PRINTF_STAMP("\t\tnumber of monomials that must be eliminated: %"
                     PRIu64 "\n", nlterm_num);
        PRINTF_STAMP("\t\tnumber of monomials that can be turned into linear "
                     "terms: %" PRIu64 "\n", lterm_num);
        PRINTF_STAMP("\t\tmaximal number of independent rows in Macaulay matrix: %"
                     PRIu64 "\n", max_indep_rows);
    }

    if(max_indep_rows <= nlterm_num) {
        max_indep_rows = 0;
    } else {
        max_indep_rows -= nlterm_num;
    }
    
    if(verbose) {
        PRINTF_STAMP("\t\tmaximal number of equations that can be extracted: %"
                     PRIu64 "\n", max_indep_rows);
    }

    if(max_indep_rows < kvar_num) {
        PRINTF_ERR_STAMP("\t\t[!] impossible to extract enough linear eqs, aborting\n");
        return false;
    }

    if(verbose) PRINTF_STAMP("\t\tcomputing indices of non-linear terms...\n");
    const uint64_t mac_term_num = mac_col_num(var_num, deg);
    const uint64_t mac_eq_num = mac_row_num(mqs->eq_num, var_num, deg);
    uint64_t* nl_indices = SMALLOC(uint64_t, nlterm_num);
    uint64_t* l_indices = SMALLOC(uint64_t, lterm_num);
    uint64_t* idx_map = SMALLOC(uint64_t, mac_term_num);
    uint64_t* clz = SMALLOC(uint64_t, mac_eq_num);
    mac_sep_indices(deg, var_num, nl_indices, l_indices, kvar_num);
    mac_idx_map(idx_map, nl_indices, l_indices, nlterm_num, lterm_num);

    uint64_t keq_suggest = mqs->opts->crossbred_mac_keq_num;
    const uint64_t mac_keq_num = MIN(mac_eq_num - nlterm_num, keq_suggest);
    if(verbose) {
        PRINTF_STAMP("\t\tkeep %" PRIu64 " equation candidates for sub-system\n",
                     mac_keq_num);
    }

    // TODO: make this a parameter
    const double exp_nfr_ratio = 0.35;
    const double exp_rmemsize = rmac_memsize(exp_nfr_ratio * nlterm_num + mac_keq_num,
                                             exp_nfr_ratio * nlterm_num + lterm_num);
    if(verbose) {
        PRINTF_STAMP("\t\texpected maximal ratio of missing pivot rows: %f\n",
                     exp_nfr_ratio);
        PRINTF_STAMP("\t\texpected maximal size of reduced Macaulay matrix: %.2f MB\n",
                     exp_rmemsize / MBFLOAT);
    }

    uint64_t* valid_eq_indices = SMALLOC(uint64_t, mac_keq_num);
    uint32_t* drow = SMALLOC(uint32_t, nlterm_num + mac_keq_num);

    // duplicate the original system into consecutive memory blocks and
    // drop square terms
    bool* const reduc_mq = SMALLOC(bool, mqs->eq_num * (mqs->xvar_num - mqs->var_num));
    simp_sys(reduc_mq, mqs);

    // number variables in the sub-system after fixing
    const uint32_t svnum = var_num - sfnum;
    const uint32_t fsub_tnum = mac_col_num(svnum, deg);

    /* =========================================================================
     *                            gpu configuration
     * =========================================================================
     */
    PRINTF_STAMP("[+] checking GPU devices...\n");
    // check GPU
    int dev_num;
    CUDA_CHECK(cudaGetDeviceCount(&dev_num));

    if(0 == dev_num) {
        EXIT_WITH_MSG("[!] no device support CUDA, aborting...\n");
    }
    if(verbose) PRINTF_STAMP("\t\tdetected %d GPU devices\n", dev_num);

    const int dev_id = mqs->opts->crossbred_dev_id; // specify device
    if(dev_id >= dev_num) {
        EXIT_WITH_MSG("[!] no CUDA device with such id: %d\n", dev_id);
    }

    cudaDeviceProp dev_prop;
    CUDA_CHECK(cudaGetDeviceProperties(&dev_prop, dev_id));
    CUDA_CHECK(cudaSetDevice(dev_id));
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaDeviceReset());
    if(verbose) {
        PRINTF_STAMP("\t\tusing GPU device %d: %s\n", dev_id, dev_prop.name);
        PRINTF_STAMP("\t\treset device...\n");
        check_gpu(&dev_prop);
    }

    const uint32_t tfnum = svnum - kfnum - kvar_num;
    if(tfnum >= 32) {
        EXIT_WITH_MSG("[!] invalid combination of parameters: \n"
                     "\t\t\tnumber of variables after preprocessing: %" PRIu64 "\n"
                     "\t\t\tnumber of variables to fix in sub-system: %" PRIu64 "\n"
                     "\t\t\tnumber of remaining variables in the sub-system: %" PRIu32 "\n"
                     "\t\t\tkeep variables: %" PRIu64 "\n"
                     "\t\t\tkernel fix variables: %" PRIu32 "\n"
                     "\t\t\tGPU thread search space: 2^%" PRIu32 " (maximal 2^31)\n",
                     var_num, sfnum, svnum, kvar_num, kfnum, tfnum);
    }

    // decide how many thread to launch and the dimension
    const uint32_t thread_num = 0x1U << kfnum;
    const uint32_t warp_num = 1 + (thread_num-1) / WARP_SIZE;
    const uint32_t block_size = WARP_PER_BLOCK * WARP_SIZE;
    const uint32_t block_num = 1 + (thread_num-1) / block_size;
    dim3 gdim(block_num);
    dim3 bdim(block_size);

    const uint32_t cmem_size = gc_cmemsize(tfnum, deg);
    if(cmem_size > MAX_SLOT_NUM * sizeof(uint32_t)) {
        EXIT_WITH_MSG("[!] require too much constant memory: %" PRIu32 " bytes\n",
                      cmem_size);
    }

    const uint32_t max_scnum = mqs->opts->crossbred_mb_size ?
                               mqs->opts->crossbred_mb_size : MAX_SC_NUM;
    const uint32_t mailbox_size = gc_mbsize(max_scnum);

    // compute the number of threads that can be launched
    const uint32_t mempw = gc_wmemsize(tfnum, kvar_num, deg);

    size_t freemem, totalmem;
    CUDA_CHECK(cudaMemGetInfo(&freemem, &totalmem));
    const uint32_t max_warp_num = (freemem - mailbox_size) / mempw;

    if(verbose) {
        PRINTF_STAMP("\t\tavailable amount of memory: %.2f MB\n", freemem / MBFLOAT);
        PRINTF_STAMP("\t\tmailbox size: %" PRIu32 "\n", max_scnum);
        PRINTF_STAMP("\t\tmemory needed for the mailbox: %.2f MB\n",
                    mailbox_size / MBFLOAT);
        PRINTF_STAMP("\t\tmemory needed per warp: %.2f KB\n", mempw / 1024.0f);
        PRINTF_STAMP("\t\tcan launch at most %" PRIu32 " warps (%"
                    PRIu32" threads) on device %d\n",
                    max_warp_num, max_warp_num * WARP_SIZE, dev_id);
        PRINTF_STAMP("\t\tlaunching %" PRIu32 " warps (2^%" PRIu32 " threads)\n",
                    warp_num, kfnum);
        PRINTF_STAMP("\t\tmemory needed for graycode enumeration: %.2f MB\n",
                    (mempw * (double) warp_num) / MBFLOAT);
        PRINTF_STAMP("\t\tdimension of grid: %u x %u x %d\n",
                    gdim.x, gdim.y, gdim.z);
        PRINTF_STAMP("\t\tdimension of block: %u x %u x %d\n",
                    bdim.x, bdim.y, bdim.z);
        PRINTF_STAMP("\t\tsearch space for each thread: 2^%" PRIu32 "\n", tfnum);
    }

    if( (warp_num * (double) mempw) >= freemem - mailbox_size ) {
        EXIT_WITH_MSG("[!] insufficient memory on GPU for graycode enumeration\n");
    }

    if(!mqs->opts->crossbred_rmac_cpu && exp_rmemsize >= (freemem - mailbox_size)) {
        PRINTF_ERR_STAMP("[!] GPU is not expected to have enough memory for "
                         "reduced Macaulay matrix\n");
    }

    /* =========================================================================
     *                              main loop
     * =========================================================================
     */
    threadpool_t* tpool = NULL;
    const uint64_t core_num = get_nprocs();
    const uint64_t thnum = (mqs->opts->thread_num) ? mqs->opts->thread_num : core_num;

    if(verbose) {
        PRINTF_STAMP("[+] initializing threadpool...\n");
        PRINTF_STAMP("\t\tdetected %" PRIu64 " cores, using %" PRIu64 " threads\n",
                     core_num, thnum);
    }

    if( NULL == (tpool = threadpool_create(thnum, thnum, 0)) ) {
        EXIT_WITH_MSG("[!] failed to initialize threadpool\n");
    }

    pthread_mutex_t bf_done_lock = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t bf_done_cond = PTHREAD_COND_INITIALIZER;
    threadpool_t* gpu_master = NULL;
    // at most one GPU master thread and one job
    if( NULL == (gpu_master = threadpool_create(1, 1, 0)) ) {
        EXIT_WITH_MSG("[!] failed to initialize GPU master thread\n");
    }

    FILE* ext_f = NULL;
    if ( mqs->opts->crossbred_mode == MQ_MODE_EXT ) {
        ext_f = fopen(mqs->opts->mq_ext_file, "w");
        if(NULL == ext_f) {
            EXIT_WITH_MSG("[!] cannot create output file: %s\n",
                          mqs->opts->mq_ext_file);
        }
    } else if ( mqs->opts->crossbred_mode == MQ_MODE_BF ) {
        ext_f = fopen(mqs->opts->mq_ext_file, "r");
        if(NULL == ext_f) {
            EXIT_WITH_MSG("[!] cannot read from output file: %s\n",
                          mqs->opts->mq_ext_file);
        }
    }

    bool solvable = false;
    bool bf_done = true;
    const uint64_t loop_num = mqf ? mqf->lnum : 1;
    // for each preprocess configuration
    for(uint64_t idx = 0; idx < loop_num; ++idx) {

        uint32_t* host_32eqs = SMALLOC(uint32_t, mq_tnum);
        uint32_t* subsys = SMALLOC(uint32_t, mac_term_num);
        bool* sys = SMALLOC(bool, mqs->eq_num * mq_tnum);

        if( mqs->opts->crossbred_mode == MQ_MODE_BF ) {
            PRINTF_STAMP("[+] loading sub-systems from file: %s\n",
                        mqs->opts->mq_ext_file);
            parse_subsys(ext_f, subsys, host_32eqs, sys);
        } else {
            if(mqf) {
                PRINTF_STAMP("[+] fixing %" PRIu64 " variables in the original"
                             " MQ system...\n", mqf->fnum);
                fix_mq(sys, reduc_mq, mqf, mqs->var_num, mqs->eq_num,
                       idx, mq_tnum, mqs->xvar_num - mqs->var_num);
            } else {
                memcpy(sys, reduc_mq, sizeof(bool) * mq_tnum * mqs->eq_num);
            }
        
            if(mqs->opts->crossbred_init_gj) {
                PRINTF_STAMP("[+] performing Gauss-Jordan elimination on initial "
                             "MQ system\n");
                gauss_jordan_elim(sys, mqs->eq_num, mq_tnum);
            }

            // pass the first 32 equations into GPU as a filter
            // and keep the rest on the host for verification
            split_sys(host_32eqs, sys, mqs->eq_num, mq_tnum, 32);

            PRINTF_STAMP("[+] computing Macaulay matrix...\n");
            Macaulay* mac = mac_create(mqs->eq_num, var_num, deg);
            if(verbose) {
                PRINTF_STAMP("\t\tdimension of Macaulay matrix: %" PRIu64 " x %" PRIu64 "\n",
                             mac->eq_num, mac->term_num);
                PRINTF_STAMP("\t\tsize of sparse Macaulay matrix: %.2f MB\n",
                             mac_memsize(mqs->eq_num, var_num, deg) / MBFLOAT);
            }

            if(NULL == mac) {
                EXIT_WITH_MSG("[!] insufficient memory\n");
            }
            // compute Macaulay matrix whose columns are in grlex order
            mac_init_perm(mac, sys, mqs->eq_num, var_num, mq_tnum, idx_map, clz);

            if(mqs->opts->crossbred_mac_stats) {
                double* tnums = SMALLOC(double, mac->eq_num);
                double average = 0.0;
                for(uint64_t i = 0; i < mac->eq_num; ++i) {
                    tnums[i] = mac->eqs[i].term_num;
                    average += tnums[i];
                }
                average /= mac->eq_num;
                PRINTF_STAMP("\t\texpected number of terms in a row: %.2f\n",
                             (double)(mq_tnum) * 0.5 );
                PRINTF_STAMP("\t\tmaximal number of terms in a row: %.2f\n",
                             max_f64(tnums, mac->eq_num));
                PRINTF_STAMP("\t\taverage number of terms in a row: %.2f (%.2f)\n",
                            average, std_dev_f64(tnums, mac->eq_num, average));
                free(tnums);
            }

            PRINTF_STAMP("[+] pivoting Macaulay matrix...\n");
            if(verbose) {
                PRINTF_STAMP("\t\tconverting non-pivot rows into dense format...\n");
                PRINTF_STAMP("\t\tconverting pivot rows into compact format...\n");
            }
            const uint64_t nfr_num = mac_pvt(mac, nlterm_num, mac_keq_num, clz, drow,
                                             mq_term_ratio);
            const double dmemsize = mac_drow_memsize(mac->term_num,
                                                     nfr_num + mac_keq_num);
            const double cmemsize = mac_crow_memsize(sum_binom(var_num, 2),
                                                     nlterm_num - nfr_num, mq_term_ratio);
            if(verbose) {
                double mpr_ratio = (double) nfr_num / nlterm_num * 100.0;
                PRINTF_STAMP("\t\tnumber of pivots not found: %" PRIu64 " (%.2f%%)\n",
                             nfr_num, mpr_ratio);
                PRINTF_STAMP("\t\textra memory used for dense rows: %.2f MB\n",
                             dmemsize / MBFLOAT );
                PRINTF_STAMP("\t\textra memory used for compact rows: %.2f MB\n",
                             cmemsize / MBFLOAT );

                if(mpr_ratio <= RMAC_MPR_THRESHOLD) {
                    PRINTF_ERR_STAMP("[!] missing pivot ratio is lower than "
                                     " the threshold\n");
                    PRINTF_ERR_STAMP("\t\tthe extracted sub-system might not "
                                     "behave like a random system\n");
                }
            }

            if(mqs->opts->crossbred_mac_stats) {
                double* snums = SMALLOC(double, nlterm_num - nfr_num);
                uint64_t cursor = 0;
                double average = 0.0;
                for(uint64_t i = 0; i < mac->eq_num; ++i) {
                    if(MAC_CROW == mac->eqs[i].term_num) {
                        snums[cursor] = (double) mac->eqs[i].row.c[0];
                        average += snums[cursor++];
                    }
                }
                average /= cursor;
                PRINTF_STAMP("\t\taverage number of slots for compact row: %.2f (%.2f)\n",
                             average, std_dev_f64(snums, cursor, average));
                free(snums);
            }

            PRINTF_STAMP("[+] computing reduced Macaulay matrix...\n");
            RMac* reduc_mac = rmac_create(nfr_num + mac_keq_num, nfr_num + lterm_num);
            mac_calc_rmac(reduc_mac, mac, lterm_num, drow, nfr_num, mac_keq_num, tpool);
            mac_free(mac);
            mac = NULL;

            if(verbose) {
                PRINTF_STAMP("\t\tdimension of reduced Macaulay matrix: %" PRIu64 " x %"
                             PRIu64 "\n", nfr_num + mac_keq_num, nfr_num + lterm_num);
                double rmemsize = rmac_memsize(nfr_num + mac_keq_num,
                                               nfr_num + lterm_num);
                PRINTF_STAMP("\t\textra memory used for reduced Macaulay "
                             "matrix: %.2f MB\n", rmemsize / MBFLOAT);
            }

            PRINTF_STAMP("[+] performing Gaussian elimination on reduced Macaulay"
                         " matrix...\n");
            if(mqs->opts->crossbred_rmac_cpu) {
                if(verbose) PRINTF_STAMP("\t\tusing CPU...\n");

                rmac_elim_cpu(reduc_mac, nfr_num, tpool);

            } else {
                if(verbose) PRINTF_STAMP("\t\tusing GPU...\n");

                // allocate memory for Gaussian elimination on GPU
                uint64_t* dev_rmac = NULL;
                CUDA_CHECK(cudaMalloc(&dev_rmac,
                           sizeof(uint64_t) * reduc_mac->eq_num * reduc_mac->slot_num));
                
                if(verbose) PRINTF_STAMP("\t\tcopying reduced Macaulay matrix "
                                         "to GPU...\n");
                CUDA_CHECK(cudaMemcpy(dev_rmac, reduc_mac->mem,
                           sizeof(uint64_t) * reduc_mac->eq_num * reduc_mac->slot_num,
                           cudaMemcpyHostToDevice));

                uint32_t* dev_row_indices = NULL;
                CUDA_CHECK(cudaMalloc(&dev_row_indices,
                                      sizeof(uint32_t) * reduc_mac->eq_num));

                uint64_t** dev_rmac_rows = NULL;
                CUDA_CHECK(cudaMalloc(&dev_rmac_rows,
                                      sizeof(uint64_t*) * reduc_mac->eq_num));

                uint32_t* pvt_indices = NULL;
                CUDA_CHECK(cudaMalloc(&pvt_indices,
                                      sizeof(uint32_t) * reduc_mac->eq_num));

                uint32_t* dev_rcount = NULL;
                CUDA_CHECK(cudaMalloc(&dev_rcount, sizeof(uint32_t)));

                uint32_t* host_rcount = NULL;
                CUDA_CHECK(cudaMallocHost(&host_rcount, sizeof(uint32_t)));

                uint32_t* host_row_indices = NULL;
                CUDA_CHECK(cudaMallocHost(&host_row_indices,
                                          sizeof(uint32_t) * reduc_mac->eq_num));

                if(verbose) PRINTF_STAMP("\t\treducing with GPU...\n");
                rmac_elim_gpu(dev_rmac, dev_rmac_rows, dev_row_indices, reduc_mac->eq_num,
                              reduc_mac->slot_num, nfr_num, pvt_indices, dev_rcount,
                              host_rcount);

                if(verbose) {
                    PRINTF_STAMP("\t\tcopying reduced Macaulay matrix back "
                                 "to the host...\n");
                }
                CUDA_CHECK(cudaMemcpy(reduc_mac->mem, dev_rmac,
                           sizeof(uint64_t) * reduc_mac->eq_num * reduc_mac->slot_num,
                           cudaMemcpyDeviceToHost));

                CUDA_CHECK(cudaMemcpy(host_row_indices, dev_row_indices,
                       sizeof(uint32_t) * reduc_mac->eq_num, cudaMemcpyDeviceToHost));

                if(verbose) PRINTF_STAMP("\t\trestoring row order on the host...\n");
                rmac_perm_rows(reduc_mac, host_row_indices);

                CUDA_CHECK(cudaFree(dev_rmac));
                CUDA_CHECK(cudaFree(pvt_indices));
                CUDA_CHECK(cudaFree(dev_rcount));
                CUDA_CHECK(cudaFreeHost(host_rcount));
                CUDA_CHECK(cudaFree(dev_rmac_rows));
                CUDA_CHECK(cudaFreeHost(host_row_indices));
                CUDA_CHECK(cudaFree(dev_row_indices));
            }

            PRINTF_STAMP("[+] packing sub-system...\n");
            uint64_t valid_eq_num = verify_subsys(valid_eq_indices, reduc_mac, nfr_num,
                                                  mac_keq_num);

            // choose at most 32 eqs and transpose the sub-system
            const uint64_t sugg_keq_num = MIN(mqs->opts->crossbred_sub_keq_num, 32);
            const uint64_t keq_num =  MIN(sugg_keq_num, valid_eq_num);
            if(verbose) {
                PRINTF_STAMP("\t\tdimension of sub-system: %" PRIu64 " x %" PRIu64 "\n",
                             valid_eq_num, lterm_num);
                PRINTF_STAMP("\t\tinitial sparsity of the sub-system: %.5f\n",
                            sparsity(reduc_mac, mac_keq_num, lterm_num));
                PRINTF_STAMP("\t\tnumber of invalid equations: %" PRIu64 "\n",
                             mac_keq_num - valid_eq_num);
                PRINTF_STAMP("\t\tcan extract %" PRIu64 " equations from"
                             " reduced Macaulay matrix\n", valid_eq_num);
                PRINTF_STAMP("\t\tsuggested to keep %" PRIu64 " equations\n", sugg_keq_num);
                PRINTF_STAMP("\t\tkeeping %" PRIu64 " equations...\n", keq_num);
            }

            if(valid_eq_num < kvar_num) {
                PRINTF_ERR_STAMP("\t\t[!] can only extract %" PRIu64
                                " equations, skipping\n", valid_eq_num);
                rmac_free(reduc_mac);
                free(sys);
                free(host_32eqs);
                free(subsys);
                continue;
            }

            // pack transposed sub-system. We keep at most 64 eqs
            pick_subsys(subsys, reduc_mac, mac_term_num, keq_num, valid_eq_num,
                        valid_eq_indices, l_indices, lterm_num);

            rmac_free(reduc_mac);
            reduc_mac = NULL;
        } // end of extracting a sub-system

        if ( mqs->opts->crossbred_mode == MQ_MODE_EXT ) {
            PRINTF_STAMP("[+] dumping sub-system into %s\n", mqs->opts->mq_ext_file);
            dump_subsys(ext_f, subsys, mac_term_num, mqf, idx, host_32eqs, mq_tnum,
                        sys, mqs->eq_num);
            free(host_32eqs);
            free(subsys);
            free(sys);
            continue;
        }

        if(verbose) {
            PRINTF_STAMP("[+] checking the status of GPU master thread...\n");
        }
        // wait until the previous batch finishes
        pthread_mutex_lock(&bf_done_lock);
        while(!bf_done) {
            pthread_cond_wait(&bf_done_cond, &bf_done_lock);
        }
        bf_done = false;
        pthread_mutex_unlock(&bf_done_lock);

        // check if previous batch yielded a solution
        // access directly because the previous batch is done so no other
        // thread will modify solvable.
        if(solvable) {
            break;
        }

        // parepare the arguments
        bf_subsys_arg* arg = SMALLOC(bf_subsys_arg, 1);
        arg->dev_id = dev_id;
        arg->mac_term_num = mac_term_num;
        arg->mq_tnum = mq_tnum;
        arg->var_num = var_num;
        arg->thread_num = thread_num;
        arg->sfnum = sfnum;
        arg->deg = deg;
        arg->kvar_num = kvar_num;
        arg->idx = idx;
        arg->mempw = mempw;
        arg->warp_num = warp_num;
        arg->mailbox_size = mailbox_size;
        arg->max_scnum = max_scnum;
        arg->cmem_size = cmem_size;
        arg->fsub_tnum = fsub_tnum;
        arg->svnum = svnum;
        arg->tfnum = tfnum;
        arg->kfnum = kfnum;
        arg->gdim = gdim;
        arg->bdim = bdim;
        arg->host_32eqs = host_32eqs;
        arg->subsys = subsys;
        arg->solution = solution;
        arg->sys = sys;
        arg->mqf = mqf;
        arg->mqs = mqs;
        arg->verbose = verbose;
        arg->bf_done_lock= &bf_done_lock;
        arg->bf_done_cond = &bf_done_cond;
        arg->bf_done = &bf_done;
        arg->solvable = &solvable;

        // start bruteforcing
        if(mqs->opts->crossbred_rmac_cpu) {
            if(verbose) {
                PRINTF_STAMP("[+] pipe data to GPU then resume Macaulay matrix "
                             "computation...\n");
            }

            if(threadpool_add(gpu_master, bf_subsys, (void*) arg, 0)) {
                EXIT_WITH_MSG("[!] failed to assign work to the GPU master thread\n");
            }
        }
        else { // sequential
            bf_subsys(arg);
        }
    } // loop over preprocess configs

    // wait for the gpu orchestrater
    if(mqs->opts->crossbred_rmac_cpu && threadpool_join(gpu_master, 0)) {
        EXIT_WITH_MSG("[!] failed to join with the GPU master thread\n");
    }

    if(verbose) PRINTF_STAMP("[+] destroying threadpool...\n");
    if( 0 != threadpool_destroy(tpool, 0) ) {
        PRINTF_ERR_STAMP("[!] error when destroying the threadpool\n");
    }

    if( 0 != threadpool_destroy(gpu_master, 0) ) {
        PRINTF_ERR_STAMP("[!] error when destroying the threadpool\n");
    }

    if (ext_f) {
        fclose(ext_f);
        ext_f = NULL;
    }

    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaDeviceReset());
    mqfix_free(mqf);
    free(l_indices);
    free(nl_indices);
    free(valid_eq_indices);
    free(clz);
    free(drow);
    free(idx_map);
    free(reduc_mq);
    pthread_mutex_destroy(&bf_done_lock);
    pthread_cond_destroy(&bf_done_cond);
    return solvable;
}

/* function: mqs_solve
 * usage: try to solve the MQ system. The algorithm to be used depends on the
 *       options passed via command line. Once a solution is found, the
 *       function prints it to the terminal and returns.
 * arguments: a pointer to structure MQSolver
 * return: 0 is a solution is found, otherwise 1.
 */
int
mqs_solve(MQSolver* mqs) {
    if(NULL == mqs->orig_sys || NULL == mqs->opts) {
        EXIT_WITH_MSG("[!] no linear system given or no options given");
    }

    bool solution[mqs->var_num];
    bool (*selec_algor) (MQSolver*, bool*) = NULL;

    /* based on the option choose algorithm accordingly */
    switch(mqs->opts->algor) {
        case MQ_ALG_FAST_EX:
            selec_algor = fast_exhaustive;
            break;
        case MQ_ALG_CROSSBRED:
            selec_algor = crossbred;
            break;
        case MQ_ALG_INVALID:
        default:
            /* do nothing */
            break;
    }

    if(NULL == selec_algor) {
        EXIT_WITH_MSG("[!] invalid algorithm choice\n");
    }

    int code = 0;
    if(selec_algor(mqs, solution)) {
        PRINTF_STAMP("[+] solution found: \n");

        uint64_t i;
        printf("\t\t\t[%d", solution[0]);
        for(i = 1; i < mqs->var_num; i++) {
            printf(", %d", solution[i]);
        }
        printf("]\n");

    } else {
        PRINTF_STAMP("[-] solution not found\n");
        code = 1;
    }
    
    return code;
}

} /* extern "C" */
