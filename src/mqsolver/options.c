/* options.c: implementation of structure Options
 */

#include "options.h"

/* function: print_usage
 * usage: print help message
 * arguments: prg_name, name of the executable
 * return: void
 */
static void
print_usage(char* prg_name) {
   printf("\n"
   "Usage: %s [OPTIONS]\n"
   "\n"
   "\n"
   "Options:                    \n"
   "                            \n"
   "  --seed=SEED           Use SEED to initialize the random number sequence.\n"
   "                            Default seed is random.\n"
   "                            \n"
   "  --challenge=FILE      Read MQ challenge from FILE.\n"
   "                            \n"
   "  --algorithm=ALGOR     Which algorithm to use. Choices are:\n"
   "                            fast_ex\n"
   "                            crossbred\n"
   "                            \n"
   "  --verbose             Print extra information.\n"
   "                            \n"
   "Options for crossbred algorithm:\n"
   "                            \n"
   "  --macaulay_deg=DEGREE Degree of the Macaulay matrix. Can be 3 or 4.\n"
   "                            \n"
   "  --keep_var=NUM        Number of variables to keep after linearization.\n"
   "                            Must be no larger than the number of variables\n"
   "                            in the sub-system after fixing variables\n"
   "                            \n"
   "  --kf_var=NUM          Log2 of the number of threads to launch, e.g.\n"
   "                            for 1024 threads, pass 10. This also decides the\n"
   "                            number of variables to be fixed by the kernel.\n"
   "                            \n"
   "  --sub_fvr=NUM         Number of variables to fix after a sub-system is obtained.\n"
   "                            Default is 0.\n"
   "                            \n"
   "  --mac_stats           Collect stats about Macaulay matrix.\n"
   "                            \n"
   "  --mq_ratio=RATIO      The maximal density of each row of the initial MQ system.\n"
   "                            Used to adjust the memory requirement for computing\n"
   "                            reduced Macaulay matrix. Default is 0.5.\n"
   "                            \n"
   "  --dev_id=NUM          Specify the id of the CUDA device to use. Default is 0.\n"
   "                            \n"
   "  --thread_num=NUM      Specify the number of CPU threads to use.\n"
   "                            Default is the number of detected hyper-threads.\n"
   "                            \n"
   "  --mq_fix=FILE             Preprocess the MQ system by fixing variables as\n"
   "                            specified by FILE. The first line of FILE contain\n"
   "                            the number of variables to fix and the number of\n"
   "                            configurations, e.g. 4 16. Each line of the rest of\n"
   "                            the file contains a configuration, which is a\n"
   "                            tuple of variable separated by space. The value for the\n"
   "                            last variable should be the last entry in the tuple,\n"
   "                            e.g. 1 1 0 fixes x_n = 0, x_{n-1} = 1, x_{n-2} = 1.\n"
   "                            \n"
   "                            Alternatively, if the file has one single line\n"
   "                            containing 3 intergers, e.g. 3 4 1,\n"
   "                            which specifies the number of variables to fix, the\n"
   "                            number of preprocess configurations, and the starting\n"
   "                            preproccess configuration.\n"
   "                            \n"
   "  --mac_keq=NUM         Number of equations to keep as sub-system row candidates.\n"
   "                            Default is 64.\n"
   "                            \n"
   "  --sub_keq=NUM         Number of equations to keep in the sub-system. Default is"
   "                            32 and maximal is 32.\n"
   "                            \n"
   "  --rmac_cpu            Use CPU instead of GPU to perform Gaussian elimination\n"
   "                            on reduced Macaulay matrix. This may be necessary\n"
   "                            when the size of reduced Macaulay matrix does not\n"
   "                            fix into GPU memory.\n"
   "                            \n"
   "  --mailbox_size=NUM    Number of solution candidates to keep in the mailbox\n"
   "                            Default is 2^15 and maximal is most 2^32-1.\n"
   "                            \n"
   "  --init_gj=true|false  Perform Gauss-Jordan elimination on the initial MQ system\n"
   "                            Default to true.\n"
   "                            \n"
   "  --mode=MODE           Mode of MQsolver. Can be:\n"
   "                            default: simply solve the MQ system\n"
   "                            extsub: extract sub-systems and output them into a file\n"
   "                            bruteforce: read extracted sub-systems from a file and\n"
   "                            bruteforce with them\n"
   "                            \n"
   "  --mq_ext_file=FILE    Output file for extsub mode.\n"
   "                            \n"
   "Examples:                   \n"
   "                            \n"
   "  %s --algorithm=crossbred --challenge=toy_example.txt --macaulay_deg=3 --keep_var=5\n"
   "    --kf_var=5              \n"
   "                            \n"
   "  %s --algorithm=crossbred --challenge=toy_example.txt --macaulay_deg=4\n"
   "    --keep_var=18 --mac_keq=32 --mailbox_size=1000\n"
   "                            \n"
   "\n", prg_name, prg_name, prg_name);
}

/* function: options_init
 * usage: init the Options structure
 * arguments: an pointer to structure Options
 * return: void
 */
void
options_init(Options* opts) {
    opts->cha_file = NULL;
    opts->mq_ext_file = NULL;

    opts->algor = MQ_ALG_INVALID;
    opts->seed = time(NULL);

    opts->crossbred_kvar_num = 0;
    opts->crossbred_mac_degree = 0;
    opts->crossbred_mac_stats = false;
    opts->crossbred_mac_keq_num = 64;
    opts->crossbred_mb_size = 0;
    opts->crossbred_sub_fnum = 0;
    opts->crossbred_dev_id = 0;
    opts->crossbred_mq_ratio = 0.5;
    opts->crossbred_sub_keq_num = 32;
    opts->crossbred_mode = MQ_MODE_DEFAULT;

    opts->mq_fix_file = NULL;
    opts->verbose = false;
    opts->crossbred_rmac_cpu = false;
    opts->crossbred_init_gj = true;
    opts->thread_num = 0;
}

/* function: options_free
 * usage: free the Options structure
 * arguments an pointer to structure Options
 * return: void
 */
void
options_free(Options* opts) {
    SFREE(opts->cha_file);
    SFREE(opts->mq_fix_file);
    SFREE(opts->mq_ext_file);
}

/* function: options_parse
 * usage: parse argv and store the options into structure Options
 *      Make sure options_init has been called.
 * arguments:
 *      1) opts: an pointer to the Options structure
 *      2) argc
 *      3) argv
 * return: void
 */

/* for getopt_long */
#define OP_SEED                 1
#define OP_CHALLENGE            2
#define OP_ALGORITHM            3
#define OP_CROSSBRED_MAC_DEG    4
#define OP_CROSSBRED_KEEP_VAR   5
#define OP_CROSSBRED_KF_VAR     6
#define OP_CROSSBRED_MAC_STATS  7
#define OP_CROSSBRED_MQ_FIX     8
#define OP_CROSSBRED_MAC_KEQ    9
#define OP_VERBOSE              10
#define OP_CROSSBRED_MB_SZ      11
#define OP_CROSSBRED_SUB_FVAR   12
#define OP_CROSSBRED_DEV_ID     13
#define OP_CROSSBRED_RMAC_CPU   14
#define OP_CROSSBRED_MQ_RATIO   15
#define OP_THREAD_NUM           16
#define OP_CROSSBRED_SUB_KEQ    17
#define OP_CROSSBRED_INIT_GJ    18
#define OP_CROSSBRED_MODE       19
#define OP_MQ_EXT_FILE          20

static char* const bool_true = "true";
static char* const bool_false = "false";

static char* const mode_default = "default";
static char* const mode_extract = "extsub";
static char* const mode_bruteforce = "bruteforce";

static struct option mq_long_opts[] = {
    { "challenge", 1, 0, OP_CHALLENGE },
    { "algorithm", 1, 0, OP_ALGORITHM },
    { "seed", 1, 0, OP_SEED },
    { "thread_num", 1, 0, OP_THREAD_NUM },
    { "verbose", 0, 0, OP_VERBOSE },

    /* options for crossbred */
    { "macaulay_deg", 1, 0, OP_CROSSBRED_MAC_DEG },
    { "keep_var", 1, 0, OP_CROSSBRED_KEEP_VAR },
    { "kf_var", 1, 0, OP_CROSSBRED_KF_VAR },
    { "mq_fix", 1, 0, OP_CROSSBRED_MQ_FIX },
    { "mac_keq", 1, 0, OP_CROSSBRED_MAC_KEQ },
    { "mailbox_size", 1, 0, OP_CROSSBRED_MB_SZ },
    { "sub_fvar", 1, 0, OP_CROSSBRED_SUB_FVAR },
    { "dev_id", 1, 0, OP_CROSSBRED_DEV_ID },
    { "rmac_cpu", 0, 0, OP_CROSSBRED_RMAC_CPU },
    { "mac_stats", 0, 0, OP_CROSSBRED_MAC_STATS },
    { "mq_ratio", 1, 0, OP_CROSSBRED_MQ_RATIO },
    { "sub_keq", 1, 0, OP_CROSSBRED_SUB_KEQ },
    { "init_gj", 1, 0, OP_CROSSBRED_INIT_GJ },
    { "mode", 1, 0, OP_CROSSBRED_MODE },
    { "mq_ext_file", 1, 0, OP_MQ_EXT_FILE },
    
    { "help", 0, 0, 'h' },
    { 0, 0, 0, 0 }
};

/* copy parsed option into a char* pointer, take at most 1024 chars */
static inline void
copy_opt(char** str, char* optarg) {
    if(NULL == ((*str) = strndup(optarg, 1024))) {
        EXIT_WITH_MSG("[!] invalid input file");
    }
}

void
options_parse(Options* opts, int argc, char** argv) {
    int c, opt_idx;
    while(-1 != (c = getopt_long(argc, argv, "h", mq_long_opts, &opt_idx))) {
        switch(c) {
            case 0:
                /* If opts option set a flag, don't do anything */
                if(mq_long_opts[opt_idx].flag == 0) {
                    PRINTF_STAMP("option %s: %s\n", mq_long_opts[opt_idx].name,
                                 optarg ? optarg : "null");
                }
                break;

            case 'h':
                print_usage(argv[0]);
                SFREE(opts);
                exit(0);
                break;

            case OP_VERBOSE:
                PRINTF_STAMP("option: verbose\n");
                opts->verbose = true;
                break;

            case OP_SEED:
                opts->seed = strtol(optarg, NULL, 0);
                PRINTF_STAMP("option seed: %u\n", opts->seed);
                break;

            case OP_CHALLENGE:
                copy_opt(&opts->cha_file, optarg);
                PRINTF_STAMP("option challenge: %s\n", opts->cha_file);
                break;

            case OP_ALGORITHM:
                opts->algor = algor2code(optarg);
                PRINTF_STAMP("option algorithm: %s\n", optarg);
                if(MQ_ALG_INVALID == opts->algor) {
                    EXIT_WITH_MSG("[!] invalid algorithm\n");
                }
                break;

            case OP_CROSSBRED_KEEP_VAR:
                opts->crossbred_kvar_num = strtoul(optarg, NULL, 0);
                PRINTF_STAMP("option keep_var: %"PRIu64"\n",
                        opts->crossbred_kvar_num);
                break;

            case OP_CROSSBRED_MAC_DEG:
                opts->crossbred_mac_degree = strtoul(optarg, NULL, 0);
                PRINTF_STAMP("option macaulay_deg: %"PRIu64"\n",
                        opts->crossbred_mac_degree);
                break;

            case OP_CROSSBRED_KF_VAR:
                opts->crossbred_kf_num = strtoul(optarg, NULL, 0);
                PRINTF_STAMP("option kf_var: %"PRIu64"\n",
                        opts->crossbred_kf_num);
                break;

            case OP_CROSSBRED_MAC_STATS:
                opts->crossbred_mac_stats = true;
                PRINTF_STAMP("option mac_stats: on\n");
                break;

            case OP_CROSSBRED_MQ_FIX:
                copy_opt(&opts->mq_fix_file, optarg);
                PRINTF_STAMP("option mq_fix: %s\n", opts->mq_fix_file);
                break;

            case OP_CROSSBRED_MAC_KEQ:
                opts->crossbred_mac_keq_num = strtoul(optarg, NULL, 0);
                PRINTF_STAMP("option mac_keq: %"PRIu64"\n",
                             opts->crossbred_mac_keq_num);
                break;

            case OP_CROSSBRED_MB_SZ:
                opts->crossbred_mb_size = strtoul(optarg, NULL, 0);
                PRINTF_STAMP("option mailbox_size: %"PRIu64"\n",
                             opts->crossbred_mb_size);
                break;

            case OP_CROSSBRED_SUB_FVAR:
                opts->crossbred_sub_fnum = strtoul(optarg, NULL, 0);
                PRINTF_STAMP("option sub_fvar: %"PRIu64"\n",
                             opts->crossbred_sub_fnum);
                break;

            case OP_CROSSBRED_DEV_ID:
                opts->crossbred_dev_id = atoi(optarg);
                PRINTF_STAMP("option dev_id: %d\n", opts->crossbred_dev_id);
                break;

            case OP_CROSSBRED_RMAC_CPU:
                opts->crossbred_rmac_cpu = true;
                PRINTF_STAMP("option rmac_cpu: true\n");
                break;

            case OP_CROSSBRED_MQ_RATIO:
                opts->crossbred_mq_ratio = strtod(optarg, NULL);
                PRINTF_STAMP("option mq_ratio: %f\n", opts->crossbred_mq_ratio);
                break;

            case OP_THREAD_NUM:
                opts->thread_num = strtoul(optarg, NULL, 0);
                PRINTF_STAMP("option thread_num: %"PRIu64"\n", opts->thread_num);
                break;

            case OP_CROSSBRED_SUB_KEQ:
                opts->crossbred_sub_keq_num = strtoul(optarg, NULL, 0);
                PRINTF_STAMP("option sub_keq: %"PRIu64"\n",
                             opts->crossbred_sub_keq_num);
                break;

            case OP_CROSSBRED_INIT_GJ:
                if (!strncmp(optarg, bool_true, strlen(bool_true))) {
                    opts->crossbred_init_gj = true;
                } else if (!strncmp(optarg, bool_false, strlen(bool_false))) {
                    opts->crossbred_init_gj = false;
                } else {
                    EXIT_WITH_MSG("[!] invalid option: %s\n", optarg);
                }

                PRINTF_STAMP("option init_gj: %s\n", optarg);
                break;

            case OP_CROSSBRED_MODE:
                if ( !strncmp(optarg, mode_extract, strlen(mode_extract)) ) {
                    opts->crossbred_mode = MQ_MODE_EXT;
                } else if ( !strncmp(optarg, mode_bruteforce, strlen(mode_bruteforce)) ) {
                    opts->crossbred_mode = MQ_MODE_BF;
                } else if ( !strncmp(optarg, mode_default, strlen(mode_default))) {
                    opts->crossbred_mode = MQ_MODE_DEFAULT;
                } else {
                    EXIT_WITH_MSG("[!] invalid option: %s\n", optarg);
                }

                PRINTF_STAMP("option mode: %s\n", optarg);
                break;

            case OP_MQ_EXT_FILE:
                copy_opt(&opts->mq_ext_file, optarg);
                PRINTF_STAMP("option mq_ext_file: %s\n", opts->mq_ext_file);
                break;

            case '?':
                /* getopt_long already printed an error message */
                break;

            default:
                EXIT_WITH_MSG("[!] unknown error\n");
        }
    }
}
