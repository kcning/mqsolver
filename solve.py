#!/usr/bin/env python3
# usage: compile and launch MQsolver according to the parameters.
# Run -h for details

import argparse
import sys
import subprocess
import os
import time
import itertools
import shutil

start_stamp = time.time()

def print_stamp(*arg, file=sys.stdout):
    stamp = time.time() - start_stamp
    print("%12.2f - " % stamp, end='', file=file)
    print(*arg, file=file)

def str2bool(b):
    if b.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif b.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


prog_name = 'MQsolver'
parser = argparse.ArgumentParser(prog=prog_name)
parser.add_argument("-s", "--seed", help="use the seed to initialize the "
                    "random number sequence. Default seed is random.", type=int)

parser.add_argument("-a", "--algorithm", help="which algorithm to use. "
                    "Can be fast_ex or crossbred. "
                    "Default is %(default)s.", type=str, default='crossbred')

parser.add_argument("MQfile", help="input MQ challenge file.",
                    type=argparse.FileType('r'))

parser.add_argument("-v", "--verbose", help="print extra information.",
                    action='store_true')

parser.add_argument("-o", "--output-file", help="dump the output, both stdout "
                    "and stderr to the file.", type=argparse.FileType('w'))

parser.add_argument('--verbose-compile', help="show compile messages",
                    action='store_true')

exe_cmd = parser.add_mutually_exclusive_group()

exe_cmd.add_argument("-C", "--compile",
                     help="compile but do not invoke %(prog)s",
                     action='store_true')

exe_cmd.add_argument("-X", "--execute", help="do no compile and use existing "
                     "%(prog)s executable from previous compilation. It is the "
                     "user's responsibility to ensure the number of variables "
                     "to keep remains the same.", action='store_true')

parser.add_argument("-H", "--host-file", help="launch %(prog)s in cluster "
                    "mode. This machine acts as the master node while the "
                    "machines in the host file server as slave nodes. The "
                    "master node copies the MQ challenge file and distributes "
                    "the workload evenly to the slave nodes. Once a slave node "
                    "finds a solution. The master node kills the rest of "
                    "the slave nodes. "
                    "The host file consists of lines, where each line is in "
                    "the following format: "
                    "<user name>@<ip address | domain name>:<port>",
                    type=argparse.FileType('r'))

cluster_mode_opts = parser.add_argument_group('cluster mode',
                                              'options for cluster mode')

cluster_mode_opts.add_argument("--cluster-stdout", help="the directory where "
                               "the stdout output from the slave nodes will be "
                               "stored. Default is %(default)s",
                               default='mq_stdout', type=str)

cluster_mode_opts.add_argument("--cluster-stderr", help="the directory where "
                               "the stderr output from the slave nodes will be "
                               "stored. Default is %(default)s",
                               default='mq_stderr', type=str)

cluster_mode_opts.add_argument("-K", "--key-file",
                               help="which ssh key file to use for "
                               "the cluster mode.", type=argparse.FileType('r'))

cluster_mode_opts.add_argument("-W", "--work-dir",
                               help="the directory on a slave node where "
                               "the %(prog)s executable locates.", type=str)

cluster_mode_opts.add_argument("-N", "--gpu-num", help="number of GPUs on a "
                               "slave node in the cluster. Default is %(default)s.",
                               type=int, default=1)

cluster_mode_opts.add_argument("-P", "--poll-period", help="period of polling "
                               "slave nodes for result. Default is %(default)s "
                               "seconds.", type=int, default=1800)

cluster_mode_opts.add_argument("--email-notice", help="send a notification to "
                               "the email address when a solution is found or "
                               "when the cluster is done. This relies on a "
                               "configured mutt email client.", type=str)

cluster_mode_opts.add_argument("--ignore-fail", help="do not abort when a "
                               "slave node in the cluster failed.",
                               action='store_true')

parser.add_argument("--no-clean", help="do not delete temporary files "
                    "and executable after %(prog)s finishes",
                    action='store_true')

crossbred_opts = parser.add_argument_group('crossbred',
                                           'options for crossbred algorithm')
crossbred_opts.add_argument('-d', '--degree', help="degree of the Macaulay "
                            "matrix. Can be 3 or 4. (Default is %(default)s).",
                            default=3, type=int)

crossbred_opts.add_argument('-k', '--keep-var-num', help="number of variables "
                            "to keep after linearization. Should be no larger "
                            "than the number of variables in the sub-system "
                            "after fixing variables.", type=int)

crossbred_opts.add_argument('-t', '--thread-fix-num', help="number of "
                            "variables to fix by a GPU kernel. The number of "
                            "GPU threads to launch is therefore 2^t. For "
                            "example, passing 10 will launch 2^10 GPU threads.",
                            type=int)

crossbred_opts.add_argument('-T', '--cpu-thread-num', help="number of CPU "
                            "threads to use. Default is the number of detected "
                            "logical cores (hyper-threads).", type=int)

crossbred_opts.add_argument('-f', '--subsys-fix-num', help="number of "
                            "variables to fix in a subsystem before launching "
                            "GPU kernel. Default is %(default)s.", type=int,
                            default=0)

mq_process = crossbred_opts.add_mutually_exclusive_group()

mq_process.add_argument('-q', '--mq-fix-num', help="number of variables to fix "
                        "in the MQ system before computing Macaulay matrix. "
                        "Default is %(default)s. This option is mandatory in "
                        "cluster mode.", type=int, default=0)

mq_process.add_argument('-Q', '--mq-fix-file', help="fix variables in the MQ "
                        "system before computing Macaulay matrix according to "
                        "the configuration file, where the first line contains "
                        "two integer, e.g. 3 2, which specifies the number of "
                        "variables to fix and the number of proprocess "
                        "configurations. Each of the remaining line is a tuple "
                        "of 1 and 0 that specifies one preprocess "
                        "configuration, e.g. "
                        "1 1 0 fixes x_n = 0, x_{n-1} = 1, x_{n-2} = 1. "
                        "Alternatively, the file can contain one single line "
                        "consisting of 3 intergers, e.g. 3 4 1, "
                        "which specifies the number of variables to fix, the "
                        "number of preprocess configurations, and the starting "
                        "preproccess configuration. "
                        "This option is ignored in cluster mode.",
                        type=str)

crossbred_opts.add_argument('-g', '--gpu-dev-id', help="id of the GPU device "
                            "to use. Default is %(default)s. This option is "
                            "ignored in cluster mode.", default=0, type=int)

crossbred_opts.add_argument('-M', '--macaulay-stats', help="Collect statistics "
                            "about Macaulay matrix.", action='store_true')

crossbred_opts.add_argument('-R', '--mq-max-ratio', help="maximal density of "
                            "each row in the matrix representing the MQ system "
                            "after preprocessing. Default is %(default)s.",
                            type=float, default=0.5)

crossbred_opts.add_argument('-e', '--subsys-candidate-keep-num', help="number "
                            "of equations to keep as sub-system candidates. "
                            "Default is %(default)s.", default=64, type=int)

crossbred_opts.add_argument('-E', '--subsys-eq-keep-num', help="number of "
                            "equations to keep in the sub-system. Default is "
                            "%(default)s and maximum is 32.", default=32, type=int)

crossbred_opts.add_argument('-c', '--cpu-reduction', help="use CPU "
                            "to perform Gaussian elimination on reduced "
                            "Macaulay matrix. This might be necessary when the "
                            "size of the reduced Macaulay matrix does not fix "
                            "into GPU memory.", action='store_true')

crossbred_opts.add_argument('-m', '--mailbox-num', help="number of solution "
                            "candidates to keep in the mailbox. Default is "
                            "%(default)s and maximum is 2^32-1.",
                            default=(2 ** 15), type=int)

crossbred_opts.add_argument('-L', '--count-linear-dep', help="count the number "
                            "of underdetermined linear systems obtained by "
                            "linearize the sub-system. This increases runtime "
                            "considerably.", action='store_true')

crossbred_opts.add_argument('-J', '--init-gj', help="specify if Gauss-Jordan "
                            "elimination should be performed on the initial MQ "
                            "system. Default to %(default)s.", default=True,
                            type=str2bool)

crossbred_opts.add_argument('--mq-ext-file', help="Output file for extsub and "
                            " bruteforcce mode. "
                            "Ignored if in default mode.", type=str)

crossbred_opts.add_argument('--mode', help="Mode of MQsolver. Can be "
                            "default: simply solve the MQ system, "
                            "extsub: (experimental) extract sub-systems and output "
                            "them into a file, "
                            "bruteforce: (experimental) read extracted sub-systems "
                            "from a file and "
                            "bruteforce with them.", type=str, default='default')

args = parser.parse_args()

build_dir = 'build'

if args.algorithm == 'crossbred':
    if not args.keep_var_num:
        print_stamp("[!] missing option: -k", file=sys.stderr)
        sys.exit(1)

    if not args.thread_fix_num:
        print_stamp("[!] missing option: -t", file=sys.stderr)
        sys.exit(1)


def check_subproc(proc, verbose=False, exit=True):
    if proc.returncode != 0:
        print_stamp("[!] subprocess failed")
        result = proc.stderr.readlines()
        for line in result:
            print_stamp('\t' + line.decode('utf-8'), file=sys.stderr)

        if exit:
            sys.exit(1)
    else:
        if verbose:
            result = proc.stdout.readlines()
            for line in result:
                print_stamp('\t' + line.decode('utf-8'), file=sys.stderr)


def format_key_opt(args):
    return "-i %s" % args.key_file.name if args.key_file else ''


def format_slave_mqcha(args):
    return "%s/%s" % (args.work_dir, os.path.basename(args.MQfile.name))


def format_scp_cmd(host, port, args):
    key_opt = format_key_opt(args)
    return "scp %s -P %s %s %s:%s" % (key_opt, port, args.MQfile.name,
                                         host, format_slave_mqcha(args))


def format_compile_cmd(host, port, args):
    key_opt = format_key_opt(args)
    dep_opt = "-DGC_DEPC_LSYS" if args.count_linear_dep else ''
    cmd = "ssh %s -p %s %s 'cd %s && mkdir -p %d && cd %d" % (key_opt, port, host,
            args.work_dir, build_dir, build_dir)
    cmd += "&& cmake -DKEEP_VAR_NUM=%d -DCMAKE_C_FLAGS=\"%s\" .. && make" % (
            args.keep_var_num, dep_opt)
    return cmd


def prep_fname(prep_fnum, id):
    return "prep%d-%d.txt" % (prep_fnum, id)


def to_binary(num, digit_num):
    return ''.join(list("{0:b}".format(num).zfill(digit_num)))


def format_gen_prep_cmd(host, port, args, id, start, chunk_size):
    key_opt = format_key_opt(args)
    prep_file = prep_fname(args.mq_fix_num, id)
    cmd = "ssh %s -p %s %s 'cd %s && ./bin/gen_efix.py %s %d %d %d'" % (key_opt,
           port, host, args.work_dir, prep_file, args.mq_fix_num,
           start, chunk_size)
    return cmd

def mq_log_fname(host, port, id):
    return "%s-%s-%d.txt" % (host, port, id)


def format_mq_cmd(host, port, args, id):
    key_opt = format_key_opt(args)
    verbose_opt = '--verbose' if args.verbose else ''
    seed_opt = '--seed=%d' if args.seed else ''
    mac_stats_opt = '--mac_stats' if args.macaulay_stats else ''
    prep_file = prep_fname(args.mq_fix_num, id)
    cpu_reduc_opt = '--rmac_cpu' if args.cpu_reduction else ''
    cpu_thread_opt = '--thread_num=%d' % args.cpu_thread_num if \
            args.cpu_thread_num else ''
    ofile = mq_log_fname(host, port, id)

    cmd = "ssh %s -p %s %s 'screen -dmS mqsolver_session%d " % (key_opt, port, host, id)
    cmd += "bash -c \"cd %s && ./mqsolver --challenge=%s " % (args.work_dir,
            format_slave_mqcha(args))
    cmd += "--algorithm=%s --macaulay_deg=%d %s --keep_var=%d --kf_var=%d " % (
            args.algorithm, args.degree, verbose_opt, args.keep_var_num,
            args.thread_fix_num)
    cmd += "--mac_keq=%d --mailbox_size=%d %s %s --sub_fvar=%d " % (
            args.subsys_candidate_keep_num, args.mailbox_num, seed_opt,
            mac_stats_opt, args.subsys_fix_num)
    cmd += "--dev_id=%d --mq_fix=%s %s --mq_ratio=%f " % (id, prep_file,
            cpu_reduc_opt, args.mq_max_ratio)
    cmd += "%s --sub_keq=%d &> %s\"'" % (cpu_thread_opt, args.subsys_eq_keep_num,
            ofile)
    return cmd


def format_scp_stats_cmd(host, port, args, id):
    key_opt = format_key_opt(args)
    ofile = mq_log_fname(host, port, id)
    return "scp %s -P %s %s:%s/%s %s/%s" % (key_opt, port, host, args.work_dir,
            ofile, args.cluster_stdout, ofile)


def format_check_cmd(host, port, args):
    key_opt = format_key_opt(args)
    return "ssh %s -p %s %s 'pgrep mqsolver'" % (key_opt, port, host)


def format_kill_cmd(host, port, args, id):
    key_opt = format_key_opt(args)
    return "ssh %s -p %s %s 'screen -X -S mqsolver_session%d quit'" % (key_opt,
            port, host, id)


def format_clean_cmd(host, port, args):
    key_opt = format_key_opt(args)
    cmd = "ssh %s -p %s %s 'cd %s && rm -rf %s && " % (key_opt,
            port, host, args.work_dir, build_dir)
    cmd += "rm -rf %s-%s-*.txt %s'" % (host, port, 
            os.path.basename(args.MQfile.name))
    return cmd


if args.host_file:
    print_stamp("[+] mode: cluster")

    if args.mq_fix_file:
        print_stamp("[!] ignore local preprocess configuration file in cluster mode",
              file=sys.stderr)
        args.mq_fix_file = None

    if not args.work_dir:
        print_stamp("[!] missing option: -W", file=sys.stderr)
        sys.exit(1)

    if args.mq_fix_num == 0:
        print_stamp("[!] invalid option: -q %d" % args.mq_fix_num)
        sys.exit(1)

    hosts = args.host_file.readlines()

    procs = []
    if not args.execute:
        print_stamp("[+] compiling %s on the cluster..." % prog_name)
        for line in hosts:
            host, port = line.strip().split(':')
            cmd = format_compile_cmd(host, port, args)
            comp_proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE)
            procs.append(comp_proc)

        for proc in procs:
            proc.wait()
            if proc.returncode != 0:
                print_stamp("[!] %s failed, return code %d" %
                      (''.join(proc.args), proc.returncode), file=sys.stderr)
                sys.exit(1)

    if args.compile:
        sys.exit(0)

    procs = []
    print_stamp("[+] copying MQ challenge file to the cluster...")
    for line in hosts:
        host, port = line.strip().split(':')
        cmd = format_scp_cmd(host, port, args)
        scp = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
        procs.append(scp)

    for proc in procs:
        proc.wait()
        if proc.returncode != 0:
            print_stamp("[!] %s failed, return code %d" %
                  (''.join(proc.args), proc.returncode), file=sys.stderr)
            sys.exit(1)

    procs = []

    print_stamp("[+] generating preprocess configuration on the cluster...")
    total_proc_num = len(hosts) * args.gpu_num
    prep_num = 2 ** args.mq_fix_num
    chunk_size = int(prep_num / total_proc_num)
    remainder = prep_num - total_proc_num * chunk_size

    print_stamp("number of processes to create: %d\n"
                "\t\tnumber of preprocess configurations: 2^%d\n"
                "\t\tdistribution: %d x %d + %d x %d" % (
                    total_proc_num, args.mq_fix_num, chunk_size,
                    total_proc_num - remainder, chunk_size+1, remainder))

    launch_info = []
    start = 0
    for line, id in itertools.product(hosts, range(0, args.gpu_num)):
        if remainder > 0:
            cur_chunk_size = chunk_size + 1
            remainder -= 1
        else:
            cur_chunk_size = chunk_size

        host, port = line.strip().split(':')
        cmd = format_gen_prep_cmd(host, port, args, id, start, cur_chunk_size)
        gen_proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
        procs.append(gen_proc)
        launch_info.append({'host': host, 'port': port, 'id': id,
                            'start': start, 'chunk_size': cur_chunk_size})
        start += cur_chunk_size

    for proc in procs:
        proc.wait()
        if proc.returncode != 0:
            print_stamp("[!] %s failed, return code %d" %
                  (''.join(proc.args), proc.returncode), file=sys.stderr)
            sys.exit(1)

    procs = []
    print_stamp("[+] invoking %s on the cluster..." % prog_name)

    if not os.path.exists(args.cluster_stdout):
        print_stamp("creating directory: %s" % args.cluster_stdout)
        os.makedirs(args.cluster_stdout, exist_ok=True)

    if not os.path.exists(args.cluster_stderr):
        print_stamp("creating directory: %s" % args.cluster_stderr)
        os.makedirs(args.cluster_stderr, exist_ok=True)

    for cfg in launch_info:
        host = cfg['host']
        port = cfg['port']
        id = cfg['id']
        cmd = format_mq_cmd(host, port, args, id)
        mq_proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
        procs.append(mq_proc)

    for proc in procs:
        proc.wait()
        if proc.returncode != 0:
            print_stamp("[!] %s failed, return code %d" %
                  (''.join(proc.args), proc.returncode), file=sys.stderr)
            sys.exit(1)

    print_stamp("[+] starting to poll the cluster periodically")
    print_stamp("polling period: %d seconds" % args.poll_period)
    done = False
    found_sol = False
    while done == False and found_sol == False:
        time.sleep(args.poll_period)

        print_stamp("[+] pulling results from the cluster...")
        procs = []
        for cfg in launch_info:
            host = cfg['host']
            port = cfg['port']
            id = cfg['id']
            cmd = format_scp_stats_cmd(host, port, args, id)
            pull_proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
            procs.append(pull_proc)

        for proc in procs:
            proc.wait()
            if proc.returncode != 0:
                print_stamp("[!] %s failed, return code %d" %
                      (''.join(proc.args), proc.returncode), file=sys.stderr)
                if not args.ignore_fail:
                    sys.exit(1)
        
        print_stamp("[+] checking results from the cluster...")
        cmd = "grep 'solution found' %s/*" % args.cluster_stdout
        check_proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
        check_proc.wait()
        if check_proc.returncode == 0:
            print_stamp("[+] a solution is found")
            found_sol = True
        elif check_proc.returncode == 1:
            print_stamp("[-] no solution found yet")
        else:
            print_stamp("[!] error while checking results, return code %d" %
                        check_proc.returncode, file=sys.stderr)

        print_stamp("[+] checking the status of the cluster")
        done = True
        for cfg in launch_info:
            host = cfg['host']
            port = cfg['port']
            cmd = format_check_cmd(host, port, args)
            check_proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE)
            check_proc.wait()
            if check_proc.returncode == 0:
                done = False
                break
            elif check_proc.returncode != 1:
                print_stamp("[!] error while checking MQsolver status, return code %d" %
                            check_proc.returncode, file=sys.stderr)

    if not done:
        procs = []
        print_stamp("[+] killing the rest of the cluster...")
        for line in hosts:
            host, port = line.strip().split(':')
            for id in range(0, args.gpu_num):
                cmd = format_kill_cmd(host, port, args, id)
                kill_proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE)
                procs.append(kill_proc)
                
        for proc in procs:
            proc.wait()
            if proc.returncode != 0:
                print_stamp("[!] %s failed, return code %d" %
                      (''.join(proc.args), proc.returncode), file=sys.stderr)

        procs = []

    if not found_sol:
        print_stamp("[-] solution not found")

    if args.email_notice:
        print_stamp("[+] sending notification to %s..." % args.email_notice)
        cmd = "mutt -s '%s finishes' %s < /dev/null" % (prog_name,
                args.email_notice)
        mail_proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
        mail_proc.wait()
        if mail_proc.returncode != 0:
            print_stamp('[!] failed to send notification to %s' % args.email_notice)

    if args.no_clean:
        sys.exit(0)

    print_stamp("[+] cleaning...")
    for line in hosts:
        host, port = line.strip().split(':')
        cmd = format_clean_cmd(host, port, args)
        clean_proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
        procs.append(clean_proc)
            
    for proc in procs:
        proc.wait()
        if proc.returncode != 0:
            print_stamp("[!] %s failed, return code %d" %
                  (''.join(proc.args), proc.returncode), file=sys.stderr)

else:
    print_stamp("[+] mode: normal")

    if not args.execute:
        extra_cflags = ""
        if args.count_linear_dep:
            extra_cflags += " -DGC_DEPC_LSYS"

        if os.path.exists(build_dir):
            print_stamp("cleaning directory: %s" % build_dir)
            shutil.rmtree(build_dir)

        print_stamp("creating directory: %s" % build_dir)
        os.makedirs(build_dir)

        print_stamp("[+] generating C code according to the parameters...")
        cmd = "cd %s && cmake -DKEEP_VAR_NUM=%d -DCMAKE_C_FLAGS=%s .." % (
                build_dir, args.keep_var_num, extra_cflags)
        with subprocess.Popen(cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE) as gen_meta:
            gen_meta.wait()
            check_subproc(gen_meta)

        print_stamp("[+] compiling %s..." % prog_name)

        compile_stdout = subprocess.PIPE
        compile_stderr = subprocess.PIPE
        if args.verbose_compile:
            compile_stdout = sys.stdout
            compile_stderr = sys.stderr

        cmd = "cd %s && make" % build_dir

        with subprocess.Popen(cmd, shell=True, stdout=compile_stdout,
                              stderr=compile_stderr) as compile_proc:
            compile_proc.wait()
            check_subproc(compile_proc)

    if args.compile:
        sys.exit(0)

    if args.mq_fix_num != 0:
        args.mq_fix_file = "prep%s.txt" % args.mq_fix_num
        print_stamp("[+] generating preprocess configuration file: %s..." %
              args.mq_fix_file)
        with subprocess.Popen(["bin/gen_efix.py", args.mq_fix_file,
                               str(args.mq_fix_num),
                               "0",
                               str(2 ** args.mq_fix_num)],
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE) as gen_prepcfg:
            gen_prepcfg.wait()
            check_subproc(gen_prepcfg)

    print_stamp("[+] invoking %s..." % prog_name)

    run_stdout = sys.stdout
    run_stderr = sys.stderr

    if args.output_file:
        print_stamp("[+] piping output to file: %s..." % args.output_file.name)
        run_stdout = args.output_file
        run_stderr = args.output_file

    with subprocess.Popen(['./%s/mqsolver' % build_dir,
                           '--challenge=%s' % args.MQfile.name,
                           '--algorithm=%s' % args.algorithm,
                           '--macaulay_deg=%d' % args.degree,
                           '--verbose' if args.verbose else '',
                           '--keep_var=%d' % args.keep_var_num,
                           '--kf_var=%d' % args.thread_fix_num,
                           '--mq_fix=%s' % args.mq_fix_file if args.mq_fix_file
                                else '',
                           '--mac_keq=%d' % args.subsys_candidate_keep_num,
                           '--seed=%d%' % args.seed if args.seed else ''
                           '--mailbox_size=%d' % args.mailbox_num,
                           '--mac_stats' if args.macaulay_stats else '',
                           '--sub_fvar=%d' % args.subsys_fix_num,
                           '--dev_id=%d' % args.gpu_dev_id,
                           '--rmac_cpu' if args.cpu_reduction else '',
                           '--mq_ratio=%f' % args.mq_max_ratio,
                           '--thread_num=%d' % args.cpu_thread_num if
                                args.cpu_thread_num else ''
                           '--sub_keq=%d' % args.subsys_eq_keep_num,
                           '--init_gj=%s' % str(args.init_gj).lower(),
                           '--mq_ext_file=%s' % args.mq_ext_file if args.mq_ext_file
                                else '',
                           '--mode=%s' % args.mode,
                           ],
                          stdout=run_stdout, stderr=run_stderr) as run_proc:

        run_proc.wait()

        if args.no_clean:
            sys.exit(run_proc.returncode)

        print_stamp("[+] cleaning...")
        cmd = "rm -rf %s %s" % (build_dir, args.mq_fix_file if
               args.mq_fix_file else '')
        with subprocess.Popen(cmd, shell=True) as clean_proc:
            clean_proc.wait()

        sys.exit(run_proc.returncode)
