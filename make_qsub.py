#!/usr/bin/env python

# python 3 compatibility
from __future__ import print_function
try:
    import builtins
except ImportError:
    import __builtin__ as builtins

# other modules
import os
import stat
import sys
import math
import json

# important paths
ROOT = "/home/users/vrastil/GIT/FastSim/"

SCRIPT_DIR = ROOT + "jobs/batch_scripts/"
OUT_DIR = ROOT + "jobs/output/"
LOG_DIR = ROOT + "jobs/logs/"
INPUT_CFG = ROOT + "jobs/input/generic_input.cfg"
SINGULARITY_IMG = ROOT + "jobs/FastSim.sif"

class Job_Param(object):
    def __init__(self, app, mem, cpus, n_cpus):
        mem, h, m = get_safe_mem_wall_time(mem, cpus, n_cpus)
        self.app = app
        self.mem = mem
        self.n_cpus = n_cpus
        self.wall_time_h = h
        self.wall_time_m = m
        self.sim_opt = {
                        "app" : "",
                        "size" : "",
                        "redshift" : "",
                        "mlt" : ""
                        }

    def add_sim_opt(self, sim_opt, key):
        self.sim_opt[key] += sim_opt

    def add_std_opt(self, sim_param):
        self.add_sim_opt("--mesh_num %i " % sim_param.Nm, "size")
        self.add_sim_opt("--mesh_num_pwr %i " % sim_param.NM, "size")
        self.add_sim_opt("--par_num %i " % sim_param.Np, "size")
        self.add_sim_opt("--box_size %f " % sim_param.box, "size")
        self.add_sim_opt("--redshift %f " % sim_param.z, "redshift")
        self.add_sim_opt("--redshift_0 %f " % sim_param.z0, "redshift")
        self.add_sim_opt("--time_step %f " % sim_param.da, "redshift")
        self.add_sim_opt("--print_every %i " % sim_param.print_every, "redshift")
        self.add_sim_opt("--print_vel_pwr %i " % sim_param.print_vel_pwr, "redshift")
        self.add_sim_opt("--mlt_runs %i " % sim_param.mlt_runs, "mlt")
        self.add_sim_opt("--pair %i " % sim_param.pair, "mlt")
        self.add_sim_opt("--smoothing_k %f " % sim_param.smoothing_k, "app")


class Sim_Param(object):
    def __init__(self, Nm=0, NM=0, Np=0, box=0, z=0, z0=0, da=0, print_every=0, smoothing_k=0, print_vel_pwr=0):
        self.Nm = Nm
        self.NM = NM
        self.Np = Np
        self.num_m = pow(Nm, 3)
        self.num_M = pow(NM, 3)
        self.num_p = pow(Np, 3)
        self.box = box
        self.z = z
        self.z0 = z0
        self.da = da
        self.print_every = print_every
        self.n_steps = int(math.ceil((1. / (1. + z0) - 1. / (1. + z)) / da))
        self.n_print = int(
            math.ceil(self.n_steps / float(print_every))) + 1 if print_every else 0
        self.rs = 0
        self.mlt_runs = 1
        self.pair = False
        self.chi_phi = 0
        self.chi_n = 0
        self.comp_chi_lin = 0
        self.smoothing_k = smoothing_k
        self.print_vel_pwr = print_vel_pwr


def memory_base(sim_param):
    mem = 0
    mem += sim_param.num_m * 3 * 8  # app_field
    mem += sim_param.num_M * 3 * 8  # power_aux
    mem += sim_param.num_p * 3 * 2 * 8  # particles with velocities
    return mem


def memory_za(sim_param):
    return memory_base(sim_param) / float(1024 * 1024 * 1024)  # convert to GB


def memory_aa(sim_param):
    mem = memory_base(sim_param)
    mem += sim_param.num_M * 8  # expotential
    return mem / float(1024 * 1024 * 1024)  # convert to GB


def memory_fp_pp(sim_param):
    mem = memory_base(sim_param)
    mem += sim_param.num_p * 4  # linked list
    M = (int)(sim_param.Nm / float(sim_param.rs))
    mem += pow(M, 3) * 4  # linked list
    return mem / float(1024 * 1024 * 1024)  # convert to GB

def memory_chi(sim_param):
    mem = memory_base(sim_param)
    mem += sim_param.num_m * 8 * 3  # chi_force
    mem += 8*(sim_param.num_m - 1) / 7 * 8 * 1  # drho
    mem += 8*(sim_param.num_m - 1) / 7 * 8 * 3  # chi_solver
    mem *= 1.3 # ??? jobs are killed due to the memory exceed limit
    return mem / float(1024 * 1024 * 1024)  # convert to GB

def get_std_input():
    Nm = int(input("Enter number of potential mesh points per dimenson: "))
    NM = int(input("Enter number of analysis mesh points per dimenson: "))
    Np = int(input("Enter number of particles per dimenson: "))
    box = float(input("Enter size of the simulation box: "))
    z = float(input("Enter initial redshift of the simulation: "))
    z0 = float(input("Enter final redshift of the simulation: "))
    da = float(input("Enter value of time-step: "))
    print_every = int(input("Enter how often there will be printing: "))

    sim_param = Sim_Param(Nm=Nm, NM=NM, Np=Np, box=box, z=z, z0=z0, da=da, print_every=print_every)
    # sim_param.rs = float(input("Enter value of short-range cutof (FP_pp): "))
    sim_param.mlt_runs = int(input("Enter number of runs: "))
    sim_param.pair = bool(input("Run pair of simulations? "))
    return sim_param

def get_chi_input(sim_param):
    sim_param.chi_phi = float(input("Enter value of screening potential (CHI): "))
    sim_param.chi_n = float(input("Enter value chameleon power-law potential exponent (CHI): "))
    sim_param.comp_chi_lin = bool(input("Run linear solver of chameleon? "))
    sim_param.chi_type = int(input("Which type of approximation is to be used for chamaleon?\n\tFrozen-potential (0)\n\tFrozen-flow (1)\n "))

def get_input(sim_param=None, with_chi=False):
    """ Create simulation parameters. If given already created sim_param, get only chameleon input. """
    if sim_param is None:
        sim_param = get_std_input()
    if with_chi:
        get_chi_input(sim_param)
    return sim_param

def get_param_from_json(data, default=None, with_chi=False):
    # load default values, DO NOT CHANGE IT
    if default is not None:
        basic = default["basic"].copy()
        mlt = default["mlt"].copy()
        chi = default["chi"].copy()
    else:
        basic, mlt, chi = {}, {}, {}
    
    # override specified parameters
    for group in ["basic", "mlt", "chi"]:
        if group in data:
            for key, value in data[group].items():
                locals()[group][key] = value

    sim_param = Sim_Param(**basic)
    sim_param.mlt_runs = mlt["mlt_runs"]
    sim_param.pair = mlt["pair"]
    if with_chi:
        sim_param.chi_phi = chi["phi"]
        sim_param.chi_n = chi["n"]
        sim_param.comp_chi_lin = chi["lin"]
        sim_param.chi_type = chi["type"]
    return sim_param

def cpu_base(sim_param, prep_Nm=1, prep_Np=1, integ_np_nsteps=1, print_np=1, print_NM=1):
    cpus = 0
    cpus += prep_Nm * pow(sim_param.Nm / 128., 3)  # preparation
    cpus += prep_Np * pow(sim_param.Np / 128., 3)  # preparation
    cpus += integ_np_nsteps * sim_param.n_steps * \
        pow(sim_param.Np / 128., 3)  # integration
    cpus += print_np * sim_param.n_print * \
        pow(sim_param.Np / 128., 3)  # printing
    cpus += print_NM * sim_param.n_print * \
        pow(sim_param.NM / 128., 3)  # printing
    return cpus


def cpu_pp(sim_param, prep_Nm=1, prep_Np=1, integ_np_nsteps=1, print_np=1, print_NM=1, integ_np_nsteps_extra=1):
    cpus = cpu_base(sim_param, prep_Nm, prep_Np,
                    integ_np_nsteps, print_np, print_NM)
    cpus += integ_np_nsteps_extra * sim_param.n_steps * \
        pow(sim_param.Np / 128., 3) * pow(sim_param.rs / 2.7, 3)
    return cpus

def cpu_chi(sim_param, prep_Nm=1, prep_Np=1, integ_np_nsteps=1, print_np=1, print_NM=1, integ_np_nsteps_extra=1):
    cpus = cpu_base(sim_param, prep_Nm, prep_Np,
                    integ_np_nsteps, print_np, print_NM)
    cpus += integ_np_nsteps_extra * sim_param.n_steps * \
        pow(sim_param.Nm / 128., 3)
    return cpus

def convert_s2hms(seconds):
    h = int(seconds) // (60 * 60)
    m = (int(seconds) // 60) % 60
    s = int(seconds) % 60
    return h, m, s


def time_2string(h, m, s):
    time = "%i:" % h
    if m < 10:
        time += "0%i:" % m
    else:
        time += "%i:" % m
    if s < 10:
        time += "0%i" % s
    else:
        time += "%i" % s
    return time

def get_eff_n_cpus(n_cpus):
    """ Effectiveness drops with square root of n_cpus,
        such that with 16 corer it is 75%. """
    eff = 13./12 - 1./12*math.sqrt(n_cpus)
    return eff*n_cpus

def print_n_cpus_t(cpus, n_cpus):
    eff_n_cpu = get_eff_n_cpus(n_cpus)
    h, m, s = convert_s2hms(cpus / eff_n_cpu)
    time = time_2string(h, m, s)
    print("\tOn %i cores the simulation will take %s of Wall time." % (n_cpus, time))


def get_n_cpus(cpus, mem, approx_str):
    h, m, s = convert_s2hms(cpus)
    time = time_2string(h, m, s)
    print("\n%s will need approximately %.1f GB of memory and shouldn`t take more than %s of CPU time." % (approx_str, mem, time))

    print_n_cpus_t(cpus, 4)
    if h > 1:
        print_n_cpus_t(cpus, 8)
    if h > 2:
        print_n_cpus_t(cpus, 16)
    if h > 24:
        print_n_cpus_t(cpus, 32)
        print_n_cpus_t(cpus, 64)
    
    return int(input("Enter number of cores: "))


def get_safe_mem_wall_time(mem, cpus, n_cpus):
    eff_n_cpu = get_eff_n_cpus(n_cpus)
    mem = int(math.ceil((mem + 0.3) * 1.03))  # extra 0.3 GB and 3%
    h, m, _ = convert_s2hms(cpus * 3.0 / eff_n_cpu)  # extra 200% of Wall time
    if (h + m) < 1:
        m = 1
    return mem, h, m  # seconds are ambiguous


def make_sbatch_koios(job):
    sbatch = "#!/bin/bash\n"
    sbatch += ("#SBATCH --nodes=1\n"
               "#SBATCH --cpus-per-task=32\n")
    sbatch += "#SBATCH --mem=%iG\n" % job.mem
    if job.wall_time_m < 10:
        sbatch += "#SBATCH --time=%i:0%i:00\n" % (
            job.wall_time_h, job.wall_time_m)
    else:
        sbatch += "#SBATCH --time=%i:%i:00\n" % (
            job.wall_time_h, job.wall_time_m)
    sbatch += "#SBATCH --job-name=%s_fastsim\n" % job.app
    sbatch += "#SBATCH --output=%s/%%x_%%j.log\n" % LOG_DIR
    sbatch += "#SBATCH --error=%s/%%x_%%j.log\n" % LOG_DIR
    sbatch += "#SBATCH --partition=long\n"
    sbatch += "#SBATCH --exclusive\n"
    ##############
    # RUN SCRIPT #
    ##############
    sbatch += "\n# preparation\n"
    sbatch += "export OMP_NUM_THREADS=32\n"
    sbatch += "\n# input/output\n"
    sbatch += "OUT_DIR='%s'\n" % OUT_DIR
    sbatch += "SINGULARITY_IMG='%s'\n" % SINGULARITY_IMG
    sbatch += "INPUT_CFG='%s'\n" % INPUT_CFG
    sbatch += "\n# parameters\n"
    for key, value in job.sim_opt.items():
        sbatch += "%s='%s'\n" % (key.upper(), value)
    sbatch += 'GENERIC="-c $INPUT_CFG --out_dir=$OUT_DIR"\n'
    sbatch += "\n# run\n"
    sbatch += "singularity exec $SINGULARITY_IMG FastSim $GENERIC"
    for key in job.sim_opt.keys():
        sbatch += " $%s" % key.upper()
    return sbatch + "\n"

def save_job_file(job_script, job_file):
    with open(job_file, 'w') as a_file:
        a_file.write(job_script)
    print("Job parameters saved to '%s'" % job_file)

def make_job_scripts(job, app):
    sbatch_koios = make_sbatch_koios(job)

    # save_job_file(qsub_meta, "Metacentrum/%s_qsub.pbs" % app)
    save_job_file(sbatch_koios, "%s/%s_sbatch.sh" % (SCRIPT_DIR, app))


def make_submit(job_files, cmd='sbatch'):
    # bash for-loop start
    submit = ("#!/bin/bash\n"
              "NUM_ALL=$1\n"
              "for NUM in `seq $NUM_ALL`; do\n")
    
    # bash for-loop 
    for job_file in job_files:
        submit += "    %s %s\n" % (cmd, job_file)

    # bash for-loop end
    submit += "done\n"
    return submit

def create_mlt_submit(job_files, sh_file="%s/mlt_CHI.sh" % SCRIPT_DIR, cmd='sbatch'):
    # create bash file
    submit = make_submit(job_files, cmd=cmd)
    save_job_file(submit, sh_file)

    # make it executable
    st = os.stat(sh_file)
    os.chmod(sh_file, st.st_mode | stat.S_IEXEC)

# already quite safe values
PREP_PAR = 0.4
PRINT_PAR = 3.0
PRINT_NM = 1.7

def paired_sim(sim_param):
    if sim_param.pair:
        return 2
    else:
        return 1

def cpu_mlt(sim_param):
    return paired_sim(sim_param) * sim_param.mlt_runs

def qsub_ZA(sim_param, app='ZA'):
    cpu_param = {
        'prep_Nm' : 2.0, 
        'prep_Np' : PREP_PAR,
        'integ_np_nsteps' : 0.25,
        'print_np' : PRINT_PAR,
        'print_NM' : PRINT_NM
    }

    cpus = cpu_mlt(sim_param) * cpu_base(sim_param, **cpu_param)
    mem = memory_za(sim_param)
    n_cpus = 32 # get_n_cpus(cpus, mem, "Zel`dovich approximation")
    ZA = Job_Param('ZA', mem, cpus, n_cpus)
    ZA.add_std_opt(sim_param)
    ZA.add_sim_opt("--comp_ZA 1 ", "app")
    make_job_scripts(ZA, app)

def qsub_TZA(sim_param, app='TZA'):
    cpu_param = {
        'prep_Nm' : 2.0, 
        'prep_Np' : PREP_PAR,
        'integ_np_nsteps' : 0.25,
        'print_np' : PRINT_PAR,
        'print_NM' : PRINT_NM
    }

    cpus = cpu_mlt(sim_param) * cpu_base(sim_param, **cpu_param)
    mem = memory_za(sim_param)
    n_cpus = 32 # get_n_cpus(cpus, mem, "Truncated Zel`dovich approximation")
    TZA = Job_Param('TZA', mem, cpus, n_cpus)
    TZA.add_std_opt(sim_param)
    TZA.add_sim_opt("--comp_TZA 1 ", "app")
    make_job_scripts(TZA, app)

def qsub_FF(sim_param, app='FF'):
    cpu_param = {
        'prep_Nm' : 2.0, 
        'prep_Np' : PREP_PAR,
        'integ_np_nsteps' : 3.0,
        'print_np' : PRINT_PAR,
        'print_NM' : PRINT_NM
    }

    cpus = cpu_mlt(sim_param) * cpu_base(sim_param, **cpu_param)
    mem = memory_za(sim_param)
    n_cpus = 32 # get_n_cpus(cpus, mem, "Frozen-flow approximation")
    FF = Job_Param('FF', mem, cpus, n_cpus)
    FF.add_std_opt(sim_param)
    FF.add_sim_opt("--comp_FF 1 ", "app")
    make_job_scripts(FF, app)

def qsub_FP(sim_param, app='FP'):
    cpu_param = {
        'prep_Nm' : 9.0, 
        'prep_Np' : PREP_PAR,
        'integ_np_nsteps' : 3.5,
        'print_np' : PRINT_PAR,
        'print_NM' : PRINT_NM
    }

    cpus = cpu_mlt(sim_param) * cpu_base(sim_param, **cpu_param)
    mem = memory_za(sim_param)
    n_cpus = 32 # get_n_cpus(cpus, mem, "Frozen-potential approximation")
    FP = Job_Param('FP', mem, cpus, n_cpus)
    FP.add_std_opt(sim_param)
    FP.add_sim_opt("--comp_FP 1 ", "app")
    make_job_scripts(FP, app)

def qsub_FP_pp(sim_param, app='FP_pp'):
    cpu_param = {
        'prep_Nm' : 20.0, 
        'prep_Np' : PREP_PAR,
        'integ_np_nsteps' : 50,
        'print_np' : PRINT_PAR,
        'print_NM' : PRINT_NM,
        'integ_np_nsteps_extra' : 0
    }

    cpus = cpu_mlt(sim_param) * cpu_pp(sim_param, **cpu_param)
    mem = memory_fp_pp(sim_param)
    n_cpus = 32 # get_n_cpus(cpus, mem, "Frozen-potential particle-particle approximation")
    FP_pp = Job_Param('FP_pp', mem, cpus, n_cpus)
    FP_pp.add_std_opt(sim_param)
    FP_pp.add_sim_opt("--cut_radius %f " % sim_param.rs, "app")
    FP_pp.add_sim_opt("--comp_FP_pp 1 ", "app")
    make_job_scripts(FP_pp, app)

def qsub_CHI_base(sim_param, with_chi=True):
    # get chameleon simulation parameters
    sim_param = get_input(sim_param=sim_param, with_chi=with_chi)

    # get required memory and wall time
    cpu_param = {
        'prep_Nm' : 9.0, 
        'prep_Np' : PREP_PAR,
        'integ_np_nsteps' : 3.0,
        'print_np' : PRINT_PAR,
        'print_NM' : PRINT_NM,
        'integ_np_nsteps_extra' : 150
    }
    cpus = cpu_mlt(sim_param) * cpu_chi(sim_param, **cpu_param)
    mem = memory_chi(sim_param)
    n_cpus = 32 # get_n_cpus(cpus, mem, "Chameleon gravity approximation")

    # creat chameleon job
    CHI = Job_Param('CHI', mem, cpus, n_cpus)
    CHI.add_std_opt(sim_param)
    if sim_param.chi_type == 0:
        CHI.add_sim_opt("--comp_chi 1 ", "app")
    elif sim_param.chi_type == 1:
        CHI.add_sim_opt("--comp_chi_ff 1 ", "app")
    CHI.add_sim_opt("--chi_n %f " % sim_param.chi_n, "app")
    CHI.add_sim_opt("--chi_phi %E " % sim_param.chi_phi, "app")
    CHI.add_sim_opt("--comp_chi_lin %i " % sim_param.comp_chi_lin, "app")

    return CHI

def qsub_CHI(sim_param, app='CHI', with_chi=False):
    CHI = qsub_CHI_base(sim_param, with_chi=with_chi)
    make_job_scripts(CHI, app)

def qsub_CHI_mlt(sim_param):
    print("Creating multiple job files for Chameleon gravity.")
    Nchi = int(input("Enter number of sets of parameters: "))
    job_files = []
    for i in range(Nchi):
        app = 'CHI_%i' % i
        qsub_CHI(sim_param, app=app, with_chi=True)
        job_files.append("%s/%s_sbatch.sh" % (SCRIPT_DIR, app))

    create_mlt_submit(job_files)

def qsub_json(data):
    default = data.pop("default", None)
    job_files = []
    for app, param_list in data.items():
        try:
            qsub_fce = getattr(sys.modules[__name__], 'qsub_%s' % app)
            for i, param in enumerate(param_list):
                sim_param = get_param_from_json(param, default=default, with_chi=(app=='CHI'))
                name = "%s_%i" % (app, i)
                qsub_fce(sim_param, app=name)
                job_files.append("%s/%s_sbatch.sh" % (SCRIPT_DIR, name))
        except AttributeError:
            print("Wrong format of json file!")
    create_mlt_submit(job_files, sh_file="%s/all_runs.sh" % SCRIPT_DIR)


if __name__ == "__main__":
    # standard creation of all job files
    if len(sys.argv) == 1:
        sim_param = get_input()
        qsub_ZA(sim_param)
        qsub_TZA(sim_param)
        qsub_FF(sim_param)
        qsub_FP(sim_param)
        # qsub_FP_pp(sim_param)
        qsub_CHI(sim_param, with_chi=True)

    # we passed json file with all the info
    elif sys.argv[1].endswith('.json'):
        with open(sys.argv[1]) as data_file:
            data = json.loads(data_file.read())
            qsub_json(data)

    # create only one job file
    else:
        try:
            qsub_fce = getattr(sys.modules[__name__], 'qsub_%s' % sys.argv[1])
            sim_param = get_input()
            qsub_fce(sim_param)
        except AttributeError:
            print("Wrong arguments!")
