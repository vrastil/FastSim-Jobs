#!/usr/bin/env python
import math
import os
import stat


class Job_Param(object):
    def __init__(self, app, mem, cpus, n_cpus):
        mem, h, m = get_safe_mem_wall_time(mem, cpus, n_cpus)
        self.app = app
        self.mem = mem
        self.n_cpus = n_cpus
        self.wall_time_h = h
        self.wall_time_m = m
        self.sim_opt = ""

    def add_sim_opt(self, sim_opt):
        self.sim_opt += sim_opt

    def add_std_opt(self, sim_param):
        self.add_sim_opt("--mesh_num %i " % sim_param.Nm)
        self.add_sim_opt("--mesh_num_pwr %i " % sim_param.NM)
        self.add_sim_opt("--par_num %i " % sim_param.Np)
        self.add_sim_opt("--box_size %f " % sim_param.box)
        self.add_sim_opt("--redshift %f " % sim_param.z)
        self.add_sim_opt("--redshift_0 %f " % sim_param.z0)
        self.add_sim_opt("--time_step %f " % sim_param.da)
        self.add_sim_opt("--print_every %i " % sim_param.print_every)
        self.add_sim_opt("--mlt_runs %i " % sim_param.mlt_runs)


class Sim_Param(object):
    def __init__(self, Nm, NM, Np, box, z, z0, da, print_every):
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
        self.chi_phi = 0
        self.chi_n = 0


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
    return mem / float(1024 * 1024 * 1024)  # convert to GB

def get_input():
    Nm = input("Enter number of potential mesh points per dimenson: ")
    NM = input("Enter number of analysis mesh points per dimenson: ")
    Np = input("Enter number of particles per dimenson: ")
    box = input("Enter size of the simulation box: ")
    z = input("Enter initial redshift of the simulation: ")
    z0 = input("Enter final redshift of the simulation: ")
    da = input("Enter value of time-step: ")
    print_every = input("Enter how often there will be printing: ")

    sim_param = Sim_Param(Nm, NM, Np, box, z, z0, da, print_every)
    sim_param.rs = input("Enter value of short-range cutof (FP_pp): ")
    sim_param.mlt_runs = input("Enter number of runs: ")
    sim_param.chi_phi = input("Enter value of screening potential (CHI): ")
    sim_param.chi_n = input("Enter value chameleon power-law potential exponent (CHI): ")
    return sim_param


def cpu_base(sim_param, prep_Nm, prep_Np, integ_np_nsteps, print_np, print_NM):
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


def cpu_pp(sim_param, prep_Nm, prep_Np, integ_np_nsteps, print_np, print_NM, integ_np_nsteps_extra):
    cpus = cpu_base(sim_param, prep_Nm, prep_Np,
                    integ_np_nsteps, print_np, print_NM)
    cpus += integ_np_nsteps_extra * sim_param.n_steps * \
        pow(sim_param.Np / 128., 3) * pow(sim_param.rs / 2.7, 3)
    return cpus

def cpu_chi(sim_param, prep_Nm, prep_Np, integ_np_nsteps, print_np, print_NM, integ_np_nsteps_extra):
    cpus = cpu_base(sim_param, prep_Nm, prep_Np,
                    integ_np_nsteps, print_np, print_NM)
    cpus += integ_np_nsteps_extra * sim_param.n_steps * \
        pow(sim_param.Nm / 128., 3)
    return cpus

def convert_s2hms(seconds):
    h = int(seconds) / (60 * 60)
    m = (int(seconds) / 60) % 60
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
    print "\tOn %i cores the simulation will take %s of Wall time." % (n_cpus, time)


def get_n_cpus(cpus, mem, approx_str):
    h, m, s = convert_s2hms(cpus)
    time = time_2string(h, m, s)
    print "\n%s will need approximately %.1f GB of memory and shouldn`t take more than %s of CPU time." % (approx_str, mem, time)

    print_n_cpus_t(cpus, 4)
    if h > 1:
        print_n_cpus_t(cpus, 8)
    if h > 2:
        print_n_cpus_t(cpus, 16)
    if h > 24:
        print_n_cpus_t(cpus, 32)
        print_n_cpus_t(cpus, 64)
    
    return input("Enter number of cores: ")


def get_safe_mem_wall_time(mem, cpus, n_cpus):
    eff_n_cpu = get_eff_n_cpus(n_cpus)
    mem = int(math.ceil((mem + 0.3) * 1.03))  # extra 0.3 GB and 3%
    h, m, _ = convert_s2hms(cpus * 1.5 / eff_n_cpu)  # extra 50% of Wall time
    if (h + m) < 1:
        m = 1
    return mem, h, m  # seconds are ambiguous


def make_qsub_meta(job):
    qsub = "#!/bin/bash\n"
    qsub += "#PBS -l select=1:ncpus=%i:mem=%igb:scratch_local=400mb\n" % (job.n_cpus, job.mem)
    if job.wall_time_m < 10:
        qsub += "#PBS -l walltime=%i:0%i:00\n" % (
            job.wall_time_h, job.wall_time_m)
    else:
        qsub += "#PBS -l walltime=%i:%i:00\n" % (
            job.wall_time_h, job.wall_time_m)
    qsub += "#PBS -N %s_fastsim\n" % job.app
    qsub += ("#PBS -j oe\n"
             "#PBS -o logs/\n"
             "#PBS -e logs/\n\n\n"
             ##############
             # RUN SCRIPT #
             ##############
             "trap 'clean_scratch' TERM EXIT\n"
             "export MYDIR=/storage/brno2/home/vrastilm\n"
             "source /software/modules/init\n"	
             "source $MYDIR/.profile\n"	
             "export OMP_NUM_THREADS=$PBS_NUM_PPN\n"
             "cd $SCRATCHDIR\n"
             "cp -r $MYDIR/Adhesion-Approximation/include ./\n"
             "cp -r $MYDIR/Adhesion-Approximation/src ./\n"	
             "cp $MYDIR/Adhesion-Approximation/Makefile ./\n"	
             "make clean\n"	
             "make -j $PBS_NUM_PPN\n\n"	
             )
    qsub += "time ./adh_app -c $MYDIR/Adhesion-Approximation/input/generic_input.cfg %s\n" % job.sim_opt
    return qsub

def make_sbatch_koios(job):
    MYDIR = "/home/users/vrastil/GIT/FastSim-Jobs"
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
    sbatch += "#SBATCH --output=%s/logs/%%x_%%j.log\n" % MYDIR
    sbatch += "#SBATCH --error=%s/logs/%%x_%%j_err.log\n" % MYDIR
    sbatch += "#SBATCH --partition=long\n"
    ##############
    # RUN SCRIPT #
    ##############
    sbatch += "export OMP_NUM_THREADS=32\n"
    sbatch += "cd %s/../FastSim-Container/debian/\n" % MYDIR
    sbatch += "singularity exec -B %s:/data fastsim.simg fastsim -c /data/input/generic_input.cfg %s --out_dir=/data/output/\n" % (MYDIR, job.sim_opt)
    return sbatch

def save_job_file(job_script, job_file):
    with open(job_file, 'w') as file:
        file.write(job_script)

def make_job_scripts(job, app):
    qsub_meta = make_qsub_meta(job)
    sbatch_koios = make_sbatch_koios(job)

    save_job_file(qsub_meta, "Metacentrum/%s_qsub.pbs" % app)
    save_job_file(sbatch_koios, "KOIOS/%s_sbatch.sh" % app)


# def make_submit():
#     submit = ("#!/bin/bash\n"
#               "NUM_ALL=$1\n"
#               "for NUM in `seq $NUM_ALL`; do\n"
#               "    qsub scripts/ZA_qsub.pbs\n"
#               "    qsub scripts/FF_qsub.pbs\n"
#               "    qsub scripts/FP_qsub.pbs\n"
#               "    qsub scripts/FP_pp_qsub.pbs\n"
#               "done\n")
#     return submit

def make_stack_qsub():
    qsub = ("#!/bin/bash\n"
            "#PBS -l select=1:ncpus=1:mem=400mb\n"
            "#PBS -l walltime=0:30:00\n"
            "#PBS -j oe\n"
            "#PBS -N cosmo_stack\n"
            "#PBS -o logs/\n"
            "#PBS -e logs/\n\n\n"
            "source /software/modules/init\n"
            "module add python27-modules-intel\n"
            "export PYTHONPATH=$PYTHONPATH:/storage/brno2/home/vrastilm/Adhesion-Approximation:/storage/brno2/home/vrastilm/CCL\n"
            "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/storage/brno2/home/vrastilm/local/lib\n\n"
            "COM+=\"from simpy.stacking import stack_all; \"\n"
            "COM+=\"out_dir='/storage/brno2/home/vrastilm/Adhesion-Approximation/output/'; \"\n"
            "COM+=\"stack_all(out_dir)\"\n\n"
            "python -c \"$COM\"\n"
            )
    return qsub

# already quite safe values
PREP_PAR = 0.4
PRINT_PAR = 3.0
PRINT_NM = 1.7


def qsub_ZA(sim_param):
    cpu_param = 2.0, PREP_PAR, 0.25, PRINT_PAR, PRINT_NM
    cpus = sim_param.mlt_runs * cpu_base(sim_param, *cpu_param)
    mem = memory_za(sim_param)
    n_cpus = get_n_cpus(cpus, mem, "Zel`dovich approximation")
    ZA = Job_Param('ZA', mem, cpus, n_cpus)
    ZA.add_std_opt(sim_param)
    ZA.add_sim_opt("--comp_ZA 1 ")
    make_job_scripts(ZA, 'ZA')

    # save_to_qsub(make_qsub(ZA), "scripts/ZA_qsub.pbs")


def qsub_FF(sim_param):
    cpu_param = 1.7, PREP_PAR, 2.0, PRINT_PAR, PRINT_NM
    cpus = sim_param.mlt_runs * cpu_base(sim_param, *cpu_param)
    mem = memory_za(sim_param)
    n_cpus = get_n_cpus(cpus, mem, "Frozen-flow approximation")
    FF = Job_Param('FF', mem, cpus, n_cpus)
    FF.add_std_opt(sim_param)
    FF.add_sim_opt("--comp_FF 1 ")
    make_job_scripts(FF, 'FF')

    # save_to_qsub(make_qsub(FF), "scripts/FF_qsub.pbs")


def qsub_FP(sim_param):
    cpu_param = 7.0, PREP_PAR, 2.1, PRINT_PAR, PRINT_NM
    cpus = sim_param.mlt_runs * cpu_base(sim_param, *cpu_param)
    mem = memory_za(sim_param)
    n_cpus = get_n_cpus(cpus, mem, "Frozen-potential approximation")
    FP = Job_Param('FP', mem, cpus, n_cpus)
    FP.add_std_opt(sim_param)
    FP.add_sim_opt("--comp_FP 1 ")
    make_job_scripts(FP, 'FP')

    # save_to_qsub(make_qsub(FP), "scripts/FP_qsub.pbs")


def qsub_FP_pp(sim_param):
    cpu_param = 16.0, PREP_PAR, 44.0, PRINT_PAR, PRINT_NM, 0
    cpus = sim_param.mlt_runs * cpu_pp(sim_param, *cpu_param)
    mem = memory_fp_pp(sim_param)
    n_cpus = get_n_cpus(
        cpus, mem, "Frozen-potential particle-particle approximation")
    FP_pp = Job_Param('FP_pp', mem, cpus, n_cpus)
    FP_pp.add_std_opt(sim_param)
    FP_pp.add_sim_opt("--cut_radius %f " % sim_param.rs)
    FP_pp.add_sim_opt("--comp_FP_pp 1 ")
    make_job_scripts(FP_pp, 'FP_pp')

    # save_to_qsub(make_qsub(FP_pp), "scripts/FP_pp_qsub.pbs")



def qsub_CHI(sim_param):
    cpu_param = 7.0, PREP_PAR, 2.1, PRINT_PAR, PRINT_NM, 120
    cpus = sim_param.mlt_runs * cpu_chi(sim_param, *cpu_param)
    mem = memory_chi(sim_param)
    n_cpus = get_n_cpus(cpus, mem, "Chameleon gravity approximation")
    CHI = Job_Param('CHI', mem, cpus, n_cpus)
    CHI.add_std_opt(sim_param)
    CHI.add_sim_opt("--comp_chi 1 ")
    CHI.add_sim_opt("--chi_n %f " % sim_param.chi_n)
    CHI.add_sim_opt("--chi_phi %E " % sim_param.chi_phi)
    make_job_scripts(CHI, 'CHI')

    # save_to_qsub(make_qsub(CHI), "scripts/CHI_qsub.pbs")


if __name__ == "__main__":
    sim_param = get_input()
    qsub_ZA(sim_param)
    qsub_FF(sim_param)
    qsub_FP(sim_param)
    qsub_FP_pp(sim_param)
    qsub_CHI(sim_param)

    # save_to_qsub(make_submit(), "scripts/submit_mlt.sh")
    # st = os.stat('scripts/submit_mlt.sh')
    # os.chmod('scripts/submit_mlt.sh', st.st_mode | stat.S_IEXEC)

    # save_to_qsub(make_stack_qsub(), "scripts/stack_qsub.pbs")
