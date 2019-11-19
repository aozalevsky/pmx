import errno
import logging
import luigi
import os
import subprocess
import pmx.scripts.workflows.SGE_tasks.SGETunedRunner as sge_runner
from luigi.contrib.sge import logger
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask #tuned for the owl cluster



def _parse_qsub_job_array_id(qsub_out):
    """Parse job id from qsub output string.

    Assume format:

        "Your job-array <job_id>.<min>-<max>:<step> ("<job_name>") has been submitted"

    """
    return int(qsub_out.split()[2].split('.')[0])

def _build_qsub_array_command(cmd, job_name, pe, n_cpu, work_dir, nsubtasks):
    """Submit array job shell command to SGE queue via `qsub`"""
    qsub_template = """echo {cmd} | qsub -t 1:{nsubtasks} -V -r y -pe {pe} {n_cpu} -N {job_name} -wd {work_dir}"""
    return qsub_template.format(
        cmd=cmd, job_name=job_name, work_dir=work_dir,
        pe=pe, n_cpu=n_cpu, nsubtasks=nsubtasks)



class SGETunedArrayJobTask(SGETunedJobTask):

    unfinished=[] #set this in __init__()

    def _run_job(self):

        # Build a qsub argument that will run sge_runner.py on the directory we've specified
        runner_path = sge_runner.__file__
        if runner_path.endswith("pyc"):
            runner_path = runner_path[:-3] + "py"
        job_str = 'python {0} "{1}" "{2}"'.format(
            runner_path, self.tmp_dir, os.getcwd())  # enclose tmp_dir in quotes to protect from special escape chars
        if self.no_tarball:
            job_str += ' --no-tarball'

            #force loading of dependencies by sourcing a custom profile
            if(os.path.isfile("~/.luigi_profile")):
                job_str = '"source ~/.luigi_profile; '+job_str+'"'
            else:
                logging.error("Tarballing of dependencies is disabled and "
                              "~/.luigi_profile does not exist. "
                              "Will not be able to load all workflow "
                              "dependencies without it. Please create it and "
                              "within activate a conda environment containing "
                              "at least python>3.6, "
                              "pmx, luigi, MDanalysis, matplotlib, and numpy.")
                raise FileNotFoundError(errno.ENOENT,
                              os.strerror(errno.ENOENT), "~/.luigi_profile")

        #tell runner that this is an array job
        job_str += ' --arrayjob'

        #force loading of conda and luigi by sourcing a custom profile
        job_str = '"source ~/.luigi_profile; '+job_str+'"'

        self.errfile=""; #no errorfile. mdrun dumps too much into stderr

        # Build qsub submit command
        submit_cmd = _build_qsub_array_command(job_str, self.job_name,
                                               self.parallel_env, self.n_cpu,
                                               self.sim_path,
                                               len(self.unfinished))
        logger.debug('qsub command: \n' + submit_cmd)

        # Submit the job and grab job ID
        output = subprocess.check_output(submit_cmd, shell=True).decode('utf-8')
        logger.debug("Submitted job to qsub with response:\n" + output)
        self.job_id = _parse_qsub_job_array_id(output)
        #logger.debug("Submitted job to qsub with response:\n" + output)

        self._track_job()

        # Now delete the temporaries, if they're there.
        if (self.tmp_dir and os.path.exists(self.tmp_dir) and not self.dont_remove_tmp_dir):
            logger.info('Removing temporary directory %s' % self.tmp_dir)
            subprocess.call(["rm", "-rf", self.tmp_dir])
