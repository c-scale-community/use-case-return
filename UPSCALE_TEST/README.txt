1) Python script (to be called from .sh script)

test_ts_yeoda_processing.py


2) Shell script which is a wrapper to send a job contianing the python script to the Slurm scheduler (which basically handles the distribution of work on the cluster)

test_ts_yeoda_processing.sh

Note: some info about Slurm on Spider see https://spiderdocs.readthedocs.io/en/latest/Pages/getting_started.html#submitting-a-job  (section 2.4)


3) To run these scripts
- check and/or update all the paths
- $ sbatch test_ts_yeoda_processing.sh
- if above works then a job with jobid is generated incl. a slurm-[jobid].out (contains the stdout of the script and is updated as the script runs)
- to check the queue for just your jobs: $ squeue -u $USER
- to check the full queue: $ squeue

* more Slurm commands see also the above Spider documentation. In addition Slurm has very extensive online documentation (but maybe too extensive) https://slurm.schedmd.com/


4) To scale up in terms of jobs it would a good idea to use Slurm job arrays. Some information can be found here (do note that not all features on these websites may work on Spider):

Some job array  examples:

https://rcpedia.stanford.edu/topicGuides/jobArrayPythonExample.html
https://crc.ku.edu/hpc/how-to/arrays
https://slurm.schedmd.com/job_array.html


5) image viewers on spider:
 - gwenview
 - display (image magick) *
 - ristretto
