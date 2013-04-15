#!/bin/bash

# This function writes a slurm script.
# We can call it with different parameter
# settings to create different experiments

function write_script
{
TASKS=$[${NPERNODE} * ${NODES}]
JOB_NAME=$(printf 'hdf5-%05d-%02d-%05d-%05d-%03d' ${NODES} ${NPERNODE} ${CB_BUFSIZE} ${STRIPESIZE} ${MEMPERTASK})
DIR_NAME=${BASEDIR}/${JOB_NAME}

if [ -f ${DIR_NAME}/slurm.out ] ; then
        echo "${DIR_NAME}/slurm.out already exists, skipping..."
        return 0
else
        echo "Creating job ${JOB_NAME}"
fi

mkdir -p ${DIR_NAME}

cat << _EOF_ > ${DIR_NAME}/job_script.bash
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=${DIR_NAME}/slurm.out
#SBATCH --error=${DIR_NAME}/slurm.err
#SBATCH --time=04:00:00
#SBATCH --nodes=${NODES}
#SBATCH --ntasks-per-core=2
#SBATCH --distribution=cyclic:cyclic
##SBATCH --partition=${QUEUE}
##SBATCH --reservation=craytestmon

#MPIIO hints
#    cb_buffer_size      = 16777216
#    romio_cb_read       = automatic
#    romio_cb_write      = automatic
#    cb_nodes            = #nodes/8
#    romio_no_indep_rw   = false
#    ind_rd_buffer_size  = 4194304
#    ind_wr_buffer_size  = 524288
#    romio_ds_read       = automatic
#    romio_ds_write      = automatic
#    direct_io           = false
#    cb_config_list      = *:1

export MPICH_MPIIO_HINTS_DISPLAY=1
export MPICH_MPIIO_HINTS="*:romio_cb_write=enable:romio_ds_write=disable:cb_buffer_size=${CB_BUFSIZE}M:romio_no_indep_rw=true"
export MPICH_MPIIO_CB_ALIGN=2

# set stripe to -1 to use All OSS nodes
rm -rf ${DIR_NAME}/temp.h5
lfs setstripe -s ${STRIPESIZE}M -i 0 -c -1 ${DIR_NAME}

cd /scratch/daint/biddisco/build/pv-plugins

#
# -j2 tells slurm to use 2 threads per CPU as they are hyperthreaded
#

/usr/bin/aprun -j 2 -n ${TASKS} ${EXECUTABLE} -numNodes ${NODES} -processesPerNode ${NPERNODE} -memoryMB ${MEMPERTASK} -iterations 4 -pieceValidation 1 -F temp.h5 -D ${DIR_NAME} -T ${DIR_NAME}
-testName ${JOB_NAME}

rm -rf ${DIR_NAME}/temp.h5
_EOF_

chmod 775 ${DIR_NAME}/job_script.bash

echo "sbatch ${DIR_NAME}/job_script.bash" >> run_jobs.bash
}

# use current directory as base dir
#pushd `dirname $0` > /dev/null
pushd `pwd` > /dev/null
BASEDIR=`pwd`
popd > /dev/null
echo "Generating jobs using base directory ${BASEDIR}"

# Create another script to submit all generated jobs to the scheduler
echo "#!/bin/bash" > run_jobs.bash
echo "BASEDIR=${BASEDIR}" >> run_jobs.bash
chmod 775 run_jobs.bash

# set some vars which are fixed in this test
QUEUE=day
EXECUTABLE=/scratch/daint/biddisco/build/pv-plugins/bin/TestH5PartParallelWriter
CB_BUFSIZE=
STRIPESIZE=
MEMPERTASK=

# Loop through all the parameter combinations generating jobs for each
for MEMPERTASK in 32 64 128
do
    for CB_BUFSIZE in 1 4 16 64 128 
    do
        for STRIPESIZE in 1 4 16 64 128 
        do  
            NPERNODE=32
            for NODES in 64 128 256 512 1024 2048
            do
                write_script
            done
            
        done
    done
done
