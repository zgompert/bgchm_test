# bgchm_test

# Simulations and analyses for hybrid index and ancestry class proportions
Simulate data with dfuse.
```bash
#!/bin/sh 
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=dfuse
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu


cd /uufs/chpc.utah.edu/common/home/gompert-group4/projects/bgchhm_methods/hq_test

## for illustration of results
~/bin/dfuse_src/dfuse -d demefile_1 -s nullsel -o df -g 200 -c 0 -G 100 -m 0.1 -l 101

## for testing
~/bin/dfuse_src/dfuse -d demefile_1 -s nullsel -r 50 -o repdf -g 200 -c 0 -G 200 -m 0.1 -l 51
```

Wrapper to fit models.
```bash
#!/bin/sh 
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=bgchm
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

max=20
## total = total number of jobs
total=50
count=0

cd /uufs/chpc.utah.edu/common/home/gompert-group4/projects/bgchhm_methods/hq_test

for ((j=1; j<=50; j++))
do
    Rscript --vanilla commands_one.R "$j" &
    ((count++))

    if ((count >= max)); then
        wait -n
        ((count--))
    fi
done

wait
```
