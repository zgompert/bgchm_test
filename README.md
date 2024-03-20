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
This script runs [command_one.R](commands_one.R), which fits the model for each data set. The rdat files from each fit model (simulated data set) are then combined and summarized with [summarizeTest.R](summarizeTest.R).


# Simulations and analyses for genomic clines
Simulations with the genomic clines model as the generative model and a focus on cline variablity are in [simClines1.R](simClines1.R). Cline fitting is organized with the following wrapper script:

```bash
#!/bin/sh 
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=bgchm
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

max=15
## total = total number of jobs
total=50
count=0

cd /uufs/chpc.utah.edu/common/home/gompert-group4/projects/bgchhm_methods/genomicl_test

for ((j=1; j<=50; j++))
do
    Rscript --vanilla commands_one_var.R "$j" &
    ((count++))

    if ((count >= max)); then
        wait -n
        ((count--))
    fi
done

wait
```

This runs [commands_one_var.R](commands_one_var.R), and the results are then summarized with [summarizeTest1.R](summarizeTest1.R).

The next test focused on the effect of allele frequency differences between parental populations (ancestry informativness). We again used the genomic clines model as the generative model, see [simClines3.R](simClines3.R).

Cline fitting across data sets was organized with:

```bash
#!/bin/sh 
#SBATCH --time=96:00:00
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

cd /uufs/chpc.utah.edu/common/home/gompert-group4/projects/bgchhm_methods/genomicl_test

for ((j=1; j<=50; j++))
do
    Rscript --vanilla commands_one_afd_range.R "$j" &
    ((count++))

    if ((count >= max)); then
        wait -n
        ((count--))
    fi
done

wait
```
This runs [commands_one_afd_range.R](commands_one_afd_range.R), with results then summarized in [summarizeTest3.R](summarizeTest3.R).

The last test involved simulations with [dfuse](https://cbuerkle.bitbucket.io/software/dfuse/) and alternative genetic architectures for hybrid fitness.

The commands for running dfuse and simulating the hybrid zones were:

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


cd /uufs/chpc.utah.edu/common/home/gompert-group4/projects/bgchhm_methods/genomicl_test

## for testing
~/bin/dfuse_src/dfuse -d demefile_15 -s unsel_2 -r 10 -o df_neu -g 5000 -c 0 -G 1000 -m 0.05 -l 251 &
~/bin/dfuse_src/dfuse -d demefile_15 -s unsel_2 -r 10 -o df_oligo -g 5000 -c .3 -G 1000 -m 0.05 -l 251 & 
~/bin/dfuse_src/dfuse -d demefile_15 -s unsel_50 -r 10 -o df_poly -g 5000 -c .005 -G 1000 -m 0.05 -l 251 &
~/bin/dfuse_src/dfuse -d demefile_15 -s unsel_50 -r 10 -o df_strongpoly -g 5000 -c .01 -G 1000 -m 0.05 -l 251 &
```
Cline fitting for the resulting data sets was controlled by:

```bash
#!/bin/sh 
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=bgchm
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

max=8
## total = total number of jobs
total=10
count=0

cd /uufs/chpc.utah.edu/common/home/gompert-group4/projects/bgchhm_methods/genomicl_test

for ((j=1; j<=10; j++))
do
    Rscript --vanilla commands_neu.R "$j" &
    Rscript --vanilla commands_oligo.R "$j" &
    Rscript --vanilla commands_poly.R "$j" &
    Rscript --vanilla commands_spoly.R "$j" &
    ((count++))

    if ((count >= max)); then
        wait -n
        ((count--))
    fi
done

wait
```
This calls four R scripts (one per genetic architecture): [commands_neu.R](commands_neu.R), [commands_oligo.R](commands_oligo.R), [commands_poly.R](commands_poly.R), and [commands_spoly.R](commands_spoly.R).
