### Connect to idun
```
ssh trygvaan@idun-login1.hpc.ntnu.no
```

### Idun my commands
```
sbatch optewe.slurm; watch squeue -u trygvaan
```


### Idun general commands
```
squeue -u <username>   # List opp jobbene til <username>  
squeue                 # List opp alle jobber på Idun  
sbatch jobbfil.slurm   # Send inn 'jobbfil.slurm' til køen, eksempel i 'test_job.slurm' fra tar-filen nedenfor.  
scancel -j 20094329    # Avbryt og kast jobb 20094329 (såfremt du eier den). Nummeret til en jobb finner man med squeue  
sinfo -Nl              # Se listen over alle nodene og hvilken tilstand de er i 
#SBATCH --constraint=v100 #for å få noder med feature "v100" på
```