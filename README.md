## Run
```
mpirun -np 4 ./main_mpi 64 64 64 500 1
```
Runs a 64x64x64 point grid with 500 timesteps and sourcetype 1 and 4 processes





### Sourcetypes
1 = Stress monopole  
2 = Force monopole  
3 = Force dipole





### Notes
S - Stress
SXX stress double derivative mhp X

Current MPI implementation divides problem space on Z-coordinates
Current MPI implementation only works for saving recievers for problem size = 64x64x64
Nz = 1024 no longer supported for save recievers

### Commands list
Comparing two files
```
diff file1 file2
```

Count number of lines that differ
```
diff -y --suppress-common-lines file1 file2 | wc -l
```
Ignore the first 2 lines
```
diff <(tail -n +3 file1) <(tail -n +3 file2)
```

Saved command for checking number of lines
```
diff receivers.csv receivers_fasit.csv

diff -y --suppress-common-lines receivers.csv receivers_fasit.csv | wc -l


diff receivers.csv receivers_mpi.csv

diff -y --suppress-common-lines receivers.csv receivers_mpi.csv | wc -l

diff -y --suppress-common-lines <(tail -n +3 receivers.csv) <(tail -n +3 receivers_fasit.csv) | wc -l
```


debuggin shorthand
```
make; mpirun -np 1 ./main_mpi 64 64 64 500 1; mpirun -np 2 ./main_mpi 64 64 64 500 1; diff -y --suppress-common-lines receivers.csv receivers_mpi.csv | wc -l
```