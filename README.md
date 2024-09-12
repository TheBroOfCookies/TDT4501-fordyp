## Run
```
mpirun -np 4 ./main_mpi 64 64 64 500 1
```
Runs a 64x64x64 point grid with 500 timesteps and sourcetype 1 and 4 processes

```
make; mpirun -np 4 ./main_mpi 64 64 64 500 1
```
Saved command for faster debugging


### Sourcetypes
1 = Stress monopole  
2 = Force monopole  
3 = Force dipole


### Tools
Comparing two files
```
diff file1 file2
```
