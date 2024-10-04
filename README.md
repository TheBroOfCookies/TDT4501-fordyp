## Run
```
mpirun -np 2 ./main_mpi 64 64 64 500 1
```
Runs with 64x64x64 points in a grid, 500 timesteps, sourcetype 1 and 2 processes

### Sourcetypes
1 = Stress monopole  
2 = Force monopole  
3 = Force dipole

### Notes
S - Stress
SXX stress double derivative with regards to X

MPI implementation divides problem space in three dimentions, first along Z-, then Y-, then the X-axis  
MPI implementation only works for saving recievers for problem size = 64x64x64 

## 3D neighbour exchange
Border exchange with diagonal neighbours in three dimentions gives us a total of 26 neighbours  
Divided by me into the following categories  
|Category| Number of neighbours|Classification|
|--------|:-------------------:|:------------:|
| Faces  |        6            | T1           |
| Edges  |        12           | T2           |
| Corners|        8            | T3           |

![Alt text](figures\neighbours\3d.png) Red = T1, Green = T2, Blue = T3

### Neighbor order  
Offset to neighbour from current rank in carthesian coordinates (using MPI_Cart_comm)  
Each neighbour is also given a number from 0 to 25  
Numbers are arranged in the following order T1->T2->T3   
#### Faces (T1)
| Nr| z | y | x |  | Nr| z | y | x |
| - |---|---|---|--|---|---|---|---|
| 0 | 1 | 0 | 0 |  | 1 |-1 | 0 | 0 |
| 2 | 0 | 1 | 0 |  | 3 |  0|-1 | 0 |
| 4 | 0 | 0 | 1 |  | 5 |  0| 0 |-1 |

#### Edges (T2)
| Nr| z | y | x |  | Nr| z | y | x |  | Nr| z | y | x |  | Nr| z | y | x |
|---|---|---|---|--|---|---|---|---|--|---|---|---|---|--|---|---|---|---|
| 6 | 1 | 1 | 0 |  | 7 | 1 |-1 | 0 |  | 8 | 1 | 0 | 1 |  | 9 | 1 | 0 |-1 |
|10 |-1 | 1 | 0 |  |11 |-1 |-1 | 0 |  |12 |-1 | 0 | 1 |  |13 |-1 | 0 |-1 |
|14 | 0 | 1 | 1 |  |15 | 0 |-1 | 1 |  |16 | 0 | 1 |-1 |  |17 | 0 |-1 |-1 |

#### Corners (T3)
| Nr| z | y | x |  | Nr| z | y | x |  | Nr| z | y | x |  | Nr| z | y | x |
|---|---|---|---|--|---|---|---|---|--|---|---|---|---|--|---|---|---|---|
|18 | 1 | 1 | 1 |  |19 | 1 | 1 |-1 |  |20 | 1 |-1 | 1 |  |21 | 1 |-1 |-1 |
|22 |-1 | 1 | 1 |  |23 |-1 | 1 |-1 |  |24 |-1 |-1 | 1 |  |25 |-1 |-1 |-1 |

![Alt text](figures\neighbours\3d_annotated_numbered.png)
Z, Y, X  (Nr)



## Commands list
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



### Debuggin clipboard
```
make; mpirun -np 1 ./main_mpi 64 64 64 500 1; mpirun -np 2 ./main_mpi 64 64 64 500 1; diff -y --suppress-common-lines receivers.csv receivers_mpi.csv | wc -l
```

```
diff receivers.csv receivers_fasit.csv

diff -y --suppress-common-lines receivers.csv receivers_fasit.csv | wc -l

diff -y --suppress-common-lines <(tail -n +3 receivers.csv) <(tail -n +3 receivers_fasit.csv) | wc -l


diff receivers.csv receivers_mpi.csv

diff -y --suppress-common-lines receivers.csv receivers_mpi.csv | wc -l
```