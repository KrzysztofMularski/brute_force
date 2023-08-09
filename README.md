# Brute force calculating whole matrix
## Usage:
```
brute_force <TRAJECTORY> [OPTIONS]...
```
## Options:
Parameter|Description
-|-
`-s SIZE`              |setting matrix size to SIZE, default is max size for current trajectory
`-t THREADS`           |setting number of omp threads to THREADS, default is 1
`-i`                   |generating matrix image in gray scale:<br>- max rmsd: white pixel<br>- min rmsd: black pixel<br>- values are normalized<br>By default no image generated.

## Examples:
```
brute_force ./trajectory.pdb
```
```
brute_force ./trajectory.pdb -t 10 -i
```
## Results
Example result grey maps can be found under directory ![results](/results) in this project.

All maps were generated with 10 omp threads.
|md5sum|Trajectory file|Max RMSD|Frames pair|The calculation took|
|---|---|---|---|---|
|03...82|100 5ns|97.3844|(3, 312)|6m 43s|
|c7...15|100 15ns|106.677|(2912, 709)|1h 4m 41s|
|89...6f|303 5ns|222.681|(101, 876)|19m 6s|
|a2...b2|303 15ns|244.776|(0, 2516)|2h 25m 12s|

Other (old):
|md5sum|Trajectory file|Max RMSD|Frames pair|The calculation took|
|---|---|---|---|---|
|86...b0|100 5ns (no -ignh)|129.495|(0, 934)|6m 57s|

