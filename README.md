# Brute force calculating whole matrix
## Usage:
```
brute_force <TRAJECTORY> [OPTIONS]...
```
## Options:
Parameter|Type:Default|Description
-|-|-
`-s SIZE`              |setting matrix size to SIZE, default is max size for current trajectory
`-t THREADS`           |setting number of omp threads to THREADS, default is 1
`-i`                   |generating matrix image in gray scale:
                       | - max rmsd: white pixel,
                       | - min rmsd: black pixel,
                       | - values are normalized.
                       | By default no image generated.

## Examples:
```
brute_force ./trajectory.pdb
```
```
brute_force ./trajectory.pdb -t 10 -i
```
