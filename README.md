# SimByDCA

### Installation 

The code present in Code_for_cluster was present in the EPFL cluster. It is helpful to generate data and to infer contact and partners. The model has been inferred with bmDCA on the cluster and arDCA on my personal computer.
If you want to reproduce the data, I advise you to copy the folder Code_for_cluster and to do :

```
cd 21_states/Code_for_cluster/cython_code/
cythonize -i generation_sequence.pyx
cythonize -i generation_sequence_arDCA.pyx
cythonize -i analyse_sequence.pyx 
```

