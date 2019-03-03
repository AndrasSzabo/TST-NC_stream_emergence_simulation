# TST-NC_stream_emergence_simulation
Cell-based computational model for simulating stream emergence of neural crest


The code is based on the Tissue Simulation Toolkit first described in the publication: 
R M H Merks, & J A Glazier. (2005). A cell-centered approach to developmental biology. Physica A, 352(1), 113–130. doi:10.1016/j.physa.2004.12.028, for which the code is available here: https://sourceforge.net/projects/tst/

To compile, use:
```
cd prog
qmake
make
```

To run simulations:
```
cd prog
tst ../runs/par
```

Parameters for the simulations are specified in a parameter file `par`. An example is provided in `runs/par`. 
