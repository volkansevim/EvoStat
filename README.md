EvoStat
=======

EvoStat is a reaction simulator that generates biochemical circuits with desired functions. 

Requires: SUNDIALS ODE Solver (libsundials-serial-dev package on Ubuntu), Boost libraries, GNU Scientific Library (GSL), Open-MPI

Run sundial-config to get required linker flags. 

Compile using `./compile.sh ode13`

Running the simulator
---------------------
 ```
 mpirun -n 8 ~/EvoStat/a.out reactants.dat reactions.dat \
  fitnessparam:10 \
  fitnessparam2:0 \
  runtime:75 \
  popsize:20 \
  generations:10 \
  seed:12345678 \
  randomize:1 \
  fitnessavgover:1 \
  extrinsicnoise:0 \
  sigma:0.0 \
  lognoisemax:1.2 \
  mureacrate:0.3 \
  maxrate:0.0017001 \
  MaxRndRate:0.001 \
  MinRndRate:0.00001 > testevolve.txt
```

