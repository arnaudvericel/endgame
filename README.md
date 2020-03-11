# endgame
1D evolution code for growing and fragmenting dust in the presence of snow lines and particle trapping.

This code computes the parallel evolution of n-dust grains in a static disc of gas, i.e. by integrating the dust radial drift velocity with the drift and viscous components.
The user has the possibility to add a few physical processes such as growth, fragmentation, snow lines, pressure bumps and size shift after sublimation/condensation.

The code requires 2 .in files : disc.in and dust.in.

disc.in states the disc properties and the physics involved, while dust.in specifies the number of grains you want to simulate as well as their initial radii, sizes and intrinsic densities.

The output is of the form of files such as "px.dat" where x is the dust grain number formated at 3 digits.
There is one output per grain throughout the simulation.

To make the executable "endgame", simply type make.

Thanos can also make some comments troughout the execution of the code. To allow Thanos to monologue sometimes, compile the code with THANOS=yes.

Enjoy!
