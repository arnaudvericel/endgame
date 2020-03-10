# endgame
1D evolution code for growing and fragmenting dust in the presence of snow lines and particle trapping

This code computes the parallel evolution of n-dust grains in a static disc of gas.
The user has the possibility to add a few physical processes such as growth, fragmentation, snow lines, pressure bumps and size shift after sublimation.

The code requires 2 .in files : disc.in and dust.in.

disc.in states the disc properties and the physics involved, while dust.in has the number of grains you want to simulate as well as their initial radius, size and intrinsic densities.

The output is of the form of files such as "px.dat" where x is the dust grain number formated at 3 digits.
There is one output per grain throughout the simulation.

To use endgame, simply type make.

To use endgame, with fun, try out make THANOS=yes.

Enjoy!
