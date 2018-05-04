This folder includes the software needed to run the sampling process.

The worker program has to be compiled first by running `make` from within the folder `hexsampler`. The compiled executable (also called `hexsampler`) then has to be copied to the target computer (which can be a remote computer) into a folder chosen by the user (usually the `home` folder).

The file `nodesSampling.txt` has to be copied into the root folder of ImageJ. It is used to tell the plug-in `Hex_Sampler` (which is the master process for the reconstruction) where the workers (`hexsampler`) can be found. The example `nodesSampling.txt` file included here just specifies locations on the local computer, but if more computers are available on your network make sure to include them as well because it will significantly speed up the sampling process. Usually you want the number of workers per computer (first number in each line, column multiplicity) to match the number of cores of this computer.

Note that compiling the software with clang instead of g++ might give you a speedup on your computer. This can be changed at the top of the Makefile.
