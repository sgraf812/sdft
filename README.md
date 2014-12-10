SDFT
====

This tiny C library computes an N-length DFT with the help of a sliding window DFT (SDFT).

How To Build
============

This project uses CMake for project file generation. I'm using the EAP of CLion for working on this project, which supports CMake out of the box, but here are the steps to perform:

  1. `$ cd /path/to/sdft`
  2. `$ mkdir build` to have a folder where the project files can reside in an out of tree build
  3. `$ cd build`
  4. `$ cmake ..`

Now, there should be appropriate project files in the `build/` directory, depending on which target cmake chose for you (or you chose), with which compilation should be straightforward (e.g. `$ make` or opening it in Visual Studio).

How To Use
==========
For instructions on how to use it, dig into test/main.c:compare_sdft_to_dft and read through the docstrings.
