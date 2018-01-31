NAME
    mqsolver

DESCRIPTION
    C implementation of the Parallel Crossbred algorithm for solving
    Fukuoka MQ challenges on GPUs.
    
DOCUMENTATION
    For the thesis, see msc_thesis.pdf
    For the paper, see https://eprint.iacr.org/2017/1181

    Note that for the security strength estimation, one should refer to the paper
    as it provides more accurate results.

DEPENDENCIES
    cmake, make, python3, gcc, CUDA

OPTIONAL DEPENDENCIES
    ssh (for cluster mode)
    ruby (for some scripts in the bin dir)

OPTIONS
    see the help message of solve.py

EXAMPLE
    ./solve.py -d 3 -k 16 -t 20 -v -o 46-92-3-16.log challenge-46-92.txt

BUILD
    Instead of using the Python wrapper, you can manually build and launch
    mqsolver as follows:

        $ mkdir build && cd build
        $ cmake -DKEEP_VAR_NUM=16 .. && make
    
    Now, an executable named 'mqsolver' should be available in the build folder
