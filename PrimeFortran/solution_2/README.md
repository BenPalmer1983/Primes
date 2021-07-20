# Fortran solution by tjol

![Algorithm](https://img.shields.io/badge/Algorithm-base-green)
![Faithfulness](https://img.shields.io/badge/Faithful-yes-green)
![Parallelism](https://img.shields.io/badge/Parallel-no-green)
![Parallelism](https://img.shields.io/badge/Parallel-yes-green)
![Bit count](https://img.shields.io/badge/Bits-1-green)
![Bit count](https://img.shields.io/badge/Bits-8-yellowgreen)
![Bit count](https://img.shields.io/badge/Bits-unknown-yellowgreen)

This Fortran solution uses a Fortran 2003 class, unlike solution 1 by johandweber.
There are four versions:

 * `prime-bitarray`, the most faithful with 1 bit per flag and manual bit
   manipulation. 
 * `prime-8bit`, the fastest with an 8 bit integer per flag. 
 * `prime-logical-array`, which uses an array of `logical`.
 * `prime-bitarray-parallel`, a parallelized version of `prime-bitarray` by BenPalmer1983

## Run instructions

    make run

## Output

    tjol-bits;10416;5.00042915;1;algorithm=base,faithful=yes,bits=1
    tjol-8bit;14761;5.00036812;1;algorithm=base,faithful=yes,bits=8
    tjol-logical;8342;5.00035095;1;algorithm=base,faithful=yes
    benpalmer1983-bits-par;32748;5.00212574;12;algorithm=base,faithful=yes,bits=1
