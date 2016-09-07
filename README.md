SQCT -- Single Qubit Circuit Toolkit
--------------------------------------

### BUILD

You will need the following libraries installed on your system: 

1. Boost 1.48
    - program_options 
    - chrono
    - timer
    - system
2. The GNU Multiple Precision Arithmetic Library (gmp and gmpxx)
3. The GNU MPFR Library (mpfr)

Also, a C++ compiler supporting C++11 is necessary.

#### Ubuntu Linux

The list of dependencies can be installed with:

    sudo apt-get install libgmp-dev libmpfr-dev libboost-program-options-dev \
    libpari-dev libboost-chrono-dev libboost-timer-dev libboost-system-dev \
    build-essential cmake git

SQCT can be downloaded from github:

    git clone https://github.com/vadym-kl/sqct.git

And needs to be compiled with `-DStaticLink=OFF`:

    cd sqct
    mkdir build
    cd build
    cmake -DStaticLink=OFF ..
    make

### USAGE

Information about program use available through --help option.

### ABOUT 

The program code based on results of http://arxiv.org/abs/1206.5236. It also implements 
the version of Solovay Kitaev algorithm described in http://arxiv.org/abs/quant-ph/0505030. 
In addition to Boost, The GNU Multiple Precision Arithmetic Library, The GNU MPFR Library the library 
mpfr::real by Christian Schneider <software(at)chschneider(dot)eu> is used for high precision
floating point arithmetic. 

### DIRECTORY STRUCTURE 

 - `sk` -- implementation of the Solovay-Kitaev algorithm
 - `es` -- exact synthesis algorithm
 - `theory` -- numerical proof of result from arXiv:1206.5236, tests of exact synthesis algorithm 
 - `appr` -- optimal round off of unitaries

