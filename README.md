# The MiniWeather Mini App

A lightweight Fortran application for simulating weather-like flows, developed
by ORNL and NVIDIA. It supports multiple parallelization strategies:

* MPI
* OpenACC Offload 
* OpenMP Threading 

Original Author: Matt Norman, Oak Ridge National Laboratory,
https://mrnorman.github.io

Author for DKRZ Levante: Jared Frazier, Leibniz Institute of Atmospheric
Physics, https://jfdev001.github.io/

# Table of Contents

- [Introduction](#introduction)
  * [Brief Description of the Code](#brief-description-of-the-code)
- [Compiling and Running the Code](#compiling-and-running-the-code)
  * [Software Dependencies](#software-dependencies)
  * [Basic Setup](#basic-setup)
  * [Building and Testing Workflow](#building-and-testing-workflow)
  * [Altering the Code's Configurations](#altering-the-codes-configurations)
  * [Running the Code](#running-the-code)
  * [Running Performance Experiments](#running-performance-experiments)
  * [Viewing the Output](#viewing-the-output)
- [Parallelization](#parallelization)
  * [Indexing](#indexing)
  * **[MPI Domain Decomposition](#mpi-domain-decomposition)**
  * **[OpenMP CPU Threading](#openmp-cpu-threading)**
  * **[OpenACC Accelerator Threading](#openacc-accelerator-threading)**
  * **[OpenMP Offload Accelerator Threading](#openmp-offload-accelerator-threading)**
- [Numerical Experiments](#numerical-experiments)
  * [Rising Thermal](#rising-thermal)
  * [Colliding Thermals](#colliding-thermals)
  * [Mountain Gravity Waves](#mountain-gravity-waves)
  * [Density Current](#density-current)
  * [Injection](#injection)
- [Physics, PDEs, and Numerical Approximations](#physics--pdes--and-numerical-approximations)
  * [The 2-D Euler Equations](#the-2-d-euler-equations)
  * [Maintaining Hydrostatic Balance](#maintaining-hydrostatic-balance)
  * [Dimensional Splitting](#dimensional-splitting)
  * [Finite-Volume Spatial Discretization](#finite-volume-spatial-discretization)
  * [Runge-Kutta Time Integration](#runge-kutta-time-integration)
  * [Hyper-viscosity](#hyper-viscosity)
- [MiniWeather Model Scaling Details](#miniweather-model-scaling-details)
- [Checking for Correctness](#checking-for-correctness)


# Introduction

The miniWeather code mimics the basic dynamics seen in atmospheric weather and
climate. The dynamics themselves are dry compressible, stratified,
non-hydrostatic flows dominated by buoyant forces that are relatively small
perturbations on a hydrostatic background state. The equations in this code
themselves form the backbone of pretty much all fluid dynamics codes, and this
particular flavor forms the base of all weather and climate modeling.

With about 500 total lines of code (and only about 200 lines that you care
about), it serves as an approachable place to learn parallelization and porting
using MPI + X, where X is OpenMP, OpenACC, CUDA, or potentially other
approaches to CPU and accelerated parallelization. The code uses periodic
boundary conditions in the x-direction and solid wall boundary conditions in
the z-direction. 

## Brief Description of the Code

### Domain Parameters

A fuller description of the science, math, and dynamics are play are given
later, but this section is reserved to describing some of the main variables
and flow in the code. The code is decomposed in two spatial dimensions, x and
z, with `nx_glob` and `nz_glob` cells in the global domain and nx and nz cells
in the local domain, using straightforward domain decomposition for MPI-level
parallelization. The global domain is of size xlen and zlen meters, and hs
“halo” cells are appended to both sides of each dimension for convenience in
forming stencils of cells for reconstruction.

### Fluid State Variables

There are four main arrays used in this code: `state`, `state_tmp`, `flux`, and `tend`, and the dimensions for each are given in the code upon declaration in the comments. Each of these arrays is described briefly below:

* `state`: This is the fluid state at the current time step, and it is the only array that persists from one time step to the next. The other four are only used within the calculations to advance the model to the next time step. The fluid state describes the average state over each cell area in the spatial domain. This variable contains four fluid states, which are the traditional mass, momenta, and thermodynamic quantities of most fluid models:
  1. Density (`ID_DENS`): The 2-D density of the fluid, <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\rho" title="\large \rho" />, in <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\text{kg}\&space;\text{m}^{-2}" title="\large \text{kg}\ \text{m}^{-2}" /> (note this is normally <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\text{kg}\&space;\text{m}^{-3}" title="\large \text{kg}\ \text{m}^{-3}" />, but this is a 2-D model, not 3-D)
  2. U-momentum (`ID_UMOM`): The momentum per unit area of the fluid in the x-direction calculated as <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\rho&space;u" title="\large \rho u" />, where u is the x-direction wind velocity. The units are <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\text{kg}\&space;\text{m}^{-1}\&space;\text{s}^{-1}" title="\large \text{kg}\ \text{m}^{-1}\ \text{s}^{-1}" />. Note that to get true momentum, you must integrate over the cell.
  2. W-momentum (`ID_WMOM`): The momentum per unit area of the fluid in the z-direction calculated as <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\rho&space;w" title="\large \rho w" />, where w is the z-direction wind velocity. The units are <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\text{kg}\&space;\text{m}^{-1}\&space;\text{s}^{-1}" title="\large \text{kg}\ \text{m}^{-1}\ \text{s}^{-1}" />. Note that to get true momentum, you must integrate over the cell.
  4. Potential Temperature (`ID_RHOT`): The product of density and potential temperature, <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\rho&space;\theta" title="\large \rho \theta" />, where <img src="https://latex.codecogs.com/svg.latex?\theta=T\left(P_{0}/P\right)^{R_{d}/c_{p}}" /> , <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;P_{0}=10^{5}\,\text{Pa}" title="\large P_{0}=10^{5}\,\text{Pa}" />, T is the true temperature, and <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;R_d" title="\large R_d" /> and<img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;c_p" title="\large c_p" /> are the dry air constant and specific heat at constant pressure for dry air, respectively. The units of this quantity are <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\text{K}\,\text{kg}\,\text{m}^{-2}" title="\large \text{K}\,\text{kg}\,\text{m}^{-2}" />.
* `state_tmp`: This is a temporary copy of the fluid state used in the Runge-Kutta integration to keep from overwriting the state at the beginning of the time step, and it has the same units and meaning.
* `flux`: This is fluid state at cell boundaries in the x- and z-directions, and the units and meanings are the same as for `state` and `state_tmp`. In the x-direction update, the values of `flux` at indices `i` and `i+1` represents the fluid state at the left- and right-hand boundaries of cell `i`. The indexing is analogous in the z-direction. The fluxes are used to exchange fluid properties with neighboring cells.
* `tend`: This is the time tendency of the fluid state <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\partial\mathbf{q}/\partial&space;t" title="\large \partial\mathbf{q}/\partial t" />, where <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\mathbf{q}" title="\large \mathbf{q}" /> is the the state vector, and as the name suggests, it has the same meaning and units as state, except per unit time (appending <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\text{s}^{-1}" title="\large \text{s}^{-1}" /> to the units). In the Finite-Volume method, the time tendency of a cell is equivalent to the divergence of the flux across a cell.

### Overall Model Flow

This code follows a traditional Finite-Volume control flow.

To compute the time tendency, given an initial state at the beginning of a time
step that contains cell-averages of the fluid state, the value of the state at
cell boundaries is reconstructed using a stencil. Then, a viscous term is added
to the fluid state at the cell boundaries to improve stability. Next, the
tendencies are computed as the divergence of the flux across the cell. Finally,
the tendencies are applied to the fluid state.

Once the time tendency is computed, the fluid PDEs are essentially now cast as
a set of ODEs varying only in time, and this is called the “semi-discretized”
form of the equations. To solve the equations in time, an ODE solver is used,
in this case, a Runge-Kutta method. Finally, at the highest level, the
equations are split into separate solves in each dimension, x and z.

# Compiling and Running the Code

## Software Dependencies

* Parallel-netcdf: https://github.com/Parallel-NetCDF/PnetCDF
  * This is a dependency for two reasons: (1) NetCDF files are easy to visualize and convenient to work with; (2) The users of this code shouldn't have to write their own parallel I/O.
* Ncview: http://meteora.ucsd.edu/~pierce/ncview_home_page.html
  * This is the easiest way to visualize NetCDF files.
* MPI
* For OpenACC: An OpenACC-capable compiler (PGI / Nvidia, Cray, GNU)
  * A free version of the PGI / Nvidia compiler can be obtained by googling for the "Community Edition"
* For OpenMP: An OpenMP offload capable compiler (Cray, XL, GNU)
* CMake: https://cmake.org

## Basic Setup

```shell
cd /work/bm1233/${USER}  
git clone git@github.com:jfdev001/miniWeather.git 
cd miniWeather
git submodule update --init --recursive
```

To find that repository on GitHub, go to  

```text
https://github.com/jfdev001/miniWeather
```

and star it so that you can easily find it later.

If you prefer, you can fork (see [github docs: fork a
repo](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo))
that repo and then clone your own fork of `miniweather`. This is also a good
approach since then you can upload (i.e., `git push`) your code to a repository
on your GitHub. The workflow would look something like the following

```shell
cd /work/bm1233/${USER}  
# assuming you have forked miniweather
git clone git@github.com:YOUR_GITHUB_USER_NAME_HERE/miniWeather.git 
cd miniWeather
git submodule update --init --recursive

# by default, the remote origin (i.e., source of your code on GitHub
# and its destination when pushing) is set to 
# git@github.com:YOUR_GITHUB_USER_NAME_HERE/miniWeather.git ...
# by adding another remote and calling it upstream, you can at any
# time inspect or pull the code provided by jfdev001 without breaking
# you local changes 
git remote add upstream git@github.com:jfdev001/miniWeather.git 
git fetch upstream  # allows you to pull code from jfdev001 in the future
```

## Building and Testing Workflow

The first thing you should do is verify that you can compile and run
`miniweather`:

```shell
cd miniWeather/fortran/build
source cmake_levante_test # or bash cmake_levante_test or ./cmake_levante_test
```

This generates a directory called `build_output/test` where all configuration
(e.g., auto-generated Makefiles) and compilation artifacts (e.g., executable
binaries like `serial`, `openmp`, and `mpi`).

You should *always* read the usage documentation for any script you run. For
nearly every script provided, you can do the following to get usage
documentation

```shell
# assuming you're in the directory of the script...
# also if you encounter permissions issues, call `chmod +x <script_name>`
./<script_name> -h
```

Pay close attention also to any examples section in usage documentation.

If you are using MobaXTerm you should have two ssh sessions open and connected
to Levante: (1) one just for looking at the uage documentation of scripts so
that you know what the inputs and outputs are, and (2) one for actually
launching scripts and/or editing files. If you are Linux or Mac, you can also
launch multiple terminals and each of them can independently ssh to Levante.
Levante also comes with `tmux` (terminal multiplexer) by default, so you could
use that instead if you're already familiar.

You can check to see what `cmake_levante_test` by typing

```shell
./cmake_levante_test -h
```

Note that the `cmake_levante_test` simply wraps the
`cmake_levante_config_and_build` script discussed in the following sections.

In the real world, there may be limited or *no* documentation for software that
you are using. Even worse, they're may be documentation but it could be *out of
date*. This is almost worse than having no documentation because you might
think the software is doing one thing while it is in reality doing something
completely unexpected. You need to be prepared to *read* through code to
determine what it does. When doing this, keep mind the following:

* What are the inputs to the code? 
    * Does it take positional arguments? Does it take flags? 
    * Does it have optional arguments? What are the defaults of those optional
    arguments? 
    * Does it read from a file?
    * What assumptions (if any) about paths (e.g., directories) does the code
    make?

* What are the outputs of the code?
    * Does it write to standard output (i.e., the terminal)?
    * Does it produce large files (e.g., netcdf)? Where do those files get 
    written?

* Does the script launch Slurm jobs or execute locally?


## Altering the Code's Configurations

To alter the configuration of the code, you can control the number of cells in
the x- and z-directions, the length of simulation time, the output frequency,
and the initial data to use by calling `cmake_levante_build_and_configure`,
see 

```shell
# assuming in build/ dir
./cmake_levante_build_and_configure -h
```

This script forwards arguments to two calls to `cmake` that configure and build
the script. The generated `cmake` configuration call might look like the below:

* `-DNX=400`: Uses 400 cells in the x-direction
* `-DNZ=200`: Uses 200 cells in the z-direction
* `-DSIM_TIME=1000`: Simulates for 1,000 seconds model time
* `-DOUT_FREQ=10`: Outputs every 10 seconds model time
* `-DDATA_SPEC=DATA_SPEC_THERMAL`: Initializes a rising thermal

It's best if you keep `NX` exactly twice the value of `NZ` since the domain is
20km x 10km. 

The data specifications are `DATA_SPEC_COLLISION`, `DATA_SPEC_THERMAL`,
`DATA_SPEC_MOUNTAIN`, `DATA_SPEC_DENSITY_CURRENT`, and `DATA_SPEC_INJECTION`,
and each are described later on.

## Running the Code

It is recommended to run the `build_output/test/serial_test` code first to get
an idea for the model outputs. Note that a file `output.nc` is always produced 
in the directory from which you call the `miniweather` executables.

As an example:

```shell
# assuming in fortran/ directory... this produces `output.nc` there
./build/build_output/test/serial_test
```

Since parameters are set in the code itself, you don't need to pass any
parameters to the executables.

This is fine for testing lightweight serial codes; however, we are interested
in parallel codes. Since you are using a compute cluster shared by many people,
jobs requiring more computational resources must be submitted to the Slurm
scheduler.

We provide a script that wraps the generation of a run script which you can
use later for running simulations. You can check out the parameters for this
script here:

```shell
./scripts/templates/make_run_scripts -h
```

This script can be used to generate Slurm scripts specific to your user for
running `miniweather` simulations. These scripts are, by convention, written to
`scripts/run` and are *not* tracked by `git`. If you wish to modify the
`.gitignore` file and remove the line containing `*.run`, `git` will not track
your generated run scripts. The run scripts will also be prefixed with the
partition that you have requested. Different partitions on levante (e.g.,
shared, compute, gpu) give the user differ compute resources. By default the
compute partition is used, and this will work with MPI and OpenMP jobs.

You should generate an example run script with the following:

```shell
EMAIL_HERE="put_your_email@gmail.com"
./scripts/templates/make_run_scripts ${EMAIL_HERE}
```

This generates `scripts/run/compute_miniweather.run`. You should inspect what
this script does with

```shell
scripts/run/compute_miniweather.run -h
```

In particular, you should run each of the `bash` examples in the usage doc to
get an understanding of what *would* be submitted to Slurm. There are lots of
outputs so make sure to read and understand them. When ready to submit jobs to
Slurm, look at the `sbatch` examples in the same usage doc.

Note that the `time` sbatch directive is set to 30 seconds. This is sufficient
for running tests, but may not be sufficient for running larger scale
simulations. If your simulation significantly exceeds the amount of time
allocated by the `time` directive, the simulation will timeout and you the
outputs to `output.nc` may be incomplete. Be aware that increasing the amount
of time you would like to run your job may result in you waiting longer for
the Slurm scheduler to actually launch your job. You should always prototype
any experiments or scripts that you write which involve Slurm such that they
request a very short amount of time (i.e., less than 1 minute).

## Running Performance Experiments

You may want to evaluate how the performance of `miniweather` is affected by
increasing the number of threads, increasing the number of MPI processes, or
doing a combination of both. You can inspect a sample bash script that prepares
and launches such experiments:

```shell
./scripts/scaling/launch_sample_scaling_experiments -h
```

You can use that script as a template for running your own experiments.

## Visualizing Performance Results

This will also depend heavily on the types of experiments that you wish to run,
however, an example python code that can be launched by:

```shell
python scripts/viz/sample_scaling_results.py
```

That script has no `-h` option supported; however, at the top of the file
is a small description of the contents of the script itself and what it's for.

You copy/modify it to accomplish your plotting goals for your experiments.

Below is an example output from the script:

<img width="999" height="799" alt="miniweather_openmp" src="https://github.com/user-attachments/assets/5f2959bf-393a-4ae2-8008-67383dffcc01" />


## Viewing the Output

The file I/O is done in the netCDF format: (https://www.unidata.ucar.edu/software/netcdf). To me, the easiest way to view the data is to use a tool called “ncview” (http://meteora.ucsd.edu/~pierce/ncview_home_page.html). To use it, you can simply type `ncview output.nc`, making sure you have X-forwarding enabled in your ssh session. Further, you can call `ncview -frames output.nc`, and it will dump out all of your frames in the native resolution you're viewing the data in, and you you can render a movie with tools like `ffmpeg`. 

# Parallelization

This code was designed to parallelize with MPI first and then OpenMP, OpenACC, OpenMP offload, or `parallel_for` next, but you can always parallelize with OpenMP or OpenACC without MPI if you want. But it is rewarding to be able to run it on multiple nodes at higher resolution for more and sharper eddies in the dynamics.

As you port the code, you'll want to change relatively little code at a time, re-compile, re-run, and look at the output to see that you're still getting the right answer. There are advantages to using a visual tool to check the answer (e.g., `ncview`), as it can sometimes give you clues as to why you're not getting the right answer. 

Note that you only need to make changes code within the first 450 source lines for C and Fortran, and each loop that needs threading is decorated with a `// THREAD ME` comment. Everything below that is initialization and I/O code that doesn't need to be parallelized (unless you want to) for C and Fortran directives-based approaches.

## Indexing

The code makes room for so-called “halo” cells in the fluid state. This is a common practice in any algorithm that uses stencil-based reconstruction to estimate variation within a domain. In this code, there are `hs` halo cells on either side of each spatial dimension, and I pretty much hard-code `hs=2`.

### Fortran

In the Fortran code's fluid state (`state`), the x- and z-dimensions are dimensioned as multi-dimensional arrays that range from `1-hs:nx+hs`. In the x-direction, `1-hs:0` belong to the MPI task to the left, cells `1:nx` belong to the current MPI task, and `nx+1:nx+hs` belong to the MPI task to the right. In the z-dimension, `1-hs:0` are artificially set to mimic a solid wall boundary condition at the bottom, and `nz+1:nz+hs` are the same for the top boundary. The cell-interface fluxes (`flux`) are dimensioned as `1:nx+1` and `1:nz+1` in the x- and z-directions, and the cell average tendencies (`tend`) are dimensioned `1:nx` and `1:nz` in the x- and z-directions. The cell of index `i` will have left- and right-hand interface fluxes of index `i` and `i+1`, respectively, and it will be evolved by the tendency at index `i`. The analog of this is also true in the z-direction.

## MPI Domain Decomposition

This code was designed to use domain decomposition, where each MPI rank “owns” a certain set of cells in the x-direction and contains two “halo” cells from the left- and right-hand MPI tasks in the x-direction as well. The domain is only decomposed in the x-direction and not the z-direction for simplicity.

**IMPORTANT**: Please be sure to set `nranks`, `myrank`, `nx`, `i_beg`, `left_rank`, and `right_rank`. These are clearly marked in the serial source code. You can set more variables, but these are used elsewhere in the code (particularly in the parallel file I/O), so they must be set.

To parallelize with MPI, there are only two places in the code that need to be altered. The first is the initialization, a subroutine / function named `init`, where you must determine the number of ranks, you process's rank, the beginning index of your rank's first cell in the x-direction, the number of x-direction cells your rank will own, and the MPI rank IDs that are to your left and your right. Because the code is periodic in the x-direction, your left and right neighboring ranks will wrap around. For instance, if your are rank `0`, your left-most rank will be `nranks-1`.

The second place is in the routine that sets the halo values in the x-direction. In this routine, you need to:

1. Create MPI data buffers (at the same place the other arrays are declared) to hold the data that needs to be sent and received, allocate them in the `init()` routine, and deallocate them in the `finalize()` routine.

2. Pack the data you need to send to your left and right MPI neighbors

3. Send the data to your left and right MPI neighbors

4. Receive the data from your left and right MPI neighbors

5. Unpack the data from your left and right neighbors and place the data into your MPI rank's halo cells. 

Once you complete this, the code will be fully parallelized in MPI. Both of the places you need to add code for MPI are marked in the serial code, and there are some extra hints in the `set_halo_values_x()` routine as well.

## OpenMP CPU Threading

For the OpenMP code, you basically need to decorate the loops with `omp parallel do` in Fortran or `omp parallel for` in C, and pay attention to any variables you need to make `private()` so that each thread has its own copy. Keep in mind that OpenMP works best on “outer” loops rather than “inner” loops. Also, for sake of performance, there are a couple of instances where it is wise to use the “collapse” clause because the outermost loop is not large enough to support the number of threads most CPUs have.

In Fortran, you can parallelize three loops with the following directive:

```fortran
!$omp parallel do collapse(3)
do ll = 1 , NUM_VARS
  do k = 1 , nz
    do i = 1 , nx
      state_out(i,k,ll) = state_init(i,k,ll) + dt * tend(i,k,ll)
    enddo
  enddo
enddo
```

This will collapse the three loops together (combining their parallelism) and then launch that parallelism among a number of CPU threads.

## OpenACC Accelerator Threading

To thread the same loops among the threads on a GPU, you will use the following in Fortran:

```fortran
!$acc parallel loop collapse(3)
do ll = 1 , NUM_VARS
  do k = 1 , nz
    do i = 1 , nx
      state_out(i,k,ll) = state_init(i,k,ll) + dt * tend(i,k,ll)
    enddo
  enddo
enddo
```

The OpenACC approach will differ depending on whether you're in Fortran or C. Just a forewarning, OpenACC is much more convenient in Fortran when it comes to data movement because in Fortran, the compiler knows how big your arrays are, and therefore the compiler can (and does) create all of the data movement for you (NOTE: This is true for PGI and Cray but not for GNU at the moment). All you have to do is optimize the data movement after the fact. for more information about the OpenACC copy directives, see:

https://github.com/mrnorman/miniWeather/wiki/A-Practical-Introduction-to-GPU-Refactoring-in-Fortran-with-Directives-for-Climate#optimizing--managing-data-movement

### Fortran Code

The OpenACC parallelization is a bit more involved when it comes to performance. But, it's a good idea to just start with the kernels themselves, since the compiler will generate all of your data statements for you on a per-kernel basis. You need to pay attention to private variables here as well. Only arrays need to be privatized. Scalars are automatically privatized for you.

Once you're getting the right answer with the kernels on the GPU, you can look at optimizing data movement by putting in data statements. I recommend putting data statements for the `state`, `tend`, `flux`, `hy_*`, and the MPI buffers (`sendbuf_l`, `sendbuf_r`, `recvbuf_l`, and `recvbuf_r`) around the main time stepping loop. Then, you need to move the data to the host before sending MPI data, back to the device once you receive MPI data, and to the host before file I/O.

## OpenMP Offload Accelerator Threading

To launch the same loops in OpenMP offload on a GPU's threads, you will use:

```fortran
!$omp target teams distribute parallel do simd collapse(3)
do ll = 1 , NUM_VARS
  do k = 1 , nz
    do i = 1 , nx
      state_out(i,k,ll) = state_init(i,k,ll) + dt * tend(i,k,ll)
    enddo
  enddo
enddo
```

Note that some compilers do different things for `simd` and `parallel for`, and therefore, there is no portable use of OpenMP offload at this point. In C / C++, this will be:

```C++
#pragma omp target teams distribute parallel for simd collapse(3)
for (ll=0; ll<NUM_VARS; ll++) {
  for (k=0; k<nz; k++) {
    for (i=0; i<nx; i++) {
      inds = ll*(nz+2*hs)*(nx+2*hs) + (k+hs)*(nx+2*hs) + i+hs;
      indt = ll*nz*nx + k*nx + i;
      state_out[inds] = state_init[inds] + dt * tend[indt];
    }
  }
}
```

### Why Doesn't MiniWeather Use CUDA?

Because if you've refactored your code to use kernel launching (i.e., CUDA), you should really be using a C++ portability framework. The code is basically identical, but it can run on many different backends from a single source.

### Why Doesn't MiniWeather Use Kokkos or RAJA?

I chose not to use the mainline C++ portability frameworks for two main reasons.

1. It's easier to compile and managed things with a C++ performance portability layer that's < 3K lines of code long, hence: [YAKL (Yet Another Kernel Launcher)](github.com/mrnorman/YAKL). 
2. Kokkos in particular would not play nicely with the rest of the code in the CMake project. Likely if a Kokkos version is added, it will need to be a completely separate project and directory.
3. With `YAKL.h` and `Array.h`, you can see for yourself what's going on when we launch kernels using `parallel_for` on different hardware backends.

# Numerical Experiments

A number of numerical experiments are in the code for you to play around with. You can set these by changing the `data_spec_int` variable. 

## Rising Thermal

```
data_spec_int = DATA_SPEC_THERML
sim_time = 1000
```

This simulates a rising thermal in a neutral atmosphere, which will look something like a “mushroom” cloud (without all of the violence).

Potential Temperature after 500 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/thermal_pt_0500.png" width=400/>

Potential Temperature after 1,000 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/thermal_pt_1000.png" width=400/>

## Colliding Thermals

```
data_spec_int = DATA_SPEC_COLLISION
sim_time = 700
```

This is similar to the rising thermal test case except with a cold bubble at the model top colliding with a warm bubble at the model bottom to produce some cool looking eddies.

Potential Temperature after 200 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/collision_pt_0200.png" width=400/>

Potential Temperature after 400 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/collision_pt_0400.png" width=400/>

Potential Temperature after 700 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/collision_pt_0700.png" width=400/>

## Mountain Gravity Waves

```
data_spec_int = DATA_SPEC_MOUNTAIN
sim_time = 1500
```

This test cases passes a horizontal wind over a faked mountain at the model bottom in a stable atmosphere to generate a train of stationary gravity waves across the model domain.

Potential Temperature after 400 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/mountain_pt_0400.png" width=400/>

Potential Temperature after 1,300 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/mountain_pt_1300.png" width=400/>

## Density Current

```
data_spec_int = DATA_SPEC_DENSITY_CURRENT
sim_time = 600
```

This test case creates a neutrally stratified atmosphere with a strong cold bubble in the middle of the domain that crashes into the ground to give the feel of a weather front (more of a downburst, I suppose).

Potential Temperature after 200 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/density_current_pt_0200.png" width=400/>

Potential Temperature after 600 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/density_current_pt_0600.png" width=400/>

## Injection

```
data_spec_int = DATA_SPEC_INJECTION
sim_time = 1200
```

A narrow jet of fast and slightly cold wind is injected into a balanced, neutral atmosphere at rest from the left domain near the model top. This has nothing to do with atmospheric flows. It's just here for looks. 

Potential Temperature after 300 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/injection_pt_0300.png" width=400/>

Potential Temperature after 1,000 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/injection_pt_1000.png" width=400/>

# Physics, PDEs, and Numerical Approximations

While the numerical approximations in this code are certainly cheap and dirty, they are a fast and easy way to get the job done in a relatively small amount of code. For instance, on 16 K20x GPUs, you can perform a "colliding thermals” simulation with 5 million grid cells (3200 x 1600) in just a minute or two.

## The 2-D Euler Equations

This app simulates the 2-D inviscid Euler equations for stratified fluid dynamics, which are defined as follows:

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\frac{\partial}{\partial&space;t}\left[\begin{array}{c}&space;\rho\\&space;\rho&space;u\\&space;\rho&space;w\\&space;\rho\theta&space;\end{array}\right]&plus;\frac{\partial}{\partial&space;x}\left[\begin{array}{c}&space;\rho&space;u\\&space;\rho&space;u^{2}&plus;p\\&space;\rho&space;uw\\&space;\rho&space;u\theta&space;\end{array}\right]&plus;\frac{\partial}{\partial&space;z}\left[\begin{array}{c}&space;\rho&space;w\\&space;\rho&space;wu\\&space;\rho&space;w^{2}&plus;p\\&space;\rho&space;w\theta&space;\end{array}\right]=\left[\begin{array}{c}&space;0\\&space;0\\&space;-\rho&space;g\\&space;0&space;\end{array}\right]" title="\large \frac{\partial}{\partial t}\left[\begin{array}{c} \rho\\ \rho u\\ \rho w\\ \rho\theta \end{array}\right]+\frac{\partial}{\partial x}\left[\begin{array}{c} \rho u\\ \rho u^{2}+p\\ \rho uw\\ \rho u\theta \end{array}\right]+\frac{\partial}{\partial z}\left[\begin{array}{c} \rho w\\ \rho wu\\ \rho w^{2}+p\\ \rho w\theta \end{array}\right]=\left[\begin{array}{c} 0\\ 0\\ -\rho g\\ 0 \end{array}\right]" />

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\rho_{H}=-\frac{1}{g}\frac{\partial&space;p}{\partial&space;z}" title="\large \rho_{H}=-\frac{1}{g}\frac{\partial p}{\partial z}" />

where <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\rho" title="\large \rho" /> is density, u, and w are winds in the x-, and z-directions, respectively, <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\theta" title="\large \theta" /> is potential temperature related to temperature, T, by <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\theta=T\left(P_{0}/P\right)^{R_{d}/c_{p}}" title="\large \theta=T\left(P_{0}/P\right)^{R_{d}/c_{p}}" />,<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;P_{0}=10^{5}\,\text{Pa}" title="\large P_{0}=10^{5}\,\text{Pa}" /> is the surface pressure, g=9.8<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\text{\,\&space;m}\,\mbox{s}^{-2}" title="\large \text{\,\ m}\,\mbox{s}^{-2}" /> is acceleration due to gravity,<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;p=C_{0}\left(\rho\theta\right)^{\gamma}" title="\large p=C_{0}\left(\rho\theta\right)^{\gamma}" /> is the pressure as determined by an alternative form of the ideal gas equation of state,<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;C_{0}=R_{d}^{\gamma}p_{0}^{-R_{d}/c_{v}}" title="\large C_{0}=R_{d}^{\gamma}p_{0}^{-R_{d}/c_{v}}" />, <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;R_{d}=287\,\mbox{J}\,\mbox{kg}^{-1}\,\mbox{K}^{-1}" title="\large R_{d}=287\,\mbox{J}\,\mbox{kg}^{-1}\,\mbox{K}^{-1}" /> is the dry gas constant, <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\gamma=c_{p}/c_{v}" title="\large \gamma=c_{p}/c_{v}" />, <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;c_{p}=1004\,\mbox{J}\,\mbox{kg}^{-1}\,\mbox{K}^{-1}" title="\large c_{p}=1004\,\mbox{J}\,\mbox{kg}^{-1}\,\mbox{K}^{-1}" /> is specific heat at constant pressure, and <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;c_{v}=717\,\mbox{J}\,\mbox{kg}^{-1}\,\mbox{K}^{-1}" title="\large c_{v}=717\,\mbox{J}\,\mbox{kg}^{-1}\,\mbox{K}^{-1}" /> is specific heat at constant volume. This can be cast in a more convenient form as:

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\frac{\partial\mathbf{q}}{\partial&space;t}&plus;\frac{\partial\mathbf{f}}{\partial&space;x}&plus;\frac{\partial\mathbf{h}}{\partial&space;z}=\mathbf{s}" title="\large \frac{\partial\mathbf{q}}{\partial t}+\frac{\partial\mathbf{f}}{\partial x}+\frac{\partial\mathbf{h}}{\partial z}=\mathbf{s}" />

where a bold font represents a vector quantity.

## Maintaining Hydrostatic Balance

The flows this code simulates are relatively small perturbations off of a “hydrostatic” balance, which balances gravity with a difference in pressure:

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\frac{dp}{dz}=-\rho&space;g" title="\large \frac{dp}{dz}=-\rho g" />

Because small violations of this balance lead to significant noise in the vertical momentum, it's best not to try to directly reconstruct this balance but rather to only reconstruct the perturbations. Therefore, hydrostasis is subtracted from the equations to give:

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\frac{\partial}{\partial&space;t}\left[\begin{array}{c}&space;\rho^{\prime}\\&space;\rho&space;u\\&space;\rho&space;w\\&space;\left(\rho\theta\right)^{\prime}&space;\end{array}\right]&plus;\frac{\partial}{\partial&space;x}\left[\begin{array}{c}&space;\rho&space;u\\&space;\rho&space;u^{2}&plus;p\\&space;\rho&space;uw\\&space;\rho&space;u\theta&space;\end{array}\right]&plus;\frac{\partial}{\partial&space;z}\left[\begin{array}{c}&space;\rho&space;w\\&space;\rho&space;wu\\&space;\rho&space;w^{2}&plus;p^{\prime}\\&space;\rho&space;w\theta&space;\end{array}\right]=\left[\begin{array}{c}&space;0\\&space;0\\&space;-\rho^{\prime}g\\&space;0&space;\end{array}\right]" title="\large \frac{\partial}{\partial t}\left[\begin{array}{c} \rho^{\prime}\\ \rho u\\ \rho w\\ \left(\rho\theta\right)^{\prime} \end{array}\right]+\frac{\partial}{\partial x}\left[\begin{array}{c} \rho u\\ \rho u^{2}+p\\ \rho uw\\ \rho u\theta \end{array}\right]+\frac{\partial}{\partial z}\left[\begin{array}{c} \rho w\\ \rho wu\\ \rho w^{2}+p^{\prime}\\ \rho w\theta \end{array}\right]=\left[\begin{array}{c} 0\\ 0\\ -\rho^{\prime}g\\ 0 \end{array}\right]" />

where a “prime” quantity represents that variable with the hydrostatic background state subtracted off (not a spatial derivative).

## Dimensional Splitting

This equation is solved using dimensional splitting for simplicity and speed. The equations are split into x- and z-direction solves that are, respectively:

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;x:\,\,\,\,\,\,\,\,\,\,\frac{\partial\mathbf{q}}{\partial&space;t}&plus;\frac{\partial\mathbf{f}}{\partial&space;x}=\mathbf{0}" title="\large x:\,\,\,\,\,\,\,\,\,\,\frac{\partial\mathbf{q}}{\partial t}+\frac{\partial\mathbf{f}}{\partial x}=\mathbf{0}" />

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;z:\,\,\,\,\,\,\,\,\,\,\frac{\partial\mathbf{q}}{\partial&space;t}&plus;\frac{\partial\mathbf{h}}{\partial&space;x}=\mathbf{s}" title="\large z:\,\,\,\,\,\,\,\,\,\,\frac{\partial\mathbf{q}}{\partial t}+\frac{\partial\mathbf{h}}{\partial x}=\mathbf{s}" />

Each time step, the order in which the dimensions are solved is reversed, giving second-order accuracy overall. 

## Finite-Volume Spatial Discretization

A Finite-Volume discretization is used in which the PDE in a given dimension is integrated over a cell domain, <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\Omega_{i}\in\left[x_{i-1/2},x_{i&plus;1/2}\right]" title="\large \Omega_{i}\in\left[x_{i-1/2},x_{i+1/2}\right]" />, where <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;x_{i\pm1/2}=x_{i}\pm\Delta&space;x" title="\large x_{i\pm1/2}=x_{i}\pm\Delta x" />, <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;x_{i}" title="\large x_{i}" /> is the cell center, and <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\Delta&space;x" title="\large \Delta x" /> is the width of the cell. The integration is the same in the z-direction. Using the Gauss divergence theorem, this turns the equation into (using the z-direction as an example):

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\frac{\partial\overline{\mathbf{q}}_{i,k}}{\partial&space;t}=-\frac{\mathbf{h}_{i,k&plus;1/2}-\mathbf{h}_{i,k-1/2}}{\Delta&space;z}&plus;\overline{\mathbf{s}}_{i,k}" title="\large \frac{\partial\overline{\mathbf{q}}_{i,k}}{\partial t}=-\frac{\mathbf{h}_{i,k+1/2}-\mathbf{h}_{i,k-1/2}}{\Delta z}+\overline{\mathbf{s}}_{i,k}" />

where <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\overline{\mathbf{q}}_{i,k}" title="\large \overline{\mathbf{q}}_{i,k}" /> and <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\overline{\mathbf{s}}_{i,k}" title="\large \overline{\mathbf{s}}_{i,k}" /> are the cell-average of the fluid state and source term over the cell of index `i,k`.

To compute the update one needs the flux vector at the cell interfaces and the cell-averaged source term. To compute the flux vector at interfaces, fourth-order-accurate polynomial interpolation is used using the four cell averages surrounding the cell interface in question.

## Runge-Kutta Time Integration

So far the PDEs have been discretized in space but are still continuous in time. To integrate in time, we use a simple three-stage, linearly third-order-accurate Runge-Kutta integrator. It is solved as follows:

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\mathbf{q}^{\star}=\mathbf{q}^{n}&plus;\frac{\Delta&space;t}{3}RHS\left(\mathbf{q}^{n}\right)" title="\large \mathbf{q}^{\star}=\mathbf{q}^{n}+\frac{\Delta t}{3}RHS\left(\mathbf{q}^{n}\right)" />

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\mathbf{q}^{\star\star}=\mathbf{q}^{n}&plus;\frac{\Delta&space;t}{2}RHS\left(\mathbf{q}^{\star}\right)" title="\large \mathbf{q}^{\star\star}=\mathbf{q}^{n}+\frac{\Delta t}{2}RHS\left(\mathbf{q}^{\star}\right)" />

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\mathbf{q}^{n&plus;1}=\mathbf{q}^{n}&plus;\Delta&space;tRHS\left(\mathbf{q}^{\star\star}\right)" title="\large \mathbf{q}^{n+1}=\mathbf{q}^{n}+\Delta tRHS\left(\mathbf{q}^{\star\star}\right)" />

When it comes to time step stability, I simply assume the maximum speed of propagation is 450\,\text{m}\,\text{s}^{-1}, which basically means that the maximum wind speed is assumed to be 100\,\text{m}\,\text{s}^{-1}, which is a safe assumption. I set the CFL value to 1.5 for this code.

## Hyper-viscosity

The centered fourth-order discretization is unstable for non-linear equations and requires extra dissipation to damp out small-wavelength energy that would otherwise blow up the simulation. This damping is accomplished with a scale-selective fourth-order so-called “hyper”-viscosity that is defined as:

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\frac{\partial\mathbf{q}}{\partial&space;t}&plus;\frac{\partial}{\partial&space;x}\left(-\kappa\frac{\partial^{3}\mathbf{q}}{\partial&space;x^{3}}\right)=\mathbf{0}" title="\large \frac{\partial\mathbf{q}}{\partial t}+\frac{\partial}{\partial x}\left(-\kappa\frac{\partial^{3}\mathbf{q}}{\partial x^{3}}\right)=\mathbf{0}" />

and this is also solved with the Finite-Volume method just like above. The hyperviscosity constant is defined as:

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\kappa=-\beta\left(\Delta&space;x\right)^{4}2^{-4}\left(\Delta&space;t\right)^{-1}" title="\large \kappa=-\beta\left(\Delta x\right)^{4}2^{-4}\left(\Delta t\right)^{-1}" />

where <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\beta\in\left[0,1\right]" title="\large \beta\in\left[0,1\right]" /> is a user-defined parameter to control the strength of the diffusion, where a higher value gives more diffusion. The parameter <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\beta" title="\large \beta" /> is not sensitive to the grid spacing, and it seems that <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\beta=0.25" title="\large \beta=0.25" /> is generally enough to get rid of <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;2\Delta&space;x" title="\large 2\Delta x" /> noise contamination.

# MiniWeather Model Scaling Details

If you want to do scaling studies with miniWeather, this section will be important to make sure you're doing an apples-to-apples comparison.

* `sim_time`: The `sim_time` parameter does not mean the wall time it takes to simulate but rather refers amount of model time simulated. As you increase `sim_time`, you should expect the walltime to increase linearly.
* `nx_glob, nz_glob`: As a rule, it's easiest if you always keep `nx_glob = nz_glob * 2` since the domain is always 20km x 10km in the x- and z-directions. As you increase `nx_glob` (and proportionally `nz_glob`) by some factor `f`, the time step automatically reduced by that same factor, `f`. Therefore, increasing `nx_glob` by 2x leads to 8x more work that needs to be done. Thus, with the same amount of parallelism, you should expect a 2x increase in `nx_glob` and `nz_glob` to increase the walltime by 8x (neglecting parallel overhead concerns).
  * More precisely, the time step is directly proportional to the minimum grid spacing. The x- and y-direction grid spacings are: `dx=20km/nx_glob` and `dz=10km/nz_glob`. So as you decrease the minimum grid spacing (by increasing `nx_glob` and/or `nz_glob`), you proportionally decrease the size of the time step and therefore proportionally increase the number of time steps you need to complete the simulation (thus proportionally increasing the expected walltime).
* The larger the problem size, `nx_glob` and `nz_glob`, the lower the relative parallel overheads will be. You can get to a point where there isn't enough work on the accelerator to keep it busy and / or enough local work to amortize parallel overheads. At this point, you'll need to increase the problem size to see better scaling. This is a typical [Amdahl's Law](https://en.wikipedia.org/wiki/Amdahl%27s_law) situation.

Remember that you can control each of these parameters through the CMake configure.

# Checking for Correctness

## Domain-Integrated Mass and Total Energy

There are two main ways to check for correctness. The easiest is to look at the domain-integrated mass and total energy printed out at the end of the simulation.

### Mass Change

In all cases for Fortran and C, the relative mass change printed out should be at machine precision (magnitude `1.e-13` or lower just to be flexible with reduced precision optimizations). If the mass changes more than this, you've introduced a bug.

For the C++ code, which uses single precision, the relative mass displayed at the end of the run should be of magnitude `1.e-9` or lower.

### Total Energy Change

In order to use total energy to check the answer, you need to set the following parameters:

* `nx_glob`: >= 100
* `nz_glob`: >= 50 and exactly half of `nx_glob`
* `sim_time = 400`
* `data_spec_int = DATA_SPEC_THERMAL`

Also, it is assumed you have not changed any other default parameters such as `xlen` and `zlen`

From there, you can scale up to any problem size or node count you wish. The relative change in total energy should always be negative, and the magnitude should always be less than `4.5e-5`. If the magnitude is larger than this, or if the value is positive, then you have introduced a bug. As you increase the problem size, the energy is always better conserved. These total energy change values are valid for single precision in C++ as well.


## NetCDF Files

Your other option is to create two baseline NetCDF files whose answers you trust: (1) with `-O0` optimizations; and (2) with `-O3` optimizations. Then, you can use the following python script to do a 3-way diff between the two baselines and the refactored code. The refactored diff should be of the same order of magnitude as the baseline compiler optimization diffs. Note that if you run for too long, non-linear chaotic amplification of the initially small differences will eventually be come too large to make for a useful comparison, so try to limit the simulation time to, say, 400 seconds or less.

The reason you have to go to all of this trouble is because of chaotic amplification of initially small differences (the same reason you can't predict weather reliably past a few days). Therefore, you can't compare snapshots to machine precision tolerance.

<details><summary>Click here to expand python script</summary>
 <p>
  
```python
import netCDF4
import sys
import numpy as np

#######################################################################################
#######################################################################################
##
## nccmp3.py: A simple python-based NetCDF 3-way comparison tool. The purpose of this
## is to show whether files 2 and 3 have more of a difference than files 1 and 2. The
## specific purpose is to compare refactored differences against presumed bit-level
## differences (like -O0 and -O3 compiler flags). This prints the relative 2-norm of
## the absolute differences between files 1 & 2 and files 2 & 3, as well as the ratio
## of the relative 2-norms between the 2-3 comparison and the 1-2 comparison.
##
## python nccmp.py file1.nc file2.nc file3.nc
##
#######################################################################################
#######################################################################################

#Complain if there aren't two arguments
if (len(sys.argv) < 4) :
  print("Usage: python nccmp.py file1.nc file2.nc")
  sys.exit(1)

#Open the files
nc1 = netCDF4.Dataset(sys.argv[1])
nc2 = netCDF4.Dataset(sys.argv[2])
nc3 = netCDF4.Dataset(sys.argv[3])

#Print column headers
print("Var Name".ljust(20)+":  "+"|1-2|".ljust(20)+"  ,  "+"|2-3|".ljust(20)+"  ,  "+"|2-3|/|1-2|")

#Loop through all variables
for v in nc1.variables.keys() :
  #Only compare floats
  if (nc2.variables[v].dtype == np.float64 or nc2.variables[v].dtype == np.float32) :
    #Grab the variables
    a1 = nc1.variables[v][:]
    a2 = nc2.variables[v][:]
    a3 = nc3.variables[v][:]
    #Compute the absolute difference vectors
    adiff12 = abs(a2-a1)
    adiff23 = abs(a3-a2)

    #Compute relative 2-norm between files 1 & 2 and files 2 & 3
    norm12 = np.sum( adiff12**2 )
    norm23 = np.sum( adiff23**2 )
    #Assume file 1 is "truth" for the normalization
    norm_denom = np.sum( a1**2 )
    #Only normalize if this denominator is != 0
    if (norm_denom != 0) :
      norm12 = norm12 / norm_denom
      norm23 = norm23 / norm_denom

    #Compute the ratio between the 2-3 norm and the 1-2 norm
    normRatio = norm23
    #If the denom is != 0, then go ahead and compute the ratio
    if (norm12 != 0) :
      normRatio = norm23 / norm12
    else :
      #If they're both zero, then just give a ratio of "1", showing they are the same
      if (norm23 == 0) :
        normRatio = 1
      #If the 1-2 norm is zero but the 2-3 norm is not, give a very large number so the user is informed
      else :
        normRatio = 1e50

    #Only print ratios that are > 2, meaning 2-3 diff is >2x more than the 1-2 diff.
    #In the future, this should be added as a command line parameter for the user to choose.
    if (normRatio > 2) :
      print(v.ljust(20)+":  %20.10e  ,  %20.10e  ,  %20.10e"%(norm12,norm23,norm23/norm12))
```
</p>
</details>
