# FORTRAN X-TFC

Extreme Theory of Functional Connections (X-TFC) algorithm implemented in FORTRAN to solve Ordinary Differential Equations (ODEs) as Initial Value Problems with minimal dependencies. It is applied at a simple example, for solving other problems modify `problem.f90`. This script contains the parameter definitions, the ODEs, and the Jacobian of the loss function. If you want to solve chemical kinetics equations look at the [WinNet](https://github.com/neuralXray/WinNet) project for a widely tested and generalized program (minimal input required, no code modifications, you do not have to worry about the Jacobian).

Modules:

* `problem`. Problem definition: parameters and equations.

* `xtf`. X-TFC method.

* `utils`. Utilities for saving the solution and the trained X-TFC network to text files.

* `backward_euler`. Implicit backward Euler method.

The problem is solved in a log spaced training grid.

`main.f90`. FORTRAN program that calls the routines defined in the above modules to solve the problem with both solvers. It increases the resolution of the solution provided by the X-TFC making use of the analytical representation obtained.

`plot.py`. Python script to plot the results.


## Installation

Commands for Debian-based Linux distributions.

Intel® oneAPI Math Kernel Library for FORTRAN.

```
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt update
sudo apt install intel-oneapi-common-vars intel-oneapi-compiler-fortran intel-oneapi-mkl intel-oneapi-mkl-devel
echo 'source /opt/intel/oneapi/setvars.sh' >> ~/.bashrc
source ~/.bashrc
```

Python packages: numpy, matplotlib and scipy. Just for the plot.

```
python3 -m venv venv
source venv/bin/activate
python -m pip install -r requirements.txt
```

Create `/data` directory.

`mkdir data`


## Compilation & Execution

Compile,

```
ifort -c -qmkl -diag-disable=10448 problem.f90 utils.f90 xtfc.f90 backward_euler.f90 && ifort -o main main.f90 -qmkl -diag-disable=10448 utils.o xtfc.o backward_euler.o problem.o
```

Execute,

`./main`

And plot.

`python plot.py`

Compile with extra checks.

```
ifort -c -qmkl -traceback -check all -diag-disable=10448 problem.f90  utils.f90 xtfc.f90 backward_euler.f90
ifort -o main main.f90 -qmkl -traceback -check all -diag-disable=10448 utils.o xtfc.o backward_euler.o problem.o
```

Compile with performance profile.

```
ifort -c -qmkl -pg -diag-disable=10448 problem.f90  utils.f90 xtfc.f90 backward_euler.f90
ifort -o main main.f90 -qmkl -pg -diag-disable=10448 utils.o xtfc.o backward_euler.o problem.o
```

Execute and then explore the profile (this problem is solved too fast, you will not see anything).

```
./main
gprof ./main gmon.out > output.txt
```


## Robertson's problem

$$
	\begin{cases}
        A     & \hspace{.2cm} \stackrel{k_1}{\longleftrightarrow} \hspace{.2cm} & B \\
        B + B & \hspace{.2cm} \stackrel{k_2}{\longleftrightarrow} \hspace{.2cm} & C + B \\
        B + C & \hspace{.2cm} \stackrel{k_3}{\longleftrightarrow} \hspace{.2cm} & A + C \\
	\end{cases}
$$

$$
	\begin{cases}
        k_1 = 4 \times 10^{-2}, \hspace{.2cm} k_2 = 3 \times 10^7, \hspace{.2cm} k_3 = 10^4 \\
        y_A(0) = 1, \hspace{.2cm} y_B(0) = 0, \hspace{.2cm} y_C(0) = 0 \\
        t \in [10^{-5}, 10^5]
	\end{cases}
$$

