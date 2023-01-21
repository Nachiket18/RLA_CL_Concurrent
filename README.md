# RLA_CL_Concurrent
Record Linkage with concurrent programming using openMP

# Setting the parameters for the program

```
cd src

```
1. Open file RLA_parallel.h
2. Go to section USER INPUTS and make necessary modifications

# Compiling
```
cd src

```
```
g++ RLA_parallel.cpp -fopenmp -O3

```
If the compilation throws error then run

```
g++ RLA_parallel.cpp -fopenmp -ltbb -O3

```

# Running the program

```
./RLA_parallel

```

