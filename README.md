# RLA_CL_Concurrent
Record Linkage with concurrent programming

# Installing oneTBB

```
# Do our experiments in /tmp
cd /tmp
# Clone oneTBB repository
git clone https://github.com/oneapi-src/oneTBB.git
cd oneTBB
# Create binary directory for out-of-source build
mkdir build && cd build
# Configure: customize CMAKE_INSTALL_PREFIX and disable TBB_TEST to avoid tests build
cmake -DTBB_TEST=OFF ..
# Build
cmake --build .
# Install
cmake --install .
# Well done! Your installed oneTBB is in /tmp/my_installed_onetbb
```

# Build instructions
```
cd src
cmake build .
make

```

# Running the program

```
./RLA_parallel
```

