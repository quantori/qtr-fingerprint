# How to build:
This command will create docker image prepared to build rdkit and run it (it will take several minutes..)
```./build_docker.sh```
After this you will be in Docker container. Please note, source files is mounted into the container, so after work your repo will be polluted with compiled files.
Right after execution you should get into docker container and see:
```(base) root@09af558870d3:/src/cpp/third_party/rdkit#``` 

Then in command line type:
```cmake -DPy_ENABLE_SHARED=1   -DRDK_INSTALL_INTREE=ON   -DRDK_INSTALL_STATIC_LIBS=OFF   -DRDK_BUILD_CPP_TESTS=ON -DRDK_BUILD_INCHI_SUPPORT=ON -DRDKIT_RDINCHILIB_BUILD=ON -DPYTHON_NUMPY_INCLUDE_PATH="$(python -c 'import numpy ; print(numpy.get_include())')"   -DBOOST_ROOT="$CONDA_PREFIX" .```
Then type
```make -j3``` 
Here -j3 means make uses 3 cores to paralelize, it speedups build. You may play with it, but left 1 core free for OS.
Then
```cd Docs/Book/C++Examples```
```cmake -DBoost_NO_BOOST_CMAKE=ON .```
```make -j3```
Try to execute examples::
```./example1```