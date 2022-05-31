# qtr-fingerprint

## Requirements

* `CMake 3.13 or higher`
* `libfreetype6-dev` and `libfontconfig1-dev` libs
* `g++ 9.4 or higher`

## Build
1. `git clone https://github.com/quantori/qtr-fingerprint.git`
2. `cd ./qtr-fingerprint/cpp`
3. `git submodule update --init --recursive`
4. `mkdir build && cd build`
5. `cmake ..`
6. `cmake --build .` or `cmake --build . --target tests` to specify build target

## Run tests
1. `cd ./qtr-fingerprint/cpp/build`
2. `./bin/tests`

## Program arguments

You can get help information about your executable just by running 

`./program --help`

You will see the list of available arguments and their descriptions.

And then you can specify each option as following:

`./app --path_to_query=data/query.mol --database_path=data/119697.sdf`

## Log configuration

In order to configure logger, you should add environment variables, 
most common with their default values are listed below:
1. `GLOG_log_dir=`
2. `GLOG_alsologtostderr=false`
3. `GLOG_logtostdout=false`
4. `GLOG_minloglevel=0`, order of levels are: `INFO, WARNING, ERROR, FATAL`
