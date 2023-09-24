# qtr-fingerprint

## Requirements

* `CMake 3.13 or higher`
* `ninja 1.7.2 or higher`
* `libfreetype6-dev`, `libfontconfig1-dev`, `libasio-dev`, `libgflags-dev` libs

```apt-get install libfreetype6-dev libfontconfig1-dev libasio-dev libgflags-dev```
to install them all
* `g++ 9.4 or higher`

## Build
1. `git clone https://github.com/quantori/qtr-fingerprint.git`
2. `cd ./qtr-fingerprint/cpp`
3. `git submodule update --init --recursive`
4. `mkdir build && cd build`
5. `cmake ..`
6. `cmake --build . -j10` or `cmake --build . --target tests` to specify build target

## Run tests
1. `cd ./qtr-fingerprint/cpp/build`
2. `./bin/tests`

If you want to run slow tests (tests on big data) ypu should specify big_data_dir_path. You could find big data files on s3 https://s3.console.aws.amazon.com/s3/buckets/sfo-general?region=us-east-2&tab=objects

`./bin/tests --big_data_dir_path=<path to dir with big data>`
## Program arguments

You can get help information about your executable just by running 

`./program --help`

You will see the list of available arguments and their descriptions.

And then you can specify each option as following:

`./app --path_to_query=data/query.mol --database_path=data/119697.sdf --search_engine_type=bingo`

## Log configuration

In order to configure logger, you should add environment variables, 
most common with their default values are listed below:
1. `GLOG_log_dir=`
2. `GLOG_alsologtostderr=false`
3. `GLOG_logtostdout=false`
4. `GLOG_minloglevel=0`, order of levels are: `INFO, WARNING, ERROR, FATAL`
