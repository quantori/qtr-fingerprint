# qtr-fingerprint

// TODO: translate to english

## Description

// TODO: The program consists of 3 targets... Program supports 2 DatabaseTypes

## Структура проекта

### Основные папки

* [cpp/preprocessing](cpp/preprocessing) - код программы предобработки данных.
* [cpp/build_db](cpp/build_db) - код программы построения поискового индекса по preprocessed data.
* [cpp/run_db](cpp/run_db) - код программы для выполнения запросов к поисковому индексу.
* [cpp/qtrlib](cpp/qtrlib) - основная библиотека со всем необходимым для Qtr алгоритма substructure search.

### Вспомогательные папки, не используемые в основной части программы

* [python](python) - экспериментальный код на языке python.
* [cpp/playground](cpp/playground) - экспериментальный код на языке c++.
* [notebooks](notebooks) - jupyter notebooks используемые для проведения исследований.

## Консольные приложения

Проект реализует несколько консольных приложений для работы с датасетом молекул.

Есть следующие типы баз данных:

* `QtrRam` - Qtr поисковый индекс, выгружающий данные в оперативную память при работе с базой.
* `BingoNoSQL` - база данных, использующая поисковый
  индекс [bingo NoSQL](https://lifescience.opensource.epam.com/bingo/bingo-nosql.html).
* `QtrDrive` - **(Not implemented yet)** Qtr поисковый индекс, хранящий данные на жёстком диске, используя memory
  mapping. Что позволяет алгоритму работать, в том числе в условиях ограниченного объёма оперативной памяти.

Также Qtr базы данных поддерживает поиск с учётом свойств молекул.
А именно, поддерживаются следующие свойства:

- `PUBCHEM_COMPONENT_COUNT`,
- `PUBCHEM_XLOGP3`,
- `PUBCHEM_ATOM_UDEF_STEREO_COUNT`,
- `PUBCHEM_HEAVY_ATOM_COUNT`,
- `PUBCHEM_CACTVS_TAUTO_COUNT`,
- `PUBCHEM_ISOTOPIC_ATOM_COUNT`,
- `PUBCHEM_CACTVS_HBOND_DONOR`,
- `PUBCHEM_CACTVS_ROTATABLE_BOND`,
- `PUBCHEM_MONOISOTOPIC_WEIGHT`,
- `PUBCHEM_CACTVS_HBOND_ACCEPTOR`,
- `PUBCHEM_ATOM_DEF_STEREO_COUNT`,
- `PUBCHEM_COMPOUND_CID`,
- `PUBCHEM_MOLECULAR_WEIGHT`,
- `PUBCHEM_BOND_DEF_STEREO_COUNT`,
- `PUBCHEM_TOTAL_CHARGE`,
- `PUBCHEM_EXACT_MASS`,
- `PUBCHEM_CACTVS_COMPLEXITY`,
- `PUBCHEM_BOND_UDEF_STEREO_COUNT`,
- `PUBCHEM_CACTVS_TPSA`,
- `PUBCHEM_COMPOUND_CANONICALIZED`

### preprocessing

_Предобработка данных: построение fingerprints, нумерация молекул и т.д._

* `--properties` - 1 if properties should be concerned during preprocessing. 0 otherwise
* `--sourceDir` - path to directory where source files are stored
* `--destDir` - destination directory where preprocessed files should be stored
* `--preprocessingType` - source files type. (`SDF` or `CSV`)
    - `CSV` - если выбран данный тип препроцессинга, то из `sourceDir` интерпретируются как таблицы формата `.csv`.
      В первом столбце записан `id` молекулы, во втором задана сама молекула в формате `SMILES`. Если установлен
      флаг `--properties`, то в колонках 2-22 должны быть заданы значения свойств молекулы.
    - `SDF` - deprecated, может работать некорректно.
* `--targetFilesType` - deprecated.

### build_db

_Построение поискового индекса._

* `--dbType` - `QtrDrive`/`QtrRam`/`BingoNoSQL`
* `--properties` - 1 if properties should be concerned during preprocessing. 0 otherwise.
  (только для `QtrRam`).
* `--dbName` - название поискового индекса
* `--sourceDir` - папка с предобработанным программой [preprocessing](#preprocessing) файлами.
* `--destDirs` - папки, в которые стоит сохранять файлы поискового индекса
  (для `BingoNoSQL` и `QtrDrive` будет использоваться только первый переданный путь, для `QtrRam`
  есть смысл указывать более одной папки, если они находятся физически на разных жёстких дисках.
  Тогда параллельность по жёстким дискам позволит быстрее работать с индексом).
* `--otherDataDir` - папка, в которую стоит сохранять файлы поискового индекса, которые нельзя хранить параллельно
  (только для `QtrRam`).
* `--parallelizeDepth` - глубина, начиная с которой поддеревья поискового индекса начинают строиться параллельно.
  (только для `QtrRam` и `QtrDrive`).
* `treeDepth` - дерево какой глубины должно быть построено (актуально только для `QtrRam` и `QtrDrive`).

### run_db

// TODO: describe run_db target and its arguments

### Benchmarking

// TODO: show steps to repeat experiment's results

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

## Usage

// TODO: describe all targets and its arguments

## Run tests

// TODO: fix instruction

1. `cd ./qtr-fingerprint/cpp/build`
2. `./bin/tests`

If you want to run slow tests (tests on big data) ypu should specify big_data_dir_path. You could find big data files on
s3 https://s3.console.aws.amazon.com/s3/buckets/sfo-general?region=us-east-2&tab=objects

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
