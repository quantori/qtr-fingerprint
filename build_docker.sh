#! /bin/bash
docker build -t rdkit_build:local . --progress=plain

docker run -it --mount type=bind,source=.,target=/src rdkit_build:local