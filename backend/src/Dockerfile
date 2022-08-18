FROM continuumio/miniconda3
RUN conda install -c conda-forge rdkit

ADD requirements.in /app/src/
RUN python3 -m pip install -r /app/src/requirements.in

# copy source code
COPY . /app/src