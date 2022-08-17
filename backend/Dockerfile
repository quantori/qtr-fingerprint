FROM continuumio/miniconda3
RUN conda install -c conda-forge rdkit
RUN apt update

WORKDIR /app/backend

ADD requirements.in /app/backend/
RUN python3 -m pip install -r /app/backend/requirements.in

# copy source code
COPY . /app/backend

ENV FLASK_APP app.py

EXPOSE 5000
