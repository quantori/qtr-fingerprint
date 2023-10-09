FROM python:3.9

RUN apt install libxrender1

# Web API
EXPOSE 8000

WORKDIR app
COPY requirements.txt ./
RUN pip install -r requirements.txt

COPY ./qtrshiny ./qtrshiny

ENTRYPOINT shiny run ./qtrshiny/app.py --host 0.0.0.0