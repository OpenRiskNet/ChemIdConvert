FROM continuumio/miniconda3
MAINTAINER Daniel Bachler <daniel@douglasconnect.com>

RUN apt-get update && apt-get install -y libfontconfig1 libxrender1

RUN conda config --add channels rdkit
RUN conda config --add channels conda-forge
RUN conda config --add channels mcs07
RUN conda install -y rdkit

COPY python-flask-server/requirements.txt /python-flask-server/requirements.txt
RUN conda install --yes --file /python-flask-server/requirements.txt
COPY python-flask-server/ /python-flask-server/

EXPOSE 8080
CMD ["python", "/python-flask-server/app.py"]
