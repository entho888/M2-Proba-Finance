FROM gitpod/workspace-python:latest

COPY requirements2.txt ./
RUN python -m pip install --no-cache-dir -r requirements2.txt

# fix non persistent pip packages: https://github.com/gitpod-io/gitpod/issues/7077
ENV PIP_USER=yes
