FROM python:3.10
LABEL maintainer="crashfrog@gmail.com"

WORKDIR /src
COPY ./ ./
RUN pip install -e --no-deps .