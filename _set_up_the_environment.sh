#!/bin/bash

mkdir docker
docker build -t ccribioinf/dockrstudio:4.2.0-TEST-DELETE_ME -f Rstudio-4.2.0-java.Dockerfile ./docker/

