# select parent image
FROM alpine:latest

# install ANT, OpenJDK, GIT and BASH
RUN apk add --no-cache openjdk11-jre apache-ant git bash

# copy the source tree to container
COPY ./ ./aracne
COPY .git/ ./aracne/.git/

# package application code
WORKDIR /aracne
RUN ant main
WORKDIR /

# add to classpath
ENV CLASSPATH=/aracne/dist/aracne.jar:${CLASSPATH}

# set the startup command to execute the jar
CMD ["java", "aracne.Aracne]
