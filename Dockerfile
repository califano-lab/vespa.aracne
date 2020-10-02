# select parent image
FROM alpine:latest

# install ANT, OpenJDK and GIT
RUN apk add --no-cache openjdk11-jre apache-ant git

# copy the source tree to container
COPY ./ ./
COPY .git/ ./.git/

# package application code
RUN ant main
 
# set the startup command to execute the jar
CMD ["java", "-jar", "dist/aracne.jar"]
