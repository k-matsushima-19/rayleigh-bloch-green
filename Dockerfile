FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive

#-----------------
# build-essential
#-----------------
RUN apt-get update && \
    apt-get install -y build-essential g++ \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

#-------
# eigen
#-------
RUN apt-get update && apt-get -y --no-install-recommends install \
 libeigen3-dev=3.3.* \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# copy
COPY . /app/

# clean
RUN rm -rf /usr/local/src/*

# Run
WORKDIR /app
ENTRYPOINT ["/bin/bash", "run.sh"]