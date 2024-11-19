FROM debian:stable-slim

# This cannot be an environmental variable. As an arg, it is only used for the build.
# If it is set as an environmental variable, it affects runtime behavior and breaks Terra.
ARG DEBIAN_FRONTEND=noninteractive

# Update, upgrade, install dependencies, and clean up 
RUN apt update && \
    apt upgrade -y && \
    apt install -y wget libssl-dev cpanminus gdebi-core python3-pandas python3-numpy && \
    cpanm Term::ReadLine::Stub && \
    wget https://github.com/BV-BRC/BV-BRC-CLI/releases/download/1.040/bvbrc-cli-1.040.deb && \
    gdebi -n bvbrc-cli-1.040.deb && \
    rm bvbrc-cli-1.040.deb && \
    apt clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
    
# Install Python
RUN apt-get update
RUN apt-get install -y python3 python3-pip
RUN apt-get install -y python3-bs4

# Copy custom scripts into the Docker image
COPY ./scripts /bin/

# Ensure the python scripts are executable
RUN chmod +x /bin/*.py

WORKDIR /data

# Specify the default command
CMD ["bash"]
