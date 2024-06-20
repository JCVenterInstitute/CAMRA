FROM debian:stable-slim

# Set environment variables for non-interactive installation
ENV DEBIAN_FRONTEND=noninteractive

# Update, upgrade, install dependencies, and clean up 
RUN apt update && \
    apt upgrade -y && \
    apt install -y wget libssl-dev cpanminus gdebi-core && \
    cpanm Term::ReadLine::Stub && \
    wget https://github.com/BV-BRC/BV-BRC-CLI/releases/download/1.040/bvbrc-cli-1.040.deb && \
    gdebi -n bvbrc-cli-1.040.deb && \
    rm bvbrc-cli-1.040.deb && \
    apt clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


WORKDIR /data

# Specify the default command
CMD ["bash"]
