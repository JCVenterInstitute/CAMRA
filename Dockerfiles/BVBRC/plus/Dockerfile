# Use the specified Docker image as a base
FROM danylmb/bvbrc:1.040

# Install Python
RUN apt-get update && \
    apt-get install -y python3 python3-pip

RUN apt-get install -y python3-bs4

# Copy custom scripts into the Docker image
COPY ./scripts /bin/

# Ensure the python scripts are executable
RUN chmod +x /bin/*.py

# Set the entrypoint to bash
ENTRYPOINT ["/bin/bash"]