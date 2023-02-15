FROM ubuntu:latest

# Add an LABEL with the author and description
LABEL author="Simonas Juzenas" \
      description="Docker image containing all required tools for drop-srna-seq pipeline"

# Update the package manager and install wget
RUN apt-get update && apt-get install -y wget

# Download and install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh

# Set the PATH to include Miniconda
ENV PATH="/miniconda/bin:${PATH}"

# Copy the environment.yml file to the container
COPY environment.yml /

# Create the conda environment
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH
ENV PATH /opt/conda/envs/drop-srna-seq-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name drop-srna-seq-1.0dev > drop-srna-seq-1.0dev.yml