FROM ubuntu:20.04

# Set non-interactive mode for apt-get
ENV DEBIAN_FRONTEND=noninteractive

# Install required dependencies
RUN apt-get update && \
    apt-get install -y python2.7 python3 sudo wget ed && \
    ln -s /usr/bin/python2.7 /usr/bin/python2 && \
    ln -s /usr/bin/python3 /usr/bin/python && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh && \
bash /miniconda.sh -b -p /opt/conda && \
rm /miniconda.sh

ENV PATH="/opt/conda/bin:$PATH"

# Configure Conda and install LEfSe (Python 2.7)
RUN conda update -n base -c defaults conda -y && \
conda create -n lefse_env python=2.7 -y && \
conda install -n lefse_env -c bioconda lefse -y && \
conda clean --all -y

# Install Python 3 environment
RUN conda create -n py3_env python=3.8 -y && \
conda clean --all -y

# Set environment paths
ENV PATH="/opt/conda/envs/lefse_env/bin:/opt/conda/envs/py3_env/bin:$PATH"

# Copy scripts inside the container
COPY Scripts/ ./Scripts/
  
# Set the default shell to bash
SHELL ["/bin/bash", "-c"]

CMD ["bash"]
