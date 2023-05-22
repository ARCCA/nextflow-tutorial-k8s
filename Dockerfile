# work from latest LTS ubuntu release
FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

# some container metadata
LABEL "MANTAINER"="ARCCA,Jose Munoz"
LABEL "EMAIL"="munozcriollojj@cardiff.ac.uk"

# The WORKDIR instruction sets the working directory for any RUN, CMD,
# ENTRYPOINT, COPY and ADD instructions that follow it in the Dockerfile.
# If the WORKDIR doesn’t exist, it will be created even if it’s not used 
# in any subsequent Dockerfile instruction.
WORKDIR /app

# run update and install necessary tools from package manager
RUN apt-get update -y && apt-get install -y \
    default-jre cpanminus wget \
    unzip autoconf libbz2-dev liblzma-dev libncurses5-dev cmake \
    libjsoncpp-dev r-base

# Install R - https://cloud.r-project.org/bin/linux/ubuntu/
# update indices
RUN apt update -qq -y
# install two helper packages we need
RUN apt install -y --no-install-recommends software-properties-common dirmngr
# add the signing key (by Michael Rutter) for these repos
# To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc 
# Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc |\
    tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu \
    $(lsb_release -cs)-cran40/"
RUN apt install -y --no-install-recommends r-base


# install python3 and pip
RUN apt-get install -y python3-pip

# Install multiqc
RUN pip install multiqc==1.14
# Install cutadapt
RUN pip install cutadapt==2.8

# Install fastqc
RUN mkdir -p /opt
COPY fastqc_v0.12.1.zip /opt
RUN cd /opt && unzip fastqc_v0.12.1.zip && rm fastqc_v0.12.1.zip
ENV PATH="${PATH}:/opt/FastQC"

# Install Trim Galore
ADD trim_galore.tar.gz /opt
ENV PATH="${PATH}:/opt/TrimGalore-0.6.10"

# Install Picard
RUN mkdir -p /opt/picard
ADD picard.jar /opt/picard
ENV PATH="${PATH}:/opt/picard"

# Install STAR
RUN mkdir -p /opt/star
COPY STAR_2.7.10b_alpha_230301_Linux_x86_64_static.zip /opt/star
RUN cd /opt/star && \
    unzip STAR_2.7.10b_alpha_230301_Linux_x86_64_static.zip && \
    rm STAR_2.7.10b_alpha_230301_Linux_x86_64_static.zip
ENV PATH="${PATH}:/opt/star"

# Install HTSLIB
ADD htslib-1.9.tar.gz /tmp
RUN cd /tmp/htslib-1.9 && autoreconf -i && \
    ./configure --prefix=/opt/htslib && make && make install && \
    rm -rf /tmp/htslib-1.9
ENV PATH="${PATH}:/opt/htslib/bin"
ENV LD_LIBRARY_PATH="/opt/htslib/lib:${LD_LIBRARY_PATH}"


# Install Samtools
ADD samtools-1.9.tar.gz /tmp
RUN cd /tmp/samtools-1.9 && autoheader && autoconf -Wno-syntax && \
    ./configure --prefix=/opt/samtools --with-htslib=/opt/htslib && \
    make && make install && rm -rf /tmp/samtools-1.9
ENV PATH="${PATH}:/opt/samtools/bin"

# Install Bamtools
ADD bamtools_v2.5.2.tar.gz /tmp
RUN cd /tmp/bamtools-2.5.2 && mkdir build && cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/opt/bamtools .. && make && make install && \
    rm -rf /tmp/bamtools-2.5.2
ENV PATH="${PATH}:/opt/bamtools/bin"

# Install Subread
ADD subread-2.0.0-Linux-x86_64.tar.gz /opt
ENV PATH="${PATH}:/opt/subread-2.0.0-Linux-x86_64/bin"

# Install Nextflow
COPY nextflow-21.10.6-all /tmp
RUN mkdir -p /opt/nextflow && cd /opt/nextflow && \
    wget -qO- https://get.nextflow.io | bash && \
    chmod +x /opt/nextflow/nextflow
ENV PATH="${PATH}:/opt/nextflow"

