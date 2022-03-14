FROM ubuntu:18.04
CMD ["bash", "c"]
# Install Python 3.8.2 & Pipenv:
ENV PYENV_ROOT="/root/.pyenv" \
    PIPENV_YES=1 \
    PATH="/root/.pyenv/shims:/root/.pyenv/bin:${PATH}" \
    PIPENV_DONT_LOAD_ENV=1 \
    LANG="en_US.UTF-8" \
    LC_ALL="C.UTF-8" 
RUN apt-get update && apt-get install -yq \
  build-essential git libreadline-dev zlib1g-dev libssl-dev libbz2-dev libsqlite3-dev libffi-dev jq curl liblzma-dev\
  && curl -L https://raw.githubusercontent.com/yyuu/pyenv-installer/master/bin/pyenv-installer | bash \
  && pyenv install 3.8.2 && pyenv global 3.8.2 && pyenv rehash && pip install --upgrade pip
RUN apt-get update && apt-get install -yq curl \
  apt-transport-https ca-certificates software-properties-common \
  && curl -fsSL https://download.docker.com/linux/ubuntu/gpg | apt-key add - \
  && add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu bionic stable" \
  && apt update && apt-cache policy docker-ce \
  && apt install -yq --no-install-recommends docker-ce
RUN apt-get update && apt-get install -yq net-tools netcat iputils-ping
RUN apt-get update && apt-get install -yq vim
# R
ENV DEBIAN_FRONTEND noninteractive
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV R_BASE_VERSION 3.5.3
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      ed less locales vim-tiny wget fonts-texgyre tzdata \
      gsl-bin libgsl-dev \
      zlib1g zlib1g-dev \
    libssl-dev \
    libxml2-dev \
    libsqlite-dev \
    libcurl4-gnutls-dev \
    libxt-dev \
    ncurses-base libncurses5-dev \
    libcairo2-dev \
    pandoc \
    pandoc-citeproc \
    lmodern \
    && rm -rf /var/lib/apt/lists/*
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && locale-gen en_US.utf8 && /usr/sbin/update-locale LANG=en_US.UTF-8
## Use Debian unstable via pinning -- new style via APT::Default-Release
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
  && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
## Now install R and littler, and create a link for littler in /usr/local/bin
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        littler r-cran-littler \
        r-base=${R_BASE_VERSION}-* \
        r-base-dev=${R_BASE_VERSION}-* \
        r-recommended=${R_BASE_VERSION}-* \
    && ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
    && ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
    && ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
    && ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
    && install.r docopt \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds
# R packages for DBGWAS:
RUN R -e "install.packages(c('devtools', 'testthat', 'roxygen2', 'docopt', 'curl', 'Biostrings', 'optparse', 'ape'))"
RUN R -e "devtools::install_version('phangorn', version = '2.2', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('https://raw.githubusercontent.com/sgearle/bugwas/master/build/bugwas_1.0.tar.gz', repos=NULL, type='source')"
RUN apt-get update && apt-get install -yq cmake zlib1g-dev
RUN wget -O ncbi-blast.tar.gz \
    https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.1/ncbi-blast-2.10.1+-x64-linux.tar.gz && \
    tar -xvf ncbi-blast.tar.gz && mv ncbi-blast-2.10.1+ ncbi-blast && rm ncbi-blast.tar.gz 
# Install DBGWAS
RUN git clone --depth=1 --recurse-submodules -j8 -b caldera https://gitlab.com/leoisl/dbgwas.git ./dbgwas && \
    rm ./dbgwas/sample_example/* && rm -rf ./dbgwas/.git
RUN cd /dbgwas && mkdir build \
  && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make \
  && cd DBGWAS/DBGWAS && make install
ENV PATH="/ncbi-blast/bin/:${PATH}"
# Install tools to pre-process ATH data
RUN R -e "install.packages('BiocManager')" && \
    R -e "BiocManager::install(c('KEGGREST', 'DEGraph'))"
RUN mkdir scr && cd /scr && \
    wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200428.zip && \
    unzip plink_linux_x86_64_20200428.zip && rm plink_linux_x86_64_20200428.zip && \
    wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip && \
    unzip snpEff_latest_core.zip && rm snpEff_latest_core.zip
RUN apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get clean && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME
ENV PATH="/scr/:/snpEff/:${PATH}"
# Install python packages for CALDERA
RUN pip3 install --upgrade python-igraph numpy scipy sklearn awscli joblib argparse pandas progressbar
ENV PATH="/usr/local/HDF_Group/HDF5/1.8.18/bin/:/root/.local/bin:/dbgwas/blast/:${PATH}"
# Install R packages for pre and post-processing
RUN R -e "install.packages(c('here', 'dplyr', 'ggplot2', 'stringr', 'igraph'))" && \
  R -e "install.packages(c('tidyr', 'cowplot', 'network', 'ggnetwork', 'openxlsx'))" && \
  R -e "install.packages(c('phylogram', 'dendextend', 'readr', 'purrr'))" && echo "" > .here

# Install pyseer and need tools around
RUN cd ${HOME} && mkdir software && cd software && wget https://github.com/simongog/sdsl-lite/archive/v2.0.3.tar.gz &&  \
  tar -xvf v2.0.3.tar.gz && rm v2.0.3.tar.gz && cd sdsl-lite-2.0.3 && ./install.sh ../
RUN cd /scr && git clone https://github.com/nvalimak/fsm-lite.git && \
  cd fsm-lite && make depend && make
RUN cd /scr && git clone https://github.com/mgalardini/pyseer.git && rm -rf pyseer/.git && \
  cd pyseer && python3 -m pip install -e .
RUN cd /scr && wget https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar && \
  tar -xvf mash-Linux64-v2.3.tar && rm mash-Linux64-v2.3.tar && mv mash-Linux64-v2.3 mash
RUN cd /scr && git clone https://github.com/HIITMetagenomics/dsm-framework.git  && \
  cd dsm-framework && rm -rf .git && make clean && make
# Install Kover and needed tools around
# RUN apt-get update && apt-get install -y --no-install-recommends \
#   make llvm libncursesw5-dev xz-utils tk-dev libxmlsec1-dev && \
#   pyenv install 2.7.18 && pyenv global 2.7.18 && \
#   python -m pip install numpy scipy pandas progressbar cython && \
#   cd /scr && git clone https://github.com/aldro61/kover.git && \
#   cd kover && ./install.sh && \
#   pyenv global 3.8.2
# ENV PATH="/scr/kover/bin/:${PATH}"

COPY ./src /src
RUN python3 -m pip install -e src/CALDERA/. && python3 -m pip install -e src/COIN/.
COPY ./Analysis /Analysis
COPY ./Raw/Amikacin /Raw/Amikacin