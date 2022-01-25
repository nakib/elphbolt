FROM debian:buster

# User name
ENV USERNAME=elf

# Select and set time zone
# Need to do this, or else docker build prompts for TZ and then stalls.
ENV TZ=Europe/Madrid
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# This section is needed only if OpenMPI is used in the OpenCoarrays install.
# For an MPICH build (which I am going to do), this is not necessary.
#
# Below is a fix for the OpenMPI bug where a bunch of weird stuff is outputted
# for number of processes (or, coarray images) > 1. The numerical results
# are unaffected.
#
# I am not sure if this fix affects performance.
#ENV OMPI_MCA_btl_vader_single_copy_mechanism=none

# Update software repo
RUN apt-get update

# Upgrade OS
RUN apt-get upgrade -yq && apt-get dist-upgrade -yq

# Install sudo
RUN apt-get install -yq sudo

# Create, give sudo powers, and switch to user "elf"
# This section is from Milan Curcic (https://github.com/modern-fortran/modern-fortran-docker)
RUN echo "$USERNAME ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers && \
    useradd --no-log-init --home-dir /home/$USERNAME --create-home --shell /bin/bash $USERNAME && \
    adduser $USERNAME sudo

# Set user name and work directory
USER $USERNAME
WORKDIR /home/$USERNAME

# Install some other packages from Ubuntu repo
RUN sudo apt-get install -yq git nano make cmake gfortran mpich  \
    liblapack-dev libsymspg-dev && \
    sudo apt-get clean -q

# Install OpenCoarrays 2.8.0
RUN rm -rf OpenCoarrays # Remove previous build, if it exists.
RUN git clone https://github.com/sourceryinstitute/OpenCoarrays && \
    mkdir OpenCoarrays/opencoarrays-install  && \
    cd OpenCoarrays/opencoarrays-install && \
    git checkout tags/2.8.0 && \
    FC="$(command -v gfortran)" CC="$(command -v gcc)" cmake .. && \
    sudo make install && \
    caf --version && \
    cafrun --version

# Install elphbolt (develop-latest)
RUN rm -rf elphbolt # Remove previous builds if it exists.
RUN git clone https://github.com/nakib/elphbolt; cd elphbolt && \
    cp Makefiles/thinkpad_gcc.make src/ && \
    cp Makefiles/Makefile src/ &&\
    cd src/ && \
    make clean; make
