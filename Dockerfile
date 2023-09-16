# Use an appropriate base image with a Fortran compiler (e.g., gfortran) and FPM installed.
FROM debian:buster

# User name
ENV USERNAME=elf

# Select and set time zone
ENV TZ=Europe/Madrid
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Update software repo
RUN apt-get update

# Upgrade OS
RUN apt-get upgrade -yq && apt-get dist-upgrade -yq

# Install sudo
RUN apt-get install -yq sudo

# Create, give sudo powers, and switch to user "elf"
RUN echo "$USERNAME ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers && \
    useradd --no-log-init --home-dir /home/$USERNAME --create-home --shell /bin/bash $USERNAME && \
    adduser $USERNAME sudo

# Set user name and work directory
USER $USERNAME
WORKDIR /home/$USERNAME

# Install some other packages from Debian repo
RUN sudo apt-get install -yq git nano make cmake gfortran mpich  \
    liblapack-dev libsymspg-dev && \
    sudo apt-get clean -q

# Install FPM
RUN sudo apt-get install -yq ruby && \
    sudo gem install fpm

# Install elphbolt using FPM
RUN git clone https://github.com/nakib/elphbolt && \
    cd elphbolt && \
    fpm install && \
    cd .. && rm -rf elphbolt

# Optionally, install other Fortran dependencies using FPM
# RUN fpm install <dependency-name>

# Copy your Fortran project files (FPM.toml, source code, etc.) into the container
# COPY . /home/$USERNAME/my_fortran_project

# Build and install your Fortran project using FPM
# RUN cd /home/$USERNAME/my_fortran_project && \
#     fpm build && \
#     fpm install
