FROM ubuntu:16.04

#HOW TO RUN
#1. Open Terminal
#2. Open XQuartz: open -a XQuartz (download from: https://www.xquartz.org)
#3. Set IP: IP=$(ifconfig en0 | grep inet | awk '$1=="inet" {print $2}') if on ethernet
#3. or IP=$(ifconfig en1 | grep inet | awk '$1=="inet" {print $2}') if on wifi
#4. Add IP: xhost + $IP
#5. Build docker image: docker build -t mongoose .
#6. Run docker image:
# docker run -ti --rm=true -e DISPLAY=$IP:0 -v /tmp/.X11-unix:/tmp/.X11-unix mongoose
#7. Use MONGOOSE
#8. After saving results, copy any files from the container to your local machine with the following command
# on your local machine's terminal:
# docker cp <container_id>:/path/to/file/in/container /path/to/local/machine
# i.e.: docker cp 665bd8e6cbad:/MongooseGUI3/docker_results.txt /Users/hostName/Documents/Leonid/MongooseGUI3
#9. Exit Docker by the command: exit



# Update apt-get and Install git
RUN apt-get update && \
    apt-get install -y apt-utils git build-essential vim curl

# Installs python dependencies
RUN apt-get install -y python3

# Installs python dependencies
RUN apt-get install -y python3-dev python3-pip python3-pyqt4 qt4-qmake libqt4-dev

# Installs GUI framework and autoconf needed for bootstrap command
RUN apt-get install -y x11-apps autoconf libtool-bin libgmp3-dev

# Pulls from MONGOOSE repo
RUN git clone https://github.com/WGS-TB/MongooseGUI3 && \
    cd MongooseGUI3 && \
    pip3 install --upgrade pip && \
    pip3 install setuptools cython python-libsbml xlrd

#Installs SIP
RUN mkdir -p /opt/sip && \
    cd /opt/sip && \
    curl -L -o sip.tar.gz http://sourceforge.net/projects/pyqt/files/sip/sip-4.17/sip-4.17.tar.gz && \
    tar -xf sip.tar.gz && \
    cd /opt/sip/sip-* && \
    python3 configure.py && \
    make && \
    make install && \
    cd /opt && \
    rm -rf /opt/sip

#installs pyqt, must install using the following
RUN mkdir -p /opt/pyqt && \
    cd /opt/pyqt && pwd && ls && \
    curl -L -o pyqt4.tar.gz http://sourceforge.net/projects/pyqt/files/PyQt4/PyQt-4.11.3/PyQt-x11-gpl-4.11.3.tar.gz && \
    tar -xf pyqt4.tar.gz && \
    cd /opt/pyqt && pwd && ls && \
    cd /opt/pyqt/PyQt-x11-gpl-4.11.3/ && \
    pwd && ls && \
    python3 configure.py -c --confirm-license --no-designer-plugin -e QtCore -e QtGui -e QtWidgets && \
    make && \
    make install && \
    cd /opt && \
    rm -rf /opt/pyqt

#installs libtool
RUN mkdir -p /opt/libtool && \
    cd /opt/libtool && \
    curl -L -o libtool-2.4.6.tar.gz http://mirror.jre655.com/GNU/libtool/libtool-2.4.6.tar.gz && \
    tar -xf libtool-2.4.6.tar.gz && \
    cd libtool-2.4.6 && \
    ./configure  && \
    make && \
    make install && \
    cd /opt

# installs gmp
RUN mkdir -p /opt/gmp && \
    cd /opt/gmp && \
    curl -L -o gmp-6.0.0a.tar.bz2 https://gmplib.org/download/gmp/gmp-6.0.0a.tar.bz2 && \
    tar -xf gmp-6.0.0a.tar.bz2 &&\
    cd gmp-6.0.0 && \
    ./configure && \
    make && \
    make check && \
    make install && \
    cd /opt

# Installs QSopt-ex
RUN git clone https://github.com/jonls/qsopt-ex.git &&\
    cd qsopt-ex/ && \
    ./bootstrap && \
    mkdir build && \
    cd build && \
    ../configure && \
    make && \
    make check && \
    make install && \
    cd ../../

# Installs Python interface to QSopt-ex
RUN git clone https://github.com/jonls/python-qsoptex.git && \
    cd python-qsoptex/ && \
    python3 setup.py install && \
    export LD_LIBRARY_PATH="/usr/local/lib/" && \
    python3 test_qsoptex.py && \
    cd ../MongooseGUI3

# creates the necessary links for qsopt_ex
RUN ldconfig
