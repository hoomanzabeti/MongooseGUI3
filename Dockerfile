FROM ubuntu:16.04

#HOW TO RUN
#1. Open Terminal
#2. Open XQuartz: open -a XQuartz
#3. Set IP: IP=$(ifconfig en0 | grep inet | awk '$1=="inet" {print $2}') if on ethernet
#3. or IP=$(ifconfig en1 | grep inet | awk '$1=="inet" {print $2}') if on wifi
#4. Add IP: xhost + $IP
#5. Build docker image: docker build -t mongoose .
#6 Run docker image:
# docker run -ti --rm=true -e DISPLAY=$IP:0 -v /tmp/.X11-unix:/tmp/.X11-unix mongoose


# adds user
RUN adduser --quiet --disabled-password qtuser

# Update apt-get and Install git
RUN apt-get update && \
    apt-get install -y apt-utils git build-essential vim curl

# Installs python dependencies
RUN apt-get install -y python3

# Installs python dependencies
RUN apt-get install -y python3-dev python3-pip python3-pyqt4 qt4-qmake libqt4-dev

# Installs GUI framework and autoconf needed for bootstrap command
RUN apt-get install -y x11-apps autoconf libgmp3-dev libtool-bin libgmp3-dev libxkbcommon-dev libxcb-xkb-dev libxslt1-dev libgstreamer-plugins-base0.10-dev

# Pulls from MONGOOSE repo
RUN git clone https://github.com/WGS-TB/MongooseGUI3 && \
    cd MongooseGUI3 && \
    pip3 install --upgrade pip && \
    pip3 install setuptools cython python-libsbml xlrd

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


# Installs QSopt-ex
RUN git clone https://github.com/jonls/qsopt-ex.git &&\
    cd qsopt-ex/ && \
    ./bootstrap && \
    mkdir build && \
    cd build && \
    ../configure && \
    make && \
    make install && \
    cd ../../

# Installs Python interface to QSopt-ex
RUN git clone https://github.com/jonls/python-qsoptex.git && \
    cd python-qsoptex/ && \
    python3 setup.py install && \
    export LD_LIBRARY_PATH="/usr/local/lib/" && \
    python3 test_qsoptex.py && \
    cd ../MongooseGUI3

# Runs GUI
CMD [ "python3","/MongooseGUI3/MONGOOSEgui.py"]
#RUN python3 /MongooseGUI3/MONGOOSEgui.py
