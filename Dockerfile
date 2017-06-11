FROM ubuntu:16.04

# Update apt-get and Install git
RUN apt-get update && \
    apt-get install -y git build-essential

RUN apt-get install -y python3 python3-dev python3-pip python3-pyqt4 libqt4-dev python3-numpy

RUN apt-get install -y x11-apps autoconf

RUN git clone https://github.com/WGS-TB/MongooseGUI3 && \
    cd MongooseGUI3 && \
    pip3 install --upgrade pip && \
    pip3 install setuptools cython python-libsbml xlrd sip

# Install libtool-bin and libgmp3-dev
RUN apt-get install -y libtool-bin libgmp3-dev

RUN git clone https://github.com/jonls/qsopt-ex.git &&\
    cd qsopt-ex/ && \
    ./bootstrap && \
    mkdir build && \
    cd build && \
    ../configure && \
    make && \
    make install && \
    cd ../../

RUN git clone https://github.com/jonls/python-qsoptex.git && \
    cd python-qsoptex/ && \
    python3 setup.py install && \
    export LD_LIBRARY_PATH="/usr/local/lib/" && \
    python3 test_qsoptex.py && \
    cd ../MongooseGUI3 && \
    ls -a && pwd && which python3

CMD [ "python3","/MongooseGUI3/MONGOOSEgui.py"]
#RUN python3 /MongooseGUI3/MONGOOSEgui.py
