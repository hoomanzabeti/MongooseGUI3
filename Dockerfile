FROM ubuntu:16.04

#HOW TO RUN
#To run this Dockerfile the following is necessary ( a Mac OS is assumed ): do not include the beginning '$', this is to simply indicate a terminal command.
#
#Open Terminal
#
#Open XQuartz: (download from: https://www.xquartz.org)
#
#    $  open -a XQuartz
#3.1: Set IP: if on ethernet
#
#    $  IP=$(ifconfig en0 | grep inet | awk '$1=="inet" {print $2}')
#3.2: or if on wifi
#
#    $  IP=$(ifconfig en1 | grep inet | awk '$1=="inet" {print $2}')
#Add IP:
#
#    $   xhost + $IP
#Next we would like to build the docker image. This step can be done one of two ways:
#5.1:
#a. Pull our github repository (https://github.com/WGS-TB/MongooseGUI3) using the command:
#
#    $  git clone https://github.com/WGS-TB/MongooseGUI3.git
#b. cd into the cloned directory
#c. Build docker image:
#
#    $  docker build -t mongoose .
#OR
#5.2:
#a.
#
#    $  docker pull ctlevn/mongoose
#At the end of either methods, run :
#
#    $  docker images
#to confirm that the docker image has been built successfully. You will see something similar to the following:
#REPOSITORY TAG IMAGE ID CREATED SIZE
#mongoose latest f3e77cc7b7b0 2 days ago 1.21 GB
#ctlevn/mongoose latest f3e77cc7b7b0 2 days ago 1.21 GB
#
#Run docker image:
#If you did step 5.1 above, perform the following command:
#
#    $  docker run -ti --rm=true -e DISPLAY=$IP:0 -v /tmp/.X11-unix:/tmp/.X11-unix mongoose
#If you did step 5.2 above, perform the following command:
#
#    $  docker run -ti --rm=true -e DISPLAY=$IP:0 -v /tmp/.X11-unix:/tmp/.X11-unix ctlevn/mongoose
#The terminal where you ran the command is now an interactive terminal for the Docker container.
#Terminal commands such as ls, cd, and etc work. To go to the MongooseGUI3 directory run the command:
#
#    $  cd /MongooseGUI3
#Now we need to create a shelve file to analyze. The XML files can be downloaded from here:
#http://cb.csail.mit.edu/cb/mongoose/models.html (middle column)
#
#Once the file is downloaded, we can copy the file from our local machine to the container:
#First open up another terminal to navigate the local machine. There should now be two terminals open,
#one for Docker and one for the local machine.
#Now in the local machines terminal run the following:
#
#    $  docker ps -a
#This will list the running processes. Observe the container ID.
#In order to copy the downloaded file from our local machine to the container, run the following:
#
#    $  docker cp /path/to/local/machine <container_id>:/MongooseGUI3
#i.e.:
#
#    $  docker cp /Users/hostName/Documents/Leonid/MongooseGUI3/AG1Model 665bd8e6cbad:/MongooseGUI3
#where: /path/to/local/machine is the path of the file on your local machine.
#<container_id> is the corresponding container id
#/MongooseGUI3 is the proper destination path in the container.
#
#To create a shelve file run the following:
#
#   $  python3 -i ModelParsing.py
#Then, within Python:
#
#model = parseSBML('<replace_with_file_name>.xml')
#model.adjustCompartments('xt', startPos = -2)
#model.biomassCoefficients[760] = 1
#s=shelve.open('ParsedModel')
#s['AB1'] = model
#s.close()
#exit()
#
#After this you will have the file called 'DockerModel' that will contain the parsed model inside it
#Then to run the GUI, perform the command:
#$ python3 MongooseGUI3
#The GUI should now start up.
#
#After saving results, copy any files from the container to your local machine with the following command
#on your local machine's terminal:
#
#   $  docker cp <container_id>:/path/to/file/in/container /path/to/local/machine
#i.e:
#
#   $  docker cp 665bd8e6cbad:/MongooseGUI3/docker_results.txt /Users/hostName/Documents/Leonid/MongooseGUI3
#Exit Docker by the command: exit



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
