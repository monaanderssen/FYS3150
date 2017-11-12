TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    lib.cpp \
    lib.cpp
INCLUDEPATH += /usr/local/include /usr/l
LIBS += -L/usr/local/lib -larmadillo -llapack -lblas



HEADERS += \
    lib.h \
    lib.h


# FROM
# https://dragly.org/2012/03/14/developing-mpi-applications-in-qt-creator.html

# MPI Settings
QMAKE_CXX = /usr/lib64/openmpi/bin/mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = /usr/lib64/openmpi/bin/mpicc

QMAKE_CFLAGS += $$system(/usr/lib64/openmpi/bin/mpicc --showme:compile)
QMAKE_LFLAGS += $$system(/usr/lib64/openmpi/bin/mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(/usr/lib64/openmpi/bin/mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(/usr/lib64/openmpi/bin/mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK -O3
QMAKE_CXXFLAGS_DEBUG += $$system(/usr/lib64/openmpi/bin/mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK -O3
