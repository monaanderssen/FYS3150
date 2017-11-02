TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    lib.cpp
INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -larmadillo -llapack -lblas

HEADERS += \
    lib.h

