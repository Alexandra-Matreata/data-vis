#-------------------------------------------------
#
# Project created by QtCreator 2017-02-28T01:27:16
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = SmokeApplet
TEMPLATE = app


SOURCES += main.cc\
        mainwindow.cc \
    smokewidget.cc \
    simulation.cc \
    legendwidget.cc

HEADERS  += mainwindow.h \
    smokewidget.h \
    simulation.h \
    legendwidget.h

FORMS    += mainwindow.ui

LIBS += -L"$$_PRO_FILE_PWD_/fftw-2.1.5/lib" -lrfftw -lfftw -lGL -lglut -lGLU -lm
