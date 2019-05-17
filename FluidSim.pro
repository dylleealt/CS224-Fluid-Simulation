TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
#CONFIG -= qt

SOURCES += main.cpp \
    Particles.cpp \
    Renderer.cpp \
    FluidSolver.cpp \
    LevelSetSolver.cpp

HEADERS += \
    Particles.h \
    Renderer.h \
    FluidSolver.h \
    LevelSetSolver.h

DISTFILES += \
    level_sets.txt \
    level_sets_2.txt

INCLUDEPATH += glm
DEPENDPATH += glm
