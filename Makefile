#############################################################################
# Makefile for building: FluidSim
# Generated by qmake (3.0) (Qt 5.7.1)
# Project:  FluidSim.pro
# Template: app
# Command: /usr/lib/x86_64-linux-gnu/qt5/bin/qmake -spec linux-g++-64 -o Makefile FluidSim.pro
#############################################################################

MAKEFILE      = Makefile

####### Compiler, tools and options

CC            = gcc
CXX           = g++
DEFINES       = -DQT_NO_DEBUG -DQT_GUI_LIB -DQT_CORE_LIB
CFLAGS        = -m64 -pipe -O2 -Wall -W -D_REENTRANT -fPIC $(DEFINES)
CXXFLAGS      = -m64 -pipe -std=c++14 -O2 -std=gnu++1y -Wall -W -D_REENTRANT -fPIC $(DEFINES)
INCPATH       = -I. -Iglm -isystem /usr/include/x86_64-linux-gnu/qt5 -isystem /usr/include/x86_64-linux-gnu/qt5/QtGui -isystem /usr/include/x86_64-linux-gnu/qt5/QtCore -I. -I/usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++-64
QMAKE         = /usr/lib/x86_64-linux-gnu/qt5/bin/qmake
DEL_FILE      = rm -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p
COPY          = cp -f
COPY_FILE     = cp -f
COPY_DIR      = cp -f -R
INSTALL_FILE  = install -m 644 -p
INSTALL_PROGRAM = install -m 755 -p
INSTALL_DIR   = cp -f -R
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
TAR           = tar -cf
COMPRESS      = gzip -9f
DISTNAME      = FluidSim1.0.0
DISTDIR = /gpfs/main/home/dtian2/course/cs2240/final/.tmp/FluidSim1.0.0
LINK          = g++
LFLAGS        = -m64 -Wl,-O1
LIBS          = $(SUBLIBS) -L/usr/X11R6/lib64 -lQt5Gui -lQt5Core -lGL -lpthread 
AR            = ar cqs
RANLIB        = 
SED           = sed
STRIP         = strip

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = main.cpp \
		Renderer.cpp \
		FluidSolver.cpp \
		LevelSetSolver.cpp 
OBJECTS       = main.o \
		Renderer.o \
		FluidSolver.o \
		LevelSetSolver.o
DIST          = level_sets.txt \
		level_sets_2.txt \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/spec_pre.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/unix.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/linux.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/sanitize.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/gcc-base.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/gcc-base-unix.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/g++-base.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/g++-unix.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/qconfig.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_bootstrap_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_clucene_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_concurrent.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_concurrent_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_core.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_core_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_dbus.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_dbus_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_designer.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_designer_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_designercomponents_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_eglfs_device_lib_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_eglfs_kms_support_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_gui.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_gui_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_help.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_help_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_network.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_network_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_opengl.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_opengl_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_openglextensions.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_openglextensions_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_platformsupport_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_printsupport.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_printsupport_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_qml.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_qmltest.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_quick.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_quickwidgets.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_sql.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_sql_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_testlib.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_testlib_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_uiplugin.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_uitools.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_uitools_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_webkit.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_webkitwidgets.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_widgets.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_widgets_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xcb_qpa_lib_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xml.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xml_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qt_functions.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qt_config.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++-64/qmake.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/spec_post.prf \
		.qmake.stash \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/exclusive_builds.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/toolchain.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/default_pre.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/resolve_config.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/default_post.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/warn_on.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qt.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/resources.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/moc.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/unix/opengl.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/unix/thread.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/file_copies.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/testcase_targets.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/exceptions.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/yacc.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/lex.prf \
		FluidSim.pro Renderer.h \
		FluidSolver.h \
		LevelSetSolver.h main.cpp \
		Renderer.cpp \
		FluidSolver.cpp \
		LevelSetSolver.cpp
QMAKE_TARGET  = FluidSim
DESTDIR       = 
TARGET        = FluidSim


first: all
####### Build rules

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

Makefile: FluidSim.pro /usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++-64/qmake.conf /usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/spec_pre.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/unix.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/linux.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/sanitize.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/gcc-base.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/gcc-base-unix.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/g++-base.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/g++-unix.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/qconfig.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_bootstrap_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_clucene_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_concurrent.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_concurrent_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_core.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_core_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_dbus.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_dbus_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_designer.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_designer_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_designercomponents_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_eglfs_device_lib_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_eglfs_kms_support_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_gui.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_gui_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_help.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_help_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_network.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_network_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_opengl.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_opengl_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_openglextensions.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_openglextensions_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_platformsupport_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_printsupport.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_printsupport_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_qml.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_qmltest.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_quick.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_quickwidgets.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_sql.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_sql_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_testlib.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_testlib_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_uiplugin.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_uitools.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_uitools_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_webkit.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_webkitwidgets.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_widgets.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_widgets_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xcb_qpa_lib_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xml.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xml_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qt_functions.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qt_config.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++-64/qmake.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/spec_post.prf \
		.qmake.stash \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/exclusive_builds.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/toolchain.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/default_pre.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/resolve_config.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/default_post.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/warn_on.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qt.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/resources.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/moc.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/unix/opengl.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/unix/thread.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/file_copies.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/testcase_targets.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/exceptions.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/yacc.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/lex.prf \
		FluidSim.pro \
		/usr/lib/x86_64-linux-gnu/libQt5Gui.prl \
		/usr/lib/x86_64-linux-gnu/libQt5Core.prl
	$(QMAKE) -spec linux-g++-64 -o Makefile FluidSim.pro
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/spec_pre.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/unix.conf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/linux.conf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/sanitize.conf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/gcc-base.conf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/gcc-base-unix.conf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/g++-base.conf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/g++-unix.conf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/qconfig.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_bootstrap_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_clucene_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_concurrent.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_concurrent_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_core.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_core_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_dbus.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_dbus_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_designer.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_designer_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_designercomponents_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_eglfs_device_lib_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_eglfs_kms_support_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_gui.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_gui_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_help.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_help_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_network.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_network_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_opengl.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_opengl_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_openglextensions.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_openglextensions_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_platformsupport_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_printsupport.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_printsupport_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_qml.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_qmltest.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_quick.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_quickwidgets.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_sql.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_sql_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_testlib.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_testlib_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_uiplugin.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_uitools.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_uitools_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_webkit.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_webkitwidgets.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_widgets.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_widgets_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xcb_qpa_lib_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xml.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xml_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qt_functions.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qt_config.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++-64/qmake.conf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/spec_post.prf:
.qmake.stash:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/exclusive_builds.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/toolchain.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/default_pre.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/resolve_config.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/default_post.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/warn_on.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qt.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/resources.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/moc.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/unix/opengl.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/unix/thread.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/file_copies.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/testcase_targets.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/exceptions.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/yacc.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/lex.prf:
FluidSim.pro:
/usr/lib/x86_64-linux-gnu/libQt5Gui.prl:
/usr/lib/x86_64-linux-gnu/libQt5Core.prl:
qmake: FORCE
	@$(QMAKE) -spec linux-g++-64 -o Makefile FluidSim.pro

qmake_all: FORCE


all: Makefile $(TARGET)

dist: distdir FORCE
	(cd `dirname $(DISTDIR)` && $(TAR) $(DISTNAME).tar $(DISTNAME) && $(COMPRESS) $(DISTNAME).tar) && $(MOVE) `dirname $(DISTDIR)`/$(DISTNAME).tar.gz . && $(DEL_FILE) -r $(DISTDIR)

distdir: FORCE
	@test -d $(DISTDIR) || mkdir -p $(DISTDIR)
	$(COPY_FILE) --parents $(DIST) $(DISTDIR)/
	$(COPY_FILE) --parents Renderer.h FluidSolver.h LevelSetSolver.h $(DISTDIR)/
	$(COPY_FILE) --parents main.cpp Renderer.cpp FluidSolver.cpp LevelSetSolver.cpp $(DISTDIR)/


clean: compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


distclean: clean 
	-$(DEL_FILE) $(TARGET) 
	-$(DEL_FILE) .qmake.stash
	-$(DEL_FILE) Makefile


####### Sub-libraries

mocclean: compiler_moc_header_clean compiler_moc_source_clean

mocables: compiler_moc_header_make_all compiler_moc_source_make_all

check: first

benchmark: first

compiler_rcc_make_all:
compiler_rcc_clean:
compiler_moc_header_make_all:
compiler_moc_header_clean:
compiler_moc_source_make_all:
compiler_moc_source_clean:
compiler_yacc_decl_make_all:
compiler_yacc_decl_clean:
compiler_yacc_impl_make_all:
compiler_yacc_impl_clean:
compiler_lex_make_all:
compiler_lex_clean:
compiler_clean: 

####### Compile

main.o: main.cpp Renderer.h \
		glm/glm.hpp \
		glm/detail/_fixes.hpp \
		glm/fwd.hpp \
		glm/detail/type_int.hpp \
		glm/detail/setup.hpp \
		glm/detail/type_float.hpp \
		glm/detail/type_vec.hpp \
		glm/detail/precision.hpp \
		glm/detail/type_mat.hpp \
		glm/vec2.hpp \
		glm/detail/type_vec2.hpp \
		glm/detail/_swizzle.hpp \
		glm/detail/_swizzle_func.hpp \
		glm/detail/type_vec2.inl \
		glm/vec3.hpp \
		glm/detail/type_vec3.hpp \
		glm/detail/type_vec3.inl \
		glm/vec4.hpp \
		glm/detail/type_vec4.hpp \
		glm/detail/type_vec4.inl \
		glm/mat2x2.hpp \
		glm/detail/type_mat2x2.hpp \
		glm/detail/type_mat2x2.inl \
		glm/mat2x3.hpp \
		glm/detail/type_mat2x3.hpp \
		glm/detail/type_mat2x3.inl \
		glm/mat2x4.hpp \
		glm/detail/type_mat2x4.hpp \
		glm/detail/type_mat2x4.inl \
		glm/mat3x2.hpp \
		glm/detail/type_mat3x2.hpp \
		glm/detail/type_mat3x2.inl \
		glm/mat3x3.hpp \
		glm/detail/type_mat3x3.hpp \
		glm/detail/type_mat3x3.inl \
		glm/mat3x4.hpp \
		glm/detail/type_mat3x4.hpp \
		glm/detail/type_mat3x4.inl \
		glm/mat4x2.hpp \
		glm/detail/type_mat4x2.hpp \
		glm/detail/type_mat4x2.inl \
		glm/mat4x3.hpp \
		glm/detail/type_mat4x3.hpp \
		glm/detail/type_mat4x3.inl \
		glm/mat4x4.hpp \
		glm/detail/type_mat4x4.hpp \
		glm/detail/type_mat4x4.inl \
		glm/trigonometric.hpp \
		glm/detail/func_trigonometric.hpp \
		glm/detail/func_trigonometric.inl \
		glm/detail/_vectorize.hpp \
		glm/detail/type_vec1.hpp \
		glm/detail/type_vec1.inl \
		glm/exponential.hpp \
		glm/detail/func_exponential.hpp \
		glm/detail/func_exponential.inl \
		glm/detail/func_vector_relational.hpp \
		glm/detail/func_vector_relational.inl \
		glm/common.hpp \
		glm/detail/func_common.hpp \
		glm/detail/func_common.inl \
		glm/packing.hpp \
		glm/detail/func_packing.hpp \
		glm/detail/func_packing.inl \
		glm/detail/type_half.hpp \
		glm/detail/type_half.inl \
		glm/geometric.hpp \
		glm/detail/func_geometric.hpp \
		glm/detail/func_geometric.inl \
		glm/matrix.hpp \
		glm/detail/func_matrix.hpp \
		glm/detail/func_matrix.inl \
		glm/vector_relational.hpp \
		glm/integer.hpp \
		glm/detail/func_integer.hpp \
		glm/detail/func_integer.inl \
		FluidSolver.h \
		LevelSetSolver.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main.o main.cpp

Renderer.o: Renderer.cpp Renderer.h \
		glm/glm.hpp \
		glm/detail/_fixes.hpp \
		glm/fwd.hpp \
		glm/detail/type_int.hpp \
		glm/detail/setup.hpp \
		glm/detail/type_float.hpp \
		glm/detail/type_vec.hpp \
		glm/detail/precision.hpp \
		glm/detail/type_mat.hpp \
		glm/vec2.hpp \
		glm/detail/type_vec2.hpp \
		glm/detail/_swizzle.hpp \
		glm/detail/_swizzle_func.hpp \
		glm/detail/type_vec2.inl \
		glm/vec3.hpp \
		glm/detail/type_vec3.hpp \
		glm/detail/type_vec3.inl \
		glm/vec4.hpp \
		glm/detail/type_vec4.hpp \
		glm/detail/type_vec4.inl \
		glm/mat2x2.hpp \
		glm/detail/type_mat2x2.hpp \
		glm/detail/type_mat2x2.inl \
		glm/mat2x3.hpp \
		glm/detail/type_mat2x3.hpp \
		glm/detail/type_mat2x3.inl \
		glm/mat2x4.hpp \
		glm/detail/type_mat2x4.hpp \
		glm/detail/type_mat2x4.inl \
		glm/mat3x2.hpp \
		glm/detail/type_mat3x2.hpp \
		glm/detail/type_mat3x2.inl \
		glm/mat3x3.hpp \
		glm/detail/type_mat3x3.hpp \
		glm/detail/type_mat3x3.inl \
		glm/mat3x4.hpp \
		glm/detail/type_mat3x4.hpp \
		glm/detail/type_mat3x4.inl \
		glm/mat4x2.hpp \
		glm/detail/type_mat4x2.hpp \
		glm/detail/type_mat4x2.inl \
		glm/mat4x3.hpp \
		glm/detail/type_mat4x3.hpp \
		glm/detail/type_mat4x3.inl \
		glm/mat4x4.hpp \
		glm/detail/type_mat4x4.hpp \
		glm/detail/type_mat4x4.inl \
		glm/trigonometric.hpp \
		glm/detail/func_trigonometric.hpp \
		glm/detail/func_trigonometric.inl \
		glm/detail/_vectorize.hpp \
		glm/detail/type_vec1.hpp \
		glm/detail/type_vec1.inl \
		glm/exponential.hpp \
		glm/detail/func_exponential.hpp \
		glm/detail/func_exponential.inl \
		glm/detail/func_vector_relational.hpp \
		glm/detail/func_vector_relational.inl \
		glm/common.hpp \
		glm/detail/func_common.hpp \
		glm/detail/func_common.inl \
		glm/packing.hpp \
		glm/detail/func_packing.hpp \
		glm/detail/func_packing.inl \
		glm/detail/type_half.hpp \
		glm/detail/type_half.inl \
		glm/geometric.hpp \
		glm/detail/func_geometric.hpp \
		glm/detail/func_geometric.inl \
		glm/matrix.hpp \
		glm/detail/func_matrix.hpp \
		glm/detail/func_matrix.inl \
		glm/vector_relational.hpp \
		glm/integer.hpp \
		glm/detail/func_integer.hpp \
		glm/detail/func_integer.inl \
		FluidSolver.h \
		LevelSetSolver.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o Renderer.o Renderer.cpp

FluidSolver.o: FluidSolver.cpp FluidSolver.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o FluidSolver.o FluidSolver.cpp

LevelSetSolver.o: LevelSetSolver.cpp LevelSetSolver.h \
		glm/glm.hpp \
		glm/detail/_fixes.hpp \
		glm/fwd.hpp \
		glm/detail/type_int.hpp \
		glm/detail/setup.hpp \
		glm/detail/type_float.hpp \
		glm/detail/type_vec.hpp \
		glm/detail/precision.hpp \
		glm/detail/type_mat.hpp \
		glm/vec2.hpp \
		glm/detail/type_vec2.hpp \
		glm/detail/_swizzle.hpp \
		glm/detail/_swizzle_func.hpp \
		glm/detail/type_vec2.inl \
		glm/vec3.hpp \
		glm/detail/type_vec3.hpp \
		glm/detail/type_vec3.inl \
		glm/vec4.hpp \
		glm/detail/type_vec4.hpp \
		glm/detail/type_vec4.inl \
		glm/mat2x2.hpp \
		glm/detail/type_mat2x2.hpp \
		glm/detail/type_mat2x2.inl \
		glm/mat2x3.hpp \
		glm/detail/type_mat2x3.hpp \
		glm/detail/type_mat2x3.inl \
		glm/mat2x4.hpp \
		glm/detail/type_mat2x4.hpp \
		glm/detail/type_mat2x4.inl \
		glm/mat3x2.hpp \
		glm/detail/type_mat3x2.hpp \
		glm/detail/type_mat3x2.inl \
		glm/mat3x3.hpp \
		glm/detail/type_mat3x3.hpp \
		glm/detail/type_mat3x3.inl \
		glm/mat3x4.hpp \
		glm/detail/type_mat3x4.hpp \
		glm/detail/type_mat3x4.inl \
		glm/mat4x2.hpp \
		glm/detail/type_mat4x2.hpp \
		glm/detail/type_mat4x2.inl \
		glm/mat4x3.hpp \
		glm/detail/type_mat4x3.hpp \
		glm/detail/type_mat4x3.inl \
		glm/mat4x4.hpp \
		glm/detail/type_mat4x4.hpp \
		glm/detail/type_mat4x4.inl \
		glm/trigonometric.hpp \
		glm/detail/func_trigonometric.hpp \
		glm/detail/func_trigonometric.inl \
		glm/detail/_vectorize.hpp \
		glm/detail/type_vec1.hpp \
		glm/detail/type_vec1.inl \
		glm/exponential.hpp \
		glm/detail/func_exponential.hpp \
		glm/detail/func_exponential.inl \
		glm/detail/func_vector_relational.hpp \
		glm/detail/func_vector_relational.inl \
		glm/common.hpp \
		glm/detail/func_common.hpp \
		glm/detail/func_common.inl \
		glm/packing.hpp \
		glm/detail/func_packing.hpp \
		glm/detail/func_packing.inl \
		glm/detail/type_half.hpp \
		glm/detail/type_half.inl \
		glm/geometric.hpp \
		glm/detail/func_geometric.hpp \
		glm/detail/func_geometric.inl \
		glm/matrix.hpp \
		glm/detail/func_matrix.hpp \
		glm/detail/func_matrix.inl \
		glm/vector_relational.hpp \
		glm/integer.hpp \
		glm/detail/func_integer.hpp \
		glm/detail/func_integer.inl
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o LevelSetSolver.o LevelSetSolver.cpp

####### Install

install:  FORCE

uninstall:  FORCE

FORCE:

