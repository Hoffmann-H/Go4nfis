######################################################################
# Automatically generated by qmake (2.01a) Fri Jan 11 08:19:38 2019
######################################################################

TEMPLATE = app
TARGET = FC-Analysis
CONFIG += console
CONFIG -= qt

INCLUDEPATH += $(ROOTSYS)/include \
                include \
                $(GO4SYS)/include \
                src

LIBS += $$system(root-config --libs) \
        $$system(root-config --glibs) \
        -L$(ROOTSYS)/lib -lASImage

QMAKE_CXXFLAGS += -O2 \
                $(shell root-config --cflags)
QMAKE_CXXFLAGS_WARN_ON += -Wno-unused-parameter
QMAKE_CFLAGS += -O2 \
                $(shell root-config --cflags)
QMAKE_LFLAGS += -O  \
                $(shell root-config --ldflags)

unix {
    UI_DIR = .ui
    MOC_DIR = .moc
    OBJECTS_DIR = .obj
}

# Input
#HEADERS += include/FC.h \
#           include/PuFC.h \
#           include/UFC.h \
#           include/Sim.h \
#           include/AnaSim.h \
#           include/Plot.h \
#           include/ToF.h

SOURCES += \
#    src/FC.cpp \
#    src/PuFC.cpp \
#    src/UFC.cpp \
#    src/Sim.cpp \
#    src/AnaSim.cpp \
#    src/Plot.cpp \
#    src/ToF.cpp \
    main.cpp \
    root/AnaSim.C \
    root/Carlson.C \
    root/Correction.C \
    root/CrossSection.C \
    root/Deposit.C \
    root/DrawPics.C \
    root/drawTL.C \
    root/FC.C \
    root/MCNPtoROOT.C \
    root/nELBEsim.C \
    root/NeutronField.C \
    root/NumberOfAtoms.C \
    root/PeakWidth.C \
    root/QDCmin.C \
    root/Runs.C \
    root/SaveToFile.C \
    root/ShadowCone.C \
    root/Stability.C \
    root/Target.C \
    root/ToF.C \
    root/VglSim.C \
    root/Transmission.C

!exists ($(ROOTSYS)/include/rootcintrule.pri):message ("The rootcintrules.pri was not found")
exists ($(ROOTSYS)/include/rootcintrule.pri) {
    include ($(ROOTSYS)/include/rootcintrule.pri)
    include ($(ROOTSYS)/include/rootcint.pri)
    CREATE_ROOT_DICT_FOR_CLASSES = ${HEADERS} LinkDef.h
    CREATE_ROOT_DICT_FOR_CLASSES -= ui_mainwindow

OTHER_FILES += \
        ../../../StyleSheets/StyleSheet.C \
}
