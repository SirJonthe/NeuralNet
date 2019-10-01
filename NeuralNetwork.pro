TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        MiniLib/MGL/mglCollision.cpp \
        MiniLib/MGL/mglColor.cpp \
        MiniLib/MGL/mglImage.cpp \
        MiniLib/MGL/mglModel.cpp \
        MiniLib/MGL/mglParticleDynamics2D.cpp \
        MiniLib/MGL/mglPlane.cpp \
        MiniLib/MGL/mglTexture.cpp \
        MiniLib/MGL/mglTransform.cpp \
        MiniLib/MML/mmlNoise.cpp \
        MiniLib/MML/mmlRandom.cpp \
        MiniLib/MPL/mplWide.cpp \
        MiniLib/MTL/mtlMathParser.cpp \
        MiniLib/MTL/mtlParser.cpp \
        MiniLib/MTL/mtlPath.cpp \
        MiniLib/MTL/mtlString.cpp \
        MiniLib/MTL/mtlType.cpp \
        layer.cpp \
        main.cpp \
        matrix.cpp \
        neuron.cpp \
    neural_network.cpp

HEADERS += \
    MiniLib/MGL/mglCollision.h \
    MiniLib/MGL/mglColor.h \
    MiniLib/MGL/mglImage.h \
    MiniLib/MGL/mglModel.h \
    MiniLib/MGL/mglParticleDynamics2D.h \
    MiniLib/MGL/mglPixel.h \
    MiniLib/MGL/mglPlane.h \
    MiniLib/MGL/mglTexture.h \
    MiniLib/MGL/mglTransform.h \
    MiniLib/MML/mmlFixed.h \
    MiniLib/MML/mmlInt.h \
    MiniLib/MML/mmlMath.h \
    MiniLib/MML/mmlMatrix.h \
    MiniLib/MML/mmlNoise.h \
    MiniLib/MML/mmlQuaternion.h \
    MiniLib/MML/mmlRandom.h \
    MiniLib/MML/mmlVector.h \
    MiniLib/MPL/mplAlloc.h \
    MiniLib/MPL/mplCommon.h \
    MiniLib/MPL/mplMath.h \
    MiniLib/MPL/mplWide.h \
    MiniLib/MTL/mtlArray.h \
    MiniLib/MTL/mtlAsset.h \
    MiniLib/MTL/mtlBinaryTree.h \
    MiniLib/MTL/mtlBits.h \
    MiniLib/MTL/mtlDuplex.h \
    MiniLib/MTL/mtlHashTable.h \
    MiniLib/MTL/mtlList.h \
    MiniLib/MTL/mtlMathParser.h \
    MiniLib/MTL/mtlMemory.h \
    MiniLib/MTL/mtlParser.h \
    MiniLib/MTL/mtlPath.h \
    MiniLib/MTL/mtlPointer.h \
    MiniLib/MTL/mtlString.h \
    MiniLib/MTL/mtlStringMap.h \
    MiniLib/MTL/mtlType.h \
    layer.h \
    matrix.h \
    neuron.h \
    neural_network.h
