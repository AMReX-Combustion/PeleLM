DEFINES += -DBL_USE_VELOCITY -DBL_NOLINEVALUES
CEXE_headers += DataServices.H AmrData.H XYPlotDataList.H AmrvisConstants.H
CEXE_sources += DataServices.cpp AmrData.cpp
FEXE_sources += FABUTIL_$(DIM)D.F
