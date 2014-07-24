LIBRARY := BMinusLCutflow
OBJDIR  := $(PWD)/obj
SRCDIR  := $(PWD)/src
INCDIR  := $(PWD)/include
EXECUTABLE := $(PWD)/BMinusLCutflow

LOCAL_LIBS_PATH += $(SH_LIB_PATH)/libProgressBar.so
LOCAL_LIBS_PATH += $(SH_LIB_PATH)/libTruthNtupleLooper.so
LOCAL_LIBS_PATH += $(SH_LIB_PATH)/libTruthRecordHelpers.so
LOCAL_LIBS_PATH += $(SH_LIB_PATH)/libHistogramHandlers.so
LOCAL_LIBS_PATH += $(SH_LIB_PATH)/libLHAPDF.so

INCLUDES += -I../
INCLUDES += -I./
INCLUDES += -I../../
INCLUDES += -I../TruthNtupleLooper/include/
INCLUDES += -I../LHAPDF/include/

include $(PWD)/../Makefile.common
