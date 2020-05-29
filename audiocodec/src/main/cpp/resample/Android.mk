LOCAL_PATH_OLD := $(LOCAL_PATH)
LOCAL_PATH := $(call my-dir)

include $(CLEAR_VARS)

LOCAL_MODULE    :=jni-resample
LOCAL_LDLIBS:=-L$(SYSROOT)/usr/lib -llog
LOCAL_C_INCLUDES    := $(LOCAL_PATH)/include
LOCAL_SRC_FILES :=src/CodecResample.cpp src/S_resample.cpp src/Mixer_RS_KW.cpp src/Mixer_RS_Linear.cpp

include $(BUILD_SHARED_LIBRARY)
LOCAL_PATH := $(LOCAL_PATH_OLD)
  