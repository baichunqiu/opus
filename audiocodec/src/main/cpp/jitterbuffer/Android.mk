LOCAL_PATH_OLD := $(LOCAL_PATH)
LOCAL_PATH := $(call my-dir)

include $(CLEAR_VARS)

LOCAL_MODULE    :=jni-jitterbuffer
LOCAL_LDLIBS:=-L$(SYSROOT)/usr/lib -llog
LOCAL_C_INCLUDES    := $(LOCAL_PATH)/include
LOCAL_SRC_FILES :=src/CircBuffer.cpp src/JitterBuffer.cpp src/DelayBuf.cpp src/Wsola.cpp

include $(BUILD_SHARED_LIBRARY)
LOCAL_PATH := $(LOCAL_PATH_OLD)
  