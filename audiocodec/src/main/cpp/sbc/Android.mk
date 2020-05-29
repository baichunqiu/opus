LOCAL_PATH_OLD := $(LOCAL_PATH)
LOCAL_PATH := $(call my-dir)

include $(CLEAR_VARS)

LOCAL_MODULE    := jni-sbc
LOCAL_LDLIBS :=-llog
LOCAL_C_INCLUDES    := $(LOCAL_PATH)/include
LOCAL_SRC_FILES := src/lib_sbc.c  src/sbc.c  src/sbc_primitives.c

include $(BUILD_SHARED_LIBRARY)
LOCAL_PATH := $(LOCAL_PATH_OLD)