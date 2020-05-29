LOCAL_PATH_OLD := $(LOCAL_PATH)
LOCAL_PATH := $(call my-dir)

include $(CLEAR_VARS)

LOCAL_MODULE    :=jni-authenticate
LOCAL_LDLIBS := -llog
LOCAL_C_INCLUDES := $(LOCAL_PATH)/include
LOCAL_SRC_FILES :=src/md5.cpp src/CheckAuthenticate.cpp src/authenticate.cpp

include $(BUILD_SHARED_LIBRARY)
LOCAL_PATH := $(LOCAL_PATH_OLD)
  