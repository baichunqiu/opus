LOCAL_PATH_OLD := $(LOCAL_PATH)
LOCAL_PATH := $(call my-dir)

include $(CLEAR_VARS)

LOCAL_MODULE := jni-beamformer
LOCAL_LDLIBS := -llog
LOCAL_C_INCLUDES := $(LOCAL_PATH)/include
LOCAL_SRC_FILES :=src/beamformer.c src/icm_fft.c src/ring_buffer.c src/BeamformerJni.cpp

include $(BUILD_SHARED_LIBRARY)
LOCAL_PATH := $(LOCAL_PATH_OLD)