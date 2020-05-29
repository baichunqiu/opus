#include $(call all-subdir-makefiles)

LOCAL_PATH_OLD := $(LOCAL_PATH)
LOCAL_PATH := $(call my-dir)
#include $(CLEAR_VARS)

# jitterbuffer
include $(LOCAL_PATH)/jitterbuffer/Android.mk

#celt codec
include $(LOCAL_PATH)/opus-1.2.1/Android.mk

#sbc codec
include $(LOCAL_PATH)/sbc/Android.mk

#beanformer
#include $(LOCAL_PATH)/beamformer/Android.mk

#resample
#include $(LOCAL_PATH)/resample/Android.mk

#authenticate
include $(LOCAL_PATH)/authenticate/Android.mk