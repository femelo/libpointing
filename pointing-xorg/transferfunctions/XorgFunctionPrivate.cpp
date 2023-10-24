/* -*- mode: c++ -*-
 *
 * pointing-xorg/transferfunctions/XorgPredictableFunction.cpp --
 *
 * Derived software
 * Authors: Nicolas Roussel
 * Copyright © INRIA
 *
 */

#define Bool bool
// #define BOOL bool
#define TRUE true
#define FALSE false

#define _X_EXPORT
#define _X_INTERNAL

#if 1
#define ErrorF(...) /* */
#define DebugAccelF(...) /* */
#else
#include <stdio.h>
static void
ErrorF(const char * f, ...) {
  va_list args ;
  va_start(args, f) ;
  vfprintf(stdout, f, args) ;
  va_end(args) ;
}
#define DebugAccelF ErrorF
#endif

#define M_PI 3.1415926535897932384626433832795

#include <stdlib.h>
#include <stdint.h>
#include <string.h>

// -------------------------------------------------------------------------
// xorg-server-23.2.1/include/input.h (partial)
// -------------------------------------------------------------------------
/* lines 52-53 */
#include <X11/Xmd.h>
#include <X11/Xfuncproto.h>

/* lines 111-117 */
/*int constants for pointer acceleration schemes*/
#define PtrAccelNoOp            0
#define PtrAccelPredictable     1
#define PtrAccelLightweight     2
#define PtrAccelDefault         PtrAccelPredictable

#define MAX_VALUATORS 36

/* lines 151-152 */
typedef struct _DeviceIntRec *DeviceIntPtr;
typedef struct _ValuatorClassRec *ValuatorClassPtr;

/* line 159 */
typedef struct _DDXTouchPointInfo *DDXTouchPointInfoPtr;

/* line 162 */
typedef struct _ValuatorMask ValuatorMask;

/* lines 182-185 */
/* pointer acceleration handling */
typedef void (*PointerAccelSchemeProc) (DeviceIntPtr /*device */ ,
                                        ValuatorMask * /*valuators */ ,
                                        CARD32 /*evtime */ );

/* lines 189-192 */
struct _ValuatorAccelerationRec;
typedef Bool (*PointerAccelSchemeInitProc) (DeviceIntPtr /*dev */ ,
                                            struct _ValuatorAccelerationRec *
                                            /*protoScheme */ );

/* lines 222-225 */
typedef struct {
  int		num, den, threshold;
  //NR unsigned char id;
} PtrCtrl;

/* lines 362-363 */
typedef void (*PtrCtrlProcPtr) (DeviceIntPtr /*device */ ,
                                PtrCtrl * /*ctrl */ );

/* line 182 */
typedef void (*DeviceCallbackProc) (DeviceIntPtr /*pDev */ );

/* lines 742-783 */
extern _X_EXPORT ValuatorMask *valuator_mask_new(int num_valuators);
extern _X_EXPORT void valuator_mask_free(ValuatorMask **mask);
extern _X_EXPORT void valuator_mask_set_range(ValuatorMask *mask,
                                              int first_valuator,
                                              int num_valuators,
                                              const int *valuators);
extern _X_EXPORT void valuator_mask_set(ValuatorMask *mask, int valuator,
                                        int data);
extern _X_EXPORT void valuator_mask_set_double(ValuatorMask *mask, int valuator,
                                               double data);
extern _X_EXPORT void valuator_mask_zero(ValuatorMask *mask);
//NR extern _X_EXPORT int valuator_mask_size(const ValuatorMask *mask);
extern _X_EXPORT int valuator_mask_isset(const ValuatorMask *mask, int bit);
//NR extern _X_EXPORT void valuator_mask_unset(ValuatorMask *mask, int bit);
extern _X_EXPORT int valuator_mask_num_valuators(const ValuatorMask *mask);
#if 0 //NR
extern _X_EXPORT void valuator_mask_copy(ValuatorMask *dest,
                                         const ValuatorMask *src);
#endif //NR
extern _X_EXPORT int valuator_mask_get(const ValuatorMask *mask, int valnum);
extern _X_EXPORT double valuator_mask_get_double(const ValuatorMask *mask,
                                                 int valnum);
#if 0 //NR
extern _X_EXPORT Bool valuator_mask_fetch(const ValuatorMask *mask,
                                          int valnum, int *val);
extern _X_EXPORT Bool valuator_mask_fetch_double(const ValuatorMask *mask,
                                                 int valnum, double *val);
extern _X_EXPORT Bool valuator_mask_has_unaccelerated(const ValuatorMask *mask);
extern _X_EXPORT void valuator_mask_set_unaccelerated(ValuatorMask *mask,
                                                      int valuator,
                                                      double accel,
                                                      double unaccel);
extern _X_EXPORT void valuator_mask_set_absolute_unaccelerated(ValuatorMask *mask,
                                                               int valuator,
                                                               int absolute,
                                                               double unaccel);
extern _X_EXPORT double valuator_mask_get_accelerated(const ValuatorMask *mask,
                                                      int valuator);
extern _X_EXPORT double valuator_mask_get_unaccelerated(const ValuatorMask *mask,
                                                        int valuator);
extern _X_EXPORT Bool valuator_mask_fetch_unaccelerated(const ValuatorMask *mask,
                                                        int valuator,
                                                        double *accel,
                                                        double *unaccel);
extern _X_HIDDEN void valuator_mask_drop_unaccelerated(ValuatorMask *mask);
#endif //NR

// -------------------------------------------------------------------------
// xorg-server-23.2.1/include/inpututils.h (partial)
// -------------------------------------------------------------------------
/* lines 38-44 */
struct _ValuatorMask {
    int8_t last_bit;            /* highest bit set in mask */
    int8_t has_unaccelerated;
    uint8_t mask[(MAX_VALUATORS + 7) / 8];
    double valuators[MAX_VALUATORS];    /* valuator data */
    double unaccelerated[MAX_VALUATORS];    /* valuator data */
};


// -------------------------------------------------------------------------
// xorg-server-23.2.1/include/os.h (partial)
// -------------------------------------------------------------------------
/* lines 707-709 */
extern _X_EXPORT void
ErrorFSigSafe(const char *f, ...) { // changed definition for C++
    _X_ATTRIBUTE_PRINTF(1, 2);
}

/* lines 713-714 */
extern _X_EXPORT void
xorg_backtrace(void);


// -------------------------------------------------------------------------
// xorg-server-23.2.1/include/misc.h (partial)
// -------------------------------------------------------------------------
/* lines 135-136 */
#define _min(a, b) (((a) < (b)) ? (a) : (b)) // changed to _min to avoid clash with cmath
#define _max(a, b) (((a) > (b)) ? (a) : (b)) // changed to _max to avoid clash with cmath

/* 419-430 */
/* Don't use this directly, use BUG_WARN or BUG_WARN_MSG instead */
#define __BUG_WARN_MSG(cond, with_msg, ...)                                \
          do { if (cond) {                                                \
              ErrorFSigSafe("BUG: triggered 'if (" #cond ")'\n");          \
              ErrorFSigSafe("BUG: %s:%u in %s()\n",                        \
                           __FILE__, __LINE__, __func__);                 \
              if (with_msg) ErrorFSigSafe(__VA_ARGS__);                    \
              xorg_backtrace();                                           \
          } } while(0)

#define BUG_WARN_MSG(cond, ...)                                           \
          __BUG_WARN_MSG(cond, 1, __VA_ARGS__)

/* lines 440-441 */
#define BUG_RETURN_VAL(cond, val) \
        do { if (cond) { __BUG_WARN_MSG(cond, 0, NULL); return (val); } } while(0)


// -------------------------------------------------------------------------
// xorg-server-23.2.1/include/inputstr.h (partial)
// -------------------------------------------------------------------------

/* lines 64-66 */
#define BitIsOn(ptr, bit) (!!(((const BYTE *) (ptr))[(bit)>>3] & (1 << ((bit) & 7))))
#define SetBit(ptr, bit)  (((BYTE *) (ptr))[(bit)>>3] |= (1 << ((bit) & 7)))
#define ClearBit(ptr, bit) (((BYTE *)(ptr))[(bit)>>3] &= ~(1 << ((bit) & 7)))

/* lines 280-286 */
typedef struct _ValuatorAccelerationRec {
    int number;
    PointerAccelSchemeProc AccelSchemeProc;
    void *accelData;            /* at disposal of AccelScheme */
    PointerAccelSchemeInitProc AccelInitProc;
    DeviceCallbackProc AccelCleanupProc;
} ValuatorAccelerationRec, *ValuatorAccelerationPtr;

/* lines 288-303 */
typedef struct _ValuatorClassRec {
#if 0 //NR
    int sourceid;
    int numMotionEvents;
    int first_motion;
    int last_motion;
    void *motion;               /* motion history buffer. Different layout
                                   for MDs and SDs! */
    WindowPtr motionHintWindow;

    AxisInfoPtr axes;
    unsigned short numAxes;
    double *axisVal;            /* always absolute, but device-coord system */
#endif //NR
    ValuatorAccelerationRec accelScheme;
    //NR int h_scroll_axis;          /* horiz smooth-scrolling axis */
    //NR int v_scroll_axis;          /* vert smooth-scrolling axis */
} ValuatorClassRec;

/* line 414 */
typedef struct _PtrFeedbackClassRec *PtrFeedbackPtr;

/* lines 428-432 */
typedef struct _PtrFeedbackClassRec {
    //NR PtrCtrlProcPtr CtrlProc;
    PtrCtrl ctrl;
    //NR PtrFeedbackPtr next;
} PtrFeedbackClassRec;

/* lines 561-634 */
typedef struct _DeviceIntRec {
#if 0
    DeviceRec public;
    DeviceIntPtr next;
    Bool startup;               /* true if needs to be turned on at
                                   server initialization time */
    DeviceProc deviceProc;      /* proc(DevicePtr, DEVICE_xx). It is
                                   used to initialize, turn on, or
                                   turn off the device */
    Bool inited;                /* TRUE if INIT returns Success */
    Bool enabled;               /* TRUE if ON returns Success */
    Bool coreEvents;            /* TRUE if device also sends core */
    GrabInfoRec deviceGrab;     /* grab on the device */
    int type;                   /* MASTER_POINTER, MASTER_KEYBOARD, SLAVE */
    Atom xinput_type;
    char *name;
    int id;
    KeyClassPtr key;
#endif //NR
    ValuatorClassPtr valuator;
#if 0 //NR
    TouchClassPtr touch;
    GestureClassPtr gesture;
    ButtonClassPtr button;
    FocusClassPtr focus;
    ProximityClassPtr proximity;
    KbdFeedbackPtr kbdfeed;
#endif //NR
    PtrFeedbackPtr ptrfeed;
#if 0 //NR
    IntegerFeedbackPtr intfeed;
    StringFeedbackPtr stringfeed;
    BellFeedbackPtr bell;
    LedFeedbackPtr leds;
    struct _XkbInterest *xkb_interest;
    char *config_info;          /* used by the hotplug layer */
    ClassesPtr unused_classes;  /* for master devices */
    int saved_master_id;        /* for slaves while grabbed */
    PrivateRec *devPrivates;
    DeviceUnwrapProc unwrapProc;
    SpriteInfoPtr spriteInfo;
    DeviceIntPtr master;        /* master device */
    DeviceIntPtr lastSlave;     /* last slave device used */
#endif //NR
    /* last valuator values recorded, not posted to client;
     * for slave devices, valuators is in device coordinates, mapped to the
     * desktop
     * for master devices, valuators is in desktop coordinates.
     * see dix/getevents.c
     * remainder supports acceleration
     */
    struct {
        double valuators[MAX_VALUATORS];
        int numValuators;
        DeviceIntPtr slave;
        ValuatorMask *scroll;
        int num_touches;        /* size of the touches array */
        DDXTouchPointInfoPtr touches;
    } last;
#if 0 //NR
    /* Input device property handling. */
    struct {
        XIPropertyPtr properties;
        XIPropertyHandlerPtr handlers;  /* NULL-terminated */
    } properties;

    /* coordinate transformation matrix for relative movement. Matrix with
     * the translation component dropped */
    struct pixman_f_transform relative_transform;
    /* scale matrix for absolute devices, this is the combined matrix of
       [1/scale] . [transform] . [scale]. See DeviceSetTransform */
    struct pixman_f_transform scale_and_transform;

    /* XTest related master device id */
    int xtest_master_id;
    DeviceSendEventsProc sendEventsProc;
#endif //NR
    struct _SyncCounter *idle_counter;
} DeviceIntRec;

/* lines 809-810 */
extern _X_EXPORT void input_lock(void);
extern _X_EXPORT void input_unlock(void);

// -------------------------------------------------------------------------
// xorg-server-23.2.1/include/ptrveloc.h
// -------------------------------------------------------------------------

/*
 *
 * Copyright © 2006-2011 Simon Thum             simon dot thum at gmx dot de
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice (including the next
 * paragraph) shall be included in all copies or substantial portions of the
 * Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

//NR #ifndef POINTERVELOCITY_H
//NR #define POINTERVELOCITY_H

//NR #include <input.h>

/* constants for acceleration profiles */

#define AccelProfileNone -1
#define AccelProfileClassic  0
#define AccelProfileDeviceSpecific 1
#define AccelProfilePolynomial 2
#define AccelProfileSmoothLinear 3
#define AccelProfileSimple 4
#define AccelProfilePower 5
#define AccelProfileLinear 6
#define AccelProfileSmoothLimited 7
#define AccelProfileLAST AccelProfileSmoothLimited

/* fwd */
struct _DeviceVelocityRec;

/**
 * profile
 * returns actual acceleration depending on velocity, acceleration control,...
 */
typedef double (*PointerAccelerationProfileFunc)
 (DeviceIntPtr dev, struct _DeviceVelocityRec * vel,
  double velocity, double threshold, double accelCoeff);

/**
 * a motion history, with just enough information to
 * calc mean velocity and decide which motion was along
 * a more or less straight line
 */
typedef struct _MotionTracker {
    double dx, dy;              /* accumulated delta for each axis */
    int time;                   /* time of creation */
    int dir;                    /* initial direction bitfield */
} MotionTracker, *MotionTrackerPtr;

/**
 * Contains all data needed to implement mouse ballistics
 */
typedef struct _DeviceVelocityRec {
    MotionTrackerPtr tracker;
    int num_tracker;
    int cur_tracker;            /* current index */
    double velocity;            /* velocity as guessed by algorithm */
    double last_velocity;       /* previous velocity estimate */
    double last_dx;             /* last time-difference */
    double last_dy;             /* phase of last/current estimate */
    double corr_mul;            /* config: multiply this into velocity */
    double const_acceleration;  /* config: (recipr.) const deceleration */
    double min_acceleration;    /* config: minimum acceleration */
    short reset_time;           /* config: reset non-visible state after # ms */
    short use_softening;        /* config: use softening of mouse values */
    double max_rel_diff;        /* config: max. relative difference */
    double max_diff;            /* config: max. difference */
    int initial_range;          /* config: max. offset used as initial velocity */
    Bool average_accel;         /* config: average acceleration over velocity */
    PointerAccelerationProfileFunc Profile;
    PointerAccelerationProfileFunc deviceSpecificProfile;
    void *profile_private;      /* extended data, see  SetAccelerationProfile() */
    struct {                    /* to be able to query this information */
        int profile_number;
    } statistics;
} DeviceVelocityRec, *DeviceVelocityPtr;

/**
 * contains the run-time data for the predictable scheme, that is, a
 * DeviceVelocityPtr and the property handlers.
 */
typedef struct _PredictableAccelSchemeRec {
    DeviceVelocityPtr vel;
    long *prop_handlers;
    int num_prop_handlers;
} PredictableAccelSchemeRec, *PredictableAccelSchemePtr;

extern _X_EXPORT void
InitVelocityData(DeviceVelocityPtr vel);

extern _X_EXPORT void
InitTrackers(DeviceVelocityPtr vel, int ntracker);

extern _X_EXPORT BOOL
ProcessVelocityData2D(DeviceVelocityPtr vel, double dx, double dy, int time);

extern _X_EXPORT double
BasicComputeAcceleration(DeviceIntPtr dev, DeviceVelocityPtr vel,
                         double velocity, double threshold, double acc);

extern _X_EXPORT void
FreeVelocityData(DeviceVelocityPtr vel);

extern _X_EXPORT int
SetAccelerationProfile(DeviceVelocityPtr vel, int profile_num);

extern _X_EXPORT DeviceVelocityPtr
GetDevicePredictableAccelData(DeviceIntPtr dev);

extern _X_EXPORT void
SetDeviceSpecificAccelerationProfile(DeviceVelocityPtr vel,
                                     PointerAccelerationProfileFunc profile);

extern _X_INTERNAL void
AccelerationDefaultCleanup(DeviceIntPtr dev);

extern _X_INTERNAL Bool
InitPredictableAccelerationScheme(DeviceIntPtr dev,
                                  struct _ValuatorAccelerationRec *protoScheme);

extern _X_INTERNAL void
acceleratePointerPredictable(DeviceIntPtr dev, ValuatorMask *val,
                             CARD32 evtime);

extern _X_INTERNAL void
acceleratePointerLightweight(DeviceIntPtr dev, ValuatorMask *val,
                             CARD32 evtime);

//NR #endif                          /* POINTERVELOCITY_H */


// -------------------------------------------------------------------------
// xorg-server-23.2.1/os/backtrace.c (partial)
// -------------------------------------------------------------------------
/* lines 391-396 */
/* Default fallback if we can't find any way to get a backtrace */
void
xorg_backtrace(void)
{
    return;
}


// -------------------------------------------------------------------------
// xorg-server-23.2.1/os/inputthread.c
// -------------------------------------------------------------------------

/* lines 537-538 */
void input_lock(void) {}
void input_unlock(void) {}

// -------------------------------------------------------------------------
// xorg-server-23.2.1/dix/inpututils.c
// -------------------------------------------------------------------------

/* lines 693-704 */
int
CountBits(const uint8_t *mask, int len) // removed const to avoid warning
{
    int i;
    int ret = 0;

    for (i = 0; i < len; i++)
        if (BitIsOn(mask, i))
            ret++;

    return ret;
}

/* lines 422-472 */
/**
 * Alloc a valuator mask large enough for num_valuators.
 */
ValuatorMask *
valuator_mask_new(int num_valuators)
{
    /* alloc a fixed size mask for now and ignore num_valuators. in the
     * flying-car future, when we can dynamically alloc the masks and are
     * not constrained by signals, we can start using num_valuators */
    ValuatorMask *mask = (ValuatorMask *) calloc(1, sizeof(ValuatorMask)); // cast to avoid warning

    if (mask == NULL)
        return NULL;

    mask->last_bit = -1;
    return mask;
}

void
valuator_mask_free(ValuatorMask **mask)
{
    free(*mask);
    *mask = NULL;
}

/**
 * Sets a range of valuators between first_valuator and num_valuators with
 * the data in the valuators array. All other values are set to 0.
 */
void
valuator_mask_set_range(ValuatorMask *mask, int first_valuator,
                        int num_valuators, const int *valuators)
{
    int i;

    valuator_mask_zero(mask);

    for (i = first_valuator;
         i < _min(first_valuator + num_valuators, MAX_VALUATORS); i++) // changed to _min to avoid clash with cmath
        valuator_mask_set(mask, i, valuators[i - first_valuator]);
}

/**
 * Reset mask to zero.
 */
void
valuator_mask_zero(ValuatorMask *mask)
{
    memset(mask, 0, sizeof(*mask));
    mask->last_bit = -1;
}

/* lines 484-549 */
/**
 * Returns the number of valuators set in the given mask.
 */
int
valuator_mask_num_valuators(const ValuatorMask *mask)
{
    return CountBits((uint8_t *)(mask->mask), _min(mask->last_bit + 1, MAX_VALUATORS)); // changed to _min to avoid clash with cmath
}

/**
 * Return true if the valuator is set in the mask, or false otherwise.
 */
int
valuator_mask_isset(const ValuatorMask *mask, int valuator)
{
    return mask->last_bit >= valuator && BitIsOn(mask->mask, valuator);
}

static inline void
_valuator_mask_set_double(ValuatorMask *mask, int valuator, double data)
{
    mask->last_bit = _max(valuator, mask->last_bit); // changed to _max to avoid clash with cmath
    SetBit(mask->mask, valuator);
    mask->valuators[valuator] = data;
}

/**
 * Set the valuator to the given floating-point data.
 */
void
valuator_mask_set_double(ValuatorMask *mask, int valuator, double data)
{
    BUG_WARN_MSG(mask->has_unaccelerated,
                 "Do not mix valuator types, zero mask first\n");
    _valuator_mask_set_double(mask, valuator, data);
}

/**
 * Set the valuator to the given integer data.
 */
void
valuator_mask_set(ValuatorMask *mask, int valuator, int data)
{
    valuator_mask_set_double(mask, valuator, data);
}

/**
 * Return the requested valuator value as a double. If the mask bit is not
 * set for the given valuator, the returned value is undefined.
 */
double
valuator_mask_get_double(const ValuatorMask *mask, int valuator)
{
    return mask->valuators[valuator];
}

#include <math.h> // added here for trunc (and other functions used later on)
using namespace std;

/**
 * Return the requested valuator value as an integer, rounding towards zero.
 * If the mask bit is not set for the given valuator, the returned value is
 * undefined.
 */
int
valuator_mask_get(const ValuatorMask *mask, int valuator)
{
    return trunc(valuator_mask_get_double(mask, valuator));
}

// -------------------------------------------------------------------------
// xorg-server-23.2.1/dix/ptrveloc.c
// -------------------------------------------------------------------------
/*
 *
 * Copyright © 2006-2009 Simon Thum             simon dot thum at gmx dot de
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice (including the next
 * paragraph) shall be included in all copies or substantial portions of the
 * Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#if 0 //NR
#ifdef HAVE_DIX_CONFIG_H
#include <dix-config.h>
#endif
#endif //NR

#if 0 //NR
#include <math.h>
#include <ptrveloc.h>
#include <exevents.h>
#include <X11/Xatom.h>
#include <os.h>
#endif //NR

//NR #include <xserver-properties.h>

/*****************************************************************************
 * Predictable pointer acceleration
 *
 * 2006-2009 by Simon Thum (simon [dot] thum [at] gmx de)
 *
 * Serves 3 complementary functions:
 * 1) provide a sophisticated ballistic velocity estimate to improve
 *    the relation between velocity (of the device) and acceleration
 * 2) make arbitrary acceleration profiles possible
 * 3) decelerate by two means (constant and adaptive) if enabled
 *
 * Important concepts are the
 *
 * - Scheme
 *      which selects the basic algorithm
 *      (see devices.c/InitPointerAccelerationScheme)
 * - Profile
 *      which returns an acceleration
 *      for a given velocity
 *
 *  The profile can be selected by the user at runtime.
 *  The classic profile is intended to cleanly perform old-style
 *  function selection (threshold =/!= 0)
 *
 ****************************************************************************/

/* fwds */
static double
SimpleSmoothProfile(DeviceIntPtr dev, DeviceVelocityPtr vel, double velocity,
                    double threshold, double acc);
static PointerAccelerationProfileFunc
GetAccelerationProfile(DeviceVelocityPtr vel, int profile_num);
static BOOL
InitializePredictableAccelerationProperties(DeviceIntPtr,
                                            DeviceVelocityPtr,
                                            PredictableAccelSchemePtr);
static BOOL
DeletePredictableAccelerationProperties(DeviceIntPtr,
                                        PredictableAccelSchemePtr);


/*#define PTRACCEL_DEBUGGING*/

#if 0 //NR
#ifdef PTRACCEL_DEBUGGING
#define DebugAccelF(...) ErrorFSigSafe("dix/ptraccel: " __VA_ARGS__)
#else
#define DebugAccelF(...)        /* */
#endif
#endif //NR

/********************************
 *  Init/Uninit
 *******************************/

/* some int which is not a profile number */
#define PROFILE_UNINITIALIZE (-100)

/**
 * Init DeviceVelocity struct so it should match the average case
 */
void
InitVelocityData(DeviceVelocityPtr vel)
{
    memset(vel, 0, sizeof(DeviceVelocityRec));

    vel->corr_mul = 10.0;       /* dots per 10 millisecond should be usable */
    vel->const_acceleration = 1.0;      /* no acceleration/deceleration  */
    vel->reset_time = 300;
    vel->use_softening = 1;
    vel->min_acceleration = 1.0;        /* don't decelerate */
    vel->max_rel_diff = 0.2;
    vel->max_diff = 1.0;
    vel->initial_range = 2;
    vel->average_accel = TRUE;
    SetAccelerationProfile(vel, AccelProfileClassic);
    InitTrackers(vel, 16);
}

/**
 * Clean up DeviceVelocityRec
 */
void
FreeVelocityData(DeviceVelocityPtr vel)
{
    free(vel->tracker);
    SetAccelerationProfile(vel, PROFILE_UNINITIALIZE);
}

/**
 * Init predictable scheme
 */
Bool
InitPredictableAccelerationScheme(DeviceIntPtr dev,
                                  ValuatorAccelerationPtr protoScheme)
{
    DeviceVelocityPtr vel;
    ValuatorAccelerationRec scheme;
    PredictableAccelSchemePtr schemeData;

    scheme = *protoScheme;
    vel = (DeviceVelocityPtr) calloc(1, sizeof(DeviceVelocityRec)); // cast to avoid warning
    schemeData = (PredictableAccelSchemePtr) calloc(1, sizeof(PredictableAccelSchemeRec));  // cast to avoid warning
    if (!vel || !schemeData) {
        free(vel);
        free(schemeData);
        return FALSE;
    }
    InitVelocityData(vel);
    schemeData->vel = vel;
    scheme.accelData = schemeData;
    if (!InitializePredictableAccelerationProperties(dev, vel, schemeData)) {
        free(vel);
        free(schemeData);
        return FALSE;
    }
    /* all fine, assign scheme to device */
    dev->valuator->accelScheme = scheme;
    return TRUE;
}

/**
 *  Uninit scheme
 */
void
AccelerationDefaultCleanup(DeviceIntPtr dev)
{
    DeviceVelocityPtr vel = GetDevicePredictableAccelData(dev);

    if (vel) {
        /* the proper guarantee would be that we're not inside of
         * AccelSchemeProc(), but that seems impossible. Schemes don't get
         * switched often anyway.
         */
        input_lock();
        dev->valuator->accelScheme.AccelSchemeProc = NULL;
        FreeVelocityData(vel);
        free(vel);
        DeletePredictableAccelerationProperties(dev,
                                                (PredictableAccelSchemePtr)
                                                dev->valuator->accelScheme.
                                                accelData);
        free(dev->valuator->accelScheme.accelData);
        dev->valuator->accelScheme.accelData = NULL;
        input_unlock();
    }
}

/*************************
 * Input property support
 ************************/
#if 0 //NR
/**
 * choose profile
 */
static int
AccelSetProfileProperty(DeviceIntPtr dev, Atom atom,
                        XIPropertyValuePtr val, BOOL checkOnly)
{
    DeviceVelocityPtr vel;
    int profile, *ptr = &profile;
    int rc;
    int nelem = 1;

    if (atom != XIGetKnownProperty(ACCEL_PROP_PROFILE_NUMBER))
        return Success;

    vel = GetDevicePredictableAccelData(dev);
    if (!vel)
        return BadValue;
    rc = XIPropToInt(val, &nelem, &ptr);

    if (checkOnly) {
        if (rc)
            return rc;

        if (GetAccelerationProfile(vel, profile) == NULL)
            return BadValue;
    }
    else
        SetAccelerationProfile(vel, profile);

    return Success;
}

static long
AccelInitProfileProperty(DeviceIntPtr dev, DeviceVelocityPtr vel)
{
    int profile = vel->statistics.profile_number;
    Atom prop_profile_number = XIGetKnownProperty(ACCEL_PROP_PROFILE_NUMBER);

    XIChangeDeviceProperty(dev, prop_profile_number, XA_INTEGER, 32,
                           PropModeReplace, 1, &profile, FALSE);
    XISetDevicePropertyDeletable(dev, prop_profile_number, FALSE);
    return XIRegisterPropertyHandler(dev, AccelSetProfileProperty, NULL, NULL);
}

/**
 * constant deceleration
 */
static int
AccelSetDecelProperty(DeviceIntPtr dev, Atom atom,
                      XIPropertyValuePtr val, BOOL checkOnly)
{
    DeviceVelocityPtr vel;
    float v, *ptr = &v;
    int rc;
    int nelem = 1;

    if (atom != XIGetKnownProperty(ACCEL_PROP_CONSTANT_DECELERATION))
        return Success;

    vel = GetDevicePredictableAccelData(dev);
    if (!vel)
        return BadValue;
    rc = XIPropToFloat(val, &nelem, &ptr);

    if (checkOnly) {
        if (rc)
            return rc;
        return (v > 0) ? Success : BadValue;
    }

    vel->const_acceleration = 1 / v;

    return Success;
}

static long
AccelInitDecelProperty(DeviceIntPtr dev, DeviceVelocityPtr vel)
{
    float fval = 1.0 / vel->const_acceleration;
    Atom prop_const_decel =
        XIGetKnownProperty(ACCEL_PROP_CONSTANT_DECELERATION);
    XIChangeDeviceProperty(dev, prop_const_decel,
                           XIGetKnownProperty(XATOM_FLOAT), 32, PropModeReplace,
                           1, &fval, FALSE);
    XISetDevicePropertyDeletable(dev, prop_const_decel, FALSE);
    return XIRegisterPropertyHandler(dev, AccelSetDecelProperty, NULL, NULL);
}

/**
 * adaptive deceleration
 */
static int
AccelSetAdaptDecelProperty(DeviceIntPtr dev, Atom atom,
                           XIPropertyValuePtr val, BOOL checkOnly)
{
    DeviceVelocityPtr veloc;
    float v, *ptr = &v;
    int rc;
    int nelem = 1;

    if (atom != XIGetKnownProperty(ACCEL_PROP_ADAPTIVE_DECELERATION))
        return Success;

    veloc = GetDevicePredictableAccelData(dev);
    if (!veloc)
        return BadValue;
    rc = XIPropToFloat(val, &nelem, &ptr);

    if (checkOnly) {
        if (rc)
            return rc;
        return (v >= 1.0f) ? Success : BadValue;
    }

    if (v >= 1.0f)
        veloc->min_acceleration = 1 / v;

    return Success;
}

static long
AccelInitAdaptDecelProperty(DeviceIntPtr dev, DeviceVelocityPtr vel)
{
    float fval = 1.0 / vel->min_acceleration;
    Atom prop_adapt_decel =
        XIGetKnownProperty(ACCEL_PROP_ADAPTIVE_DECELERATION);

    XIChangeDeviceProperty(dev, prop_adapt_decel,
                           XIGetKnownProperty(XATOM_FLOAT), 32, PropModeReplace,
                           1, &fval, FALSE);
    XISetDevicePropertyDeletable(dev, prop_adapt_decel, FALSE);
    return XIRegisterPropertyHandler(dev, AccelSetAdaptDecelProperty, NULL,
                                     NULL);
}

/**
 * velocity scaling
 */
static int
AccelSetScaleProperty(DeviceIntPtr dev, Atom atom,
                      XIPropertyValuePtr val, BOOL checkOnly)
{
    DeviceVelocityPtr vel;
    float v, *ptr = &v;
    int rc;
    int nelem = 1;

    if (atom != XIGetKnownProperty(ACCEL_PROP_VELOCITY_SCALING))
        return Success;

    vel = GetDevicePredictableAccelData(dev);
    if (!vel)
        return BadValue;
    rc = XIPropToFloat(val, &nelem, &ptr);

    if (checkOnly) {
        if (rc)
            return rc;

        return (v > 0) ? Success : BadValue;
    }

    if (v > 0)
        vel->corr_mul = v;

    return Success;
}

static long
AccelInitScaleProperty(DeviceIntPtr dev, DeviceVelocityPtr vel)
{
    float fval = vel->corr_mul;
    Atom prop_velo_scale = XIGetKnownProperty(ACCEL_PROP_VELOCITY_SCALING);

    XIChangeDeviceProperty(dev, prop_velo_scale,
                           XIGetKnownProperty(XATOM_FLOAT), 32, PropModeReplace,
                           1, &fval, FALSE);
    XISetDevicePropertyDeletable(dev, prop_velo_scale, FALSE);
    return XIRegisterPropertyHandler(dev, AccelSetScaleProperty, NULL, NULL);
}
#endif //NR

static BOOL
InitializePredictableAccelerationProperties(DeviceIntPtr dev,
                                            DeviceVelocityPtr vel,
                                            PredictableAccelSchemePtr
                                            schemeData)
{
    int num_handlers = 4;

    if (!vel)
        return FALSE;

    schemeData->prop_handlers = (long *) calloc(num_handlers, sizeof(long)); // cast to avoid warning
    if (!schemeData->prop_handlers)
        return FALSE;
    schemeData->num_prop_handlers = num_handlers;
#if 0 //NR
    schemeData->prop_handlers[0] = AccelInitProfileProperty(dev, vel);
    schemeData->prop_handlers[1] = AccelInitDecelProperty(dev, vel);
    schemeData->prop_handlers[2] = AccelInitAdaptDecelProperty(dev, vel);
    schemeData->prop_handlers[3] = AccelInitScaleProperty(dev, vel);
#endif //NR
    return TRUE;
}

BOOL
DeletePredictableAccelerationProperties(DeviceIntPtr dev,
                                        PredictableAccelSchemePtr scheme)
{
#if 0 //NR
    DeviceVelocityPtr vel;
    Atom prop;
    int i;

    prop = XIGetKnownProperty(ACCEL_PROP_VELOCITY_SCALING);
    XIDeleteDeviceProperty(dev, prop, FALSE);
    prop = XIGetKnownProperty(ACCEL_PROP_ADAPTIVE_DECELERATION);
    XIDeleteDeviceProperty(dev, prop, FALSE);
    prop = XIGetKnownProperty(ACCEL_PROP_CONSTANT_DECELERATION);
    XIDeleteDeviceProperty(dev, prop, FALSE);
    prop = XIGetKnownProperty(ACCEL_PROP_PROFILE_NUMBER);
    XIDeleteDeviceProperty(dev, prop, FALSE);

    vel = GetDevicePredictableAccelData(dev);
    if (vel) {
        for (i = 0; i < scheme->num_prop_handlers; i++)
            if (scheme->prop_handlers[i])
                XIUnregisterPropertyHandler(dev, scheme->prop_handlers[i]);
    }

    free(scheme->prop_handlers);
    scheme->prop_handlers = NULL;
    scheme->num_prop_handlers = 0;
#endif //NR
    return TRUE;
}

/*********************
 * Tracking logic
 ********************/

void
InitTrackers(DeviceVelocityPtr vel, int ntracker)
{
    if (ntracker < 1) {
        ErrorF("invalid number of trackers\n");
        return;
    }
    free(vel->tracker);
    vel->tracker = (MotionTrackerPtr) calloc(ntracker, sizeof(MotionTracker));
    vel->num_tracker = ntracker;
}

enum directions {
    N = (1 << 0),
    NE = (1 << 1),
    E = (1 << 2),
    SE = (1 << 3),
    S = (1 << 4),
    SW = (1 << 5),
    W = (1 << 6),
    NW = (1 << 7),
    UNDEFINED = 0xFF
};

/**
 * return a bit field of possible directions.
 * There's no reason against widening to more precise directions (<45 degrees),
 * should it not perform well. All this is needed for is sort out non-linear
 * motion, so precision isn't paramount. However, one should not flag direction
 * too narrow, since it would then cut the linear segment to zero size way too
 * often.
 *
 * @return A bitmask for N, NE, S, SE, etc. indicating the directions for
 * this movement.
 */
static int
DoGetDirection(int dx, int dy)
{
    int dir = 0;

    /* on insignificant mickeys, flag 135 degrees */
    if (abs(dx) < 2 && abs(dy) < 2) {
        /* first check diagonal cases */
        if (dx > 0 && dy > 0)
            dir = E | SE | S;
        else if (dx > 0 && dy < 0)
            dir = N | NE | E;
        else if (dx < 0 && dy < 0)
            dir = W | NW | N;
        else if (dx < 0 && dy > 0)
            dir = W | SW | S;
        /* check axis-aligned directions */
        else if (dx > 0)
            dir = NE | E | SE;
        else if (dx < 0)
            dir = NW | W | SW;
        else if (dy > 0)
            dir = SE | S | SW;
        else if (dy < 0)
            dir = NE | N | NW;
        else
            dir = UNDEFINED;    /* shouldn't happen */
    }
    else {                      /* compute angle and set appropriate flags */
        double r;
        int i1, i2;

        r = atan2((double) dy, (double) dx); // cast to avoid warning
        /* find direction.
         *
         * Add 360° to avoid r become negative since C has no well-defined
         * modulo for such cases. Then divide by 45° to get the octant
         * number,  e.g.
         *          0 <= r <= 1 is [0-45]°
         *          1 <= r <= 2 is [45-90]°
         *          etc.
         * But we add extra 90° to match up with our N, S, etc. defines up
         * there, rest stays the same.
         */
        r = (r + (M_PI * 2.5)) / (M_PI / 4);
        /* this intends to flag 2 directions (45 degrees),
         * except on very well-aligned mickeys. */
        i1 = (int) (r + 0.1) % 8;
        i2 = (int) (r + 0.9) % 8;
        if (i1 < 0 || i1 > 7 || i2 < 0 || i2 > 7)
            dir = UNDEFINED;    /* shouldn't happen */
        else
            dir = (1 << i1 | 1 << i2);
    }
    return dir;
}

#define DIRECTION_CACHE_RANGE 5
#define DIRECTION_CACHE_SIZE (DIRECTION_CACHE_RANGE*2+1)

/* cache DoGetDirection().
 * To avoid excessive use of direction calculation, cache the values for
 * [-5..5] for both x/y. Anything outside of that is calculated on the fly.
 *
 * @return A bitmask for N, NE, S, SE, etc. indicating the directions for
 * this movement.
 */
static int
GetDirection(int dx, int dy)
{
    static int cache[DIRECTION_CACHE_SIZE][DIRECTION_CACHE_SIZE];
    int dir;

    if (abs(dx) <= DIRECTION_CACHE_RANGE && abs(dy) <= DIRECTION_CACHE_RANGE) {
        /* cacheable */
        dir = cache[DIRECTION_CACHE_RANGE + dx][DIRECTION_CACHE_RANGE + dy];
        if (dir == 0) {
            dir = DoGetDirection(dx, dy);
            cache[DIRECTION_CACHE_RANGE + dx][DIRECTION_CACHE_RANGE + dy] = dir;
        }
    }
    else {
        /* non-cacheable */
        dir = DoGetDirection(dx, dy);
    }

    return dir;
}

#undef DIRECTION_CACHE_RANGE
#undef DIRECTION_CACHE_SIZE

/* convert offset (age) to array index */
#define TRACKER_INDEX(s, d) (((s)->num_tracker + (s)->cur_tracker - (d)) % (s)->num_tracker)
#define TRACKER(s, d) &(s)->tracker[TRACKER_INDEX(s,d)]

/**
 * Add the delta motion to each tracker, then reset the latest tracker to
 * 0/0 and set it as the current one.
 */
static inline void
FeedTrackers(DeviceVelocityPtr vel, double dx, double dy, int cur_t)
{
    int n;

    for (n = 0; n < vel->num_tracker; n++) {
        vel->tracker[n].dx += dx;
        vel->tracker[n].dy += dy;
    }
    n = (vel->cur_tracker + 1) % vel->num_tracker;
    vel->tracker[n].dx = 0.0;
    vel->tracker[n].dy = 0.0;
    vel->tracker[n].time = cur_t;
    vel->tracker[n].dir = GetDirection(dx, dy);
    DebugAccelF("motion [dx: %f dy: %f dir:%d diff: %d]\n",
                dx, dy, vel->tracker[n].dir,
                cur_t - vel->tracker[vel->cur_tracker].time);
    vel->cur_tracker = n;
}

/**
 * calc velocity for given tracker, with
 * velocity scaling.
 * This assumes linear motion.
 */
static double
CalcTracker(const MotionTracker * tracker, int cur_t)
{
    double dist = sqrt(tracker->dx * tracker->dx + tracker->dy * tracker->dy);
    int dtime = cur_t - tracker->time;

    if (dtime > 0)
        return dist / dtime;
    else
        return 0;               /* synonymous for NaN, since we're not C99 */
}

/* find the most plausible velocity. That is, the most distant
 * (in time) tracker which isn't too old, the movement vector was
 * in the same octant, and where the velocity is within an
 * acceptable range to the initial velocity.
 *
 * @return The tracker's velocity or 0 if the above conditions are unmet
 */
static double
QueryTrackers(DeviceVelocityPtr vel, int cur_t)
{
    int offset, dir = UNDEFINED, used_offset = -1, age_ms;

    /* initial velocity: a low-offset, valid velocity */
    double initial_velocity = 0, result = 0, velocity_diff;
    double velocity_factor = vel->corr_mul * vel->const_acceleration;   /* premultiply */

    /* loop from current to older data */
    for (offset = 1; offset < vel->num_tracker; offset++) {
        MotionTracker *tracker = TRACKER(vel, offset);
        double tracker_velocity;

        age_ms = cur_t - tracker->time;

        /* bail out if data is too old and protect from overrun */
        if (age_ms >= vel->reset_time || age_ms < 0) {
            DebugAccelF("query: tracker too old (reset after %d, age is %d)\n",
                        vel->reset_time, age_ms);
            break;
        }

        /*
         * this heuristic avoids using the linear-motion velocity formula
         * in CalcTracker() on motion that isn't exactly linear. So to get
         * even more precision we could subdivide as a final step, so possible
         * non-linearities are accounted for.
         */
        dir &= tracker->dir;
        if (dir == 0) {         /* we've changed octant of movement (e.g. NE → NW) */
            DebugAccelF("query: no longer linear\n");
            /* instead of breaking it we might also inspect the partition after,
             * but actual improvement with this is probably rare. */
            break;
        }

        tracker_velocity = CalcTracker(tracker, cur_t) * velocity_factor;

        if ((initial_velocity == 0 || offset <= vel->initial_range) &&
            tracker_velocity != 0) {
            /* set initial velocity and result */
            result = initial_velocity = tracker_velocity;
            used_offset = offset;
        }
        else if (initial_velocity != 0 && tracker_velocity != 0) {
            velocity_diff = fabs(initial_velocity - tracker_velocity);

            if (velocity_diff > vel->max_diff &&
                velocity_diff / (initial_velocity + tracker_velocity) >=
                vel->max_rel_diff) {
                /* we're not in range, quit - it won't get better. */
                DebugAccelF("query: tracker too different:"
                            " old %2.2f initial %2.2f diff: %2.2f\n",
                            tracker_velocity, initial_velocity, velocity_diff);
                break;
            }
            /* we're in range with the initial velocity,
             * so this result is likely better
             * (it contains more information). */
            result = tracker_velocity;
            used_offset = offset;
        }
    }
    if (offset == vel->num_tracker) {
        DebugAccelF("query: last tracker in effect\n");
        used_offset = vel->num_tracker - 1;
    }
    if (used_offset >= 0) {
#ifdef PTRACCEL_DEBUGGING
        MotionTracker *tracker = TRACKER(vel, used_offset);

        DebugAccelF("result: offset %i [dx: %f dy: %f diff: %i]\n",
                    used_offset, tracker->dx, tracker->dy,
                    cur_t - tracker->time);
#endif
    }
    return result;
}

#undef TRACKER_INDEX
#undef TRACKER

/**
 * Perform velocity approximation based on 2D 'mickeys' (mouse motion delta).
 * return true if non-visible state reset is suggested
 */
BOOL
ProcessVelocityData2D(DeviceVelocityPtr vel, double dx, double dy, int time)
{
    double velocity;

    vel->last_velocity = vel->velocity;

    FeedTrackers(vel, dx, dy, time);

    velocity = QueryTrackers(vel, time);

    DebugAccelF("velocity is %f\n", velocity);

    vel->velocity = velocity;
    return velocity == 0;
}

/**
 * this flattens significant ( > 1) mickeys a little bit for more steady
 * constant-velocity response
 */
static inline double
ApplySimpleSoftening(double prev_delta, double delta)
{
    double result = delta;

    if (delta < -1.0 || delta > 1.0) {
        if (delta > prev_delta)
            result -= 0.5;
        else if (delta < prev_delta)
            result += 0.5;
    }
    return result;
}

/**
 * Soften the delta based on previous deltas stored in vel.
 *
 * @param[in,out] fdx Delta X, modified in-place.
 * @param[in,out] fdx Delta Y, modified in-place.
 */
static void
ApplySoftening(DeviceVelocityPtr vel, double *fdx, double *fdy)
{
    if (vel->use_softening) {
        *fdx = ApplySimpleSoftening(vel->last_dx, *fdx);
        *fdy = ApplySimpleSoftening(vel->last_dy, *fdy);
    }
}

static void
ApplyConstantDeceleration(DeviceVelocityPtr vel, double *fdx, double *fdy)
{
    *fdx *= vel->const_acceleration;
    *fdy *= vel->const_acceleration;
}

/*
 * compute the acceleration for given velocity and enforce min_acceleration
 */
double
BasicComputeAcceleration(DeviceIntPtr dev,
                         DeviceVelocityPtr vel,
                         double velocity, double threshold, double acc)
{

    double result;

    result = vel->Profile(dev, vel, velocity, threshold, acc);

    /* enforce min_acceleration */
    if (result < vel->min_acceleration)
        result = vel->min_acceleration;
    return result;
}

/**
 * Compute acceleration. Takes into account averaging, nv-reset, etc.
 * If the velocity has changed, an average is taken of 6 velocity factors:
 * current velocity, last velocity and 4 times the average between the two.
 */
static double
ComputeAcceleration(DeviceIntPtr dev,
                    DeviceVelocityPtr vel, double threshold, double acc)
{
    double result;

    if (vel->velocity <= 0) {
        DebugAccelF("profile skipped\n");
        /*
         * If we have no idea about device velocity, don't pretend it.
         */
        return 1;
    }

    if (vel->average_accel && vel->velocity != vel->last_velocity) {
        /* use simpson's rule to average acceleration between
         * current and previous velocity.
         * Though being the more natural choice, it causes a minor delay
         * in comparison, so it can be disabled. */
        result =
            BasicComputeAcceleration(dev, vel, vel->velocity, threshold, acc);
        result +=
            BasicComputeAcceleration(dev, vel, vel->last_velocity, threshold,
                                     acc);
        result +=
            4.0 * BasicComputeAcceleration(dev, vel,
                                            (vel->last_velocity +
                                             vel->velocity) / 2,
                                            threshold,
                                            acc);
        result /= 6.0;
        DebugAccelF("profile average [%.2f ... %.2f] is %.3f\n",
                    vel->velocity, vel->last_velocity, result);
    }
    else {
        result = BasicComputeAcceleration(dev, vel,
                                          vel->velocity, threshold, acc);
        DebugAccelF("profile sample [%.2f] is %.3f\n",
                    vel->velocity, result);
    }

    return result;
}

/*****************************************
 *  Acceleration functions and profiles
 ****************************************/

/**
 * Polynomial function similar previous one, but with f(1) = 1
 */
static double
PolynomialAccelerationProfile(DeviceIntPtr dev,
                              DeviceVelocityPtr vel,
                              double velocity, double ignored, double acc)
{
    return pow(velocity, (acc - 1.0) * 0.5);
}

/**
 * returns acceleration for velocity.
 * This profile selects the two functions like the old scheme did
 */
static double
ClassicProfile(DeviceIntPtr dev,
               DeviceVelocityPtr vel,
               double velocity, double threshold, double acc)
{
    if (threshold > 0) {
        return SimpleSmoothProfile(dev, vel, velocity, threshold, acc);
    }
    else {
        return PolynomialAccelerationProfile(dev, vel, velocity, 0, acc);
    }
}

/**
 * Power profile
 * This has a completely smooth transition curve, i.e. no jumps in the
 * derivatives.
 *
 * This has the expense of overall response dependency on min-acceleration.
 * In effect, min_acceleration mimics const_acceleration in this profile.
 */
static double
PowerProfile(DeviceIntPtr dev,
             DeviceVelocityPtr vel,
             double velocity, double threshold, double acc)
{
    double vel_dist;

    acc = (acc - 1.0) * 0.1 + 1.0;     /* without this, acc of 2 is unuseable */

    if (velocity <= threshold)
        return vel->min_acceleration;
    vel_dist = velocity - threshold;
    return (pow(acc, vel_dist)) * vel->min_acceleration;
}

/**
 * just a smooth function in [0..1] -> [0..1]
 *  - point symmetry at 0.5
 *  - f'(0) = f'(1) = 0
 *  - starts faster than a sinoid
 *  - smoothness C1 (Cinf if you dare to ignore endpoints)
 */
static inline double
CalcPenumbralGradient(double x)
{
    x *= 2.0;
    x -= 1.0;
    return 0.5 + (x * sqrt(1.0 - x * x) + asin(x)) / M_PI;
}

/**
 * acceleration function similar to classic accelerated/unaccelerated,
 * but with smooth transition in between (and towards zero for adaptive dec.).
 */
static double
SimpleSmoothProfile(DeviceIntPtr dev,
                    DeviceVelocityPtr vel,
                    double velocity, double threshold, double acc)
{
    if (velocity < 1.0f)
        return CalcPenumbralGradient(0.5 + velocity * 0.5) * 2.0f - 1.0f;
    if (threshold < 1.0f)
        threshold = 1.0f;
    if (velocity <= threshold)
        return 1;
    velocity /= threshold;
    if (velocity >= acc)
        return acc;
    else
        return 1.0f + (CalcPenumbralGradient(velocity / acc) * (acc - 1.0f));
}

/**
 * This profile uses the first half of the penumbral gradient as a start
 * and then scales linearly.
 */
static double
SmoothLinearProfile(DeviceIntPtr dev,
                    DeviceVelocityPtr vel,
                    double velocity, double threshold, double acc)
{
    double res, nv;

    if (acc > 1.0)
        acc -= 1.0;            /*this is so acc = 1 is no acceleration */
    else
        return 1.0;

    nv = (velocity - threshold) * acc * 0.5;

    if (nv < 0) {
        res = 0;
    }
    else if (nv < 2) {
        res = CalcPenumbralGradient(nv * 0.25) * 2.0;
    }
    else {
        nv -= 2.0;
        res = nv * 2.0 / M_PI  /* steepness of gradient at 0.5 */
            + 1.0;             /* gradient crosses 2|1 */
    }
    res += vel->min_acceleration;
    return res;
}

/**
 * From 0 to threshold, the response graduates smoothly from min_accel to
 * acceleration. Beyond threshold it is exactly the specified acceleration.
 */
static double
SmoothLimitedProfile(DeviceIntPtr dev,
                     DeviceVelocityPtr vel,
                     double velocity, double threshold, double acc)
{
    double res;

    if (velocity >= threshold || threshold == 0.0)
        return acc;

    velocity /= threshold;      /* should be [0..1[ now */

    res = CalcPenumbralGradient(velocity) * (acc - vel->min_acceleration);

    return vel->min_acceleration + res;
}

static double
LinearProfile(DeviceIntPtr dev,
              DeviceVelocityPtr vel,
              double velocity, double threshold, double acc)
{
    return acc * velocity;
}

static double
NoProfile(DeviceIntPtr dev,
          DeviceVelocityPtr vel, double velocity, double threshold, double acc)
{
    return 1.0;
}

static PointerAccelerationProfileFunc
GetAccelerationProfile(DeviceVelocityPtr vel, int profile_num)
{
    switch (profile_num) {
    case AccelProfileClassic:
        return ClassicProfile;
    case AccelProfileDeviceSpecific:
        return vel->deviceSpecificProfile;
    case AccelProfilePolynomial:
        return PolynomialAccelerationProfile;
    case AccelProfileSmoothLinear:
        return SmoothLinearProfile;
    case AccelProfileSimple:
        return SimpleSmoothProfile;
    case AccelProfilePower:
        return PowerProfile;
    case AccelProfileLinear:
        return LinearProfile;
    case AccelProfileSmoothLimited:
        return SmoothLimitedProfile;
    case AccelProfileNone:
        return NoProfile;
    default:
        return NULL;
    }
}

/**
 * Set the profile by number.
 * Intended to make profiles exchangeable at runtime.
 * If you created a profile, give it a number here and in the header to
 * make it selectable. In case some profile-specific init is needed, here
 * would be a good place, since FreeVelocityData() also calls this with
 * PROFILE_UNINITIALIZE.
 *
 * returns FALSE if profile number is unavailable, TRUE otherwise.
 */
int
SetAccelerationProfile(DeviceVelocityPtr vel, int profile_num)
{
    PointerAccelerationProfileFunc profile;

    profile = GetAccelerationProfile(vel, profile_num);

    if (profile == NULL && profile_num != PROFILE_UNINITIALIZE)
        return FALSE;

    /* Here one could free old profile-private data */
    free(vel->profile_private);
    vel->profile_private = NULL;
    /* Here one could init profile-private data */
    vel->Profile = profile;
    vel->statistics.profile_number = profile_num;
    return TRUE;
}

/**********************************************
 * driver interaction
 **********************************************/

/**
 * device-specific profile
 *
 * The device-specific profile is intended as a hook for a driver
 * which may want to provide an own acceleration profile.
 * It should not rely on profile-private data, instead
 * it should do init/uninit in the driver (ie. with DEVICE_INIT and friends).
 * Users may override or choose it.
 */
void
SetDeviceSpecificAccelerationProfile(DeviceVelocityPtr vel,
                                     PointerAccelerationProfileFunc profile)
{
    if (vel)
        vel->deviceSpecificProfile = profile;
}

/**
 * Use this function to obtain a DeviceVelocityPtr for a device. Will return NULL if
 * the predictable acceleration scheme is not in effect.
 */
DeviceVelocityPtr
GetDevicePredictableAccelData(DeviceIntPtr dev)
{
    BUG_RETURN_VAL(!dev, NULL);

    if (dev->valuator &&
        dev->valuator->accelScheme.AccelSchemeProc ==
        acceleratePointerPredictable &&
        dev->valuator->accelScheme.accelData != NULL) {

        return ((PredictableAccelSchemePtr)
                dev->valuator->accelScheme.accelData)->vel;
    }
    return NULL;
}

/********************************
 *  acceleration schemes
 *******************************/

/**
 * Modifies valuators in-place.
 * This version employs a velocity approximation algorithm to
 * enable fine-grained predictable acceleration profiles.
 */
void
acceleratePointerPredictable(DeviceIntPtr dev, ValuatorMask *val, CARD32 evtime)
{
    double dx = 0, dy = 0;
    DeviceVelocityPtr velocitydata = GetDevicePredictableAccelData(dev);
    Bool soften = TRUE;

    if (valuator_mask_num_valuators(val) == 0 || !velocitydata)
        return;

    if (velocitydata->statistics.profile_number == AccelProfileNone &&
        velocitydata->const_acceleration == 1.0) {
        return;                 /*we're inactive anyway, so skip the whole thing. */
    }

    if (valuator_mask_isset(val, 0)) {
        dx = valuator_mask_get_double(val, 0);
    }

    if (valuator_mask_isset(val, 1)) {
        dy = valuator_mask_get_double(val, 1);
    }

    if (dx != 0.0 || dy != 0.0) {
        /* reset non-visible state? */
        if (ProcessVelocityData2D(velocitydata, dx, dy, evtime)) {
            soften = FALSE;
        }

        if (dev->ptrfeed && dev->ptrfeed->ctrl.num) {
            double mult;

            /* invoke acceleration profile to determine acceleration */
            mult = ComputeAcceleration(dev, velocitydata,
                                       dev->ptrfeed->ctrl.threshold,
                                       (double) dev->ptrfeed->ctrl.num /
                                       (double) dev->ptrfeed->ctrl.den);

            DebugAccelF("mult is %f\n", mult);
            if (mult != 1.0 || velocitydata->const_acceleration != 1.0) {
                if (mult > 1.0 && soften)
                    ApplySoftening(velocitydata, &dx, &dy);
                ApplyConstantDeceleration(velocitydata, &dx, &dy);

                if (dx != 0.0)
                    valuator_mask_set_double(val, 0, mult * dx);
                if (dy != 0.0)
                    valuator_mask_set_double(val, 1, mult * dy);
                DebugAccelF("delta x:%.3f y:%.3f\n", mult * dx, mult * dy);
            }
        }
    }
    /* remember last motion delta (for softening/slow movement treatment) */
    velocitydata->last_dx = dx;
    velocitydata->last_dy = dy;
}

/**
 * Originally a part of xf86PostMotionEvent; modifies valuators
 * in-place. Retained mostly for embedded scenarios.
 */
void
acceleratePointerLightweight(DeviceIntPtr dev,
                             ValuatorMask *val, CARD32 ignored)
{
    double mult = 0.0, tmpf;
    double dx = 0.0, dy = 0.0;

    if (valuator_mask_isset(val, 0)) {
        dx = valuator_mask_get(val, 0);
    }

    if (valuator_mask_isset(val, 1)) {
        dy = valuator_mask_get(val, 1);
    }

    if (valuator_mask_num_valuators(val) == 0)
        return;

    if (dev->ptrfeed && dev->ptrfeed->ctrl.num) {
        /* modeled from xf86Events.c */
        if (dev->ptrfeed->ctrl.threshold) {
            if ((fabs(dx) + fabs(dy)) >= dev->ptrfeed->ctrl.threshold) {
                if (dx != 0.0) {
                    tmpf = (dx * (double) (dev->ptrfeed->ctrl.num)) /
                        (double) (dev->ptrfeed->ctrl.den);
                    valuator_mask_set_double(val, 0, tmpf);
                }

                if (dy != 0.0) {
                    tmpf = (dy * (double) (dev->ptrfeed->ctrl.num)) /
                        (double) (dev->ptrfeed->ctrl.den);
                    valuator_mask_set_double(val, 1, tmpf);
                }
            }
        }
        else {
            mult = pow(dx * dx + dy * dy,
                       ((double) (dev->ptrfeed->ctrl.num) /
                        (double) (dev->ptrfeed->ctrl.den) - 1.0) / 2.0) / 2.0;
            if (dx != 0.0)
                valuator_mask_set_double(val, 0, mult * dx);
            if (dy != 0.0)
                valuator_mask_set_double(val, 1, mult * dy);
        }
    }
}


// -------------------------------------------------------------------------
// xorg-server-23.2.1/dix/devices.c (partial)
// -------------------------------------------------------------------------

/* lines 1380-1387 */
/* global list of acceleration schemes */
ValuatorAccelerationRec pointerAccelerationScheme[] = {
    {PtrAccelNoOp, NULL, NULL, NULL, NULL},
    {PtrAccelPredictable, acceleratePointerPredictable, NULL,
     InitPredictableAccelerationScheme, AccelerationDefaultCleanup},
    {PtrAccelLightweight, acceleratePointerLightweight, NULL, NULL, NULL},
    {-1, NULL, NULL, NULL, NULL}        /* terminator */
};

/* lines 1389-1430 */
/**
 * install an acceleration scheme. returns TRUE on success, and should not
 * change anything if unsuccessful.
 */
Bool
InitPointerAccelerationScheme(DeviceIntPtr dev, int scheme)
{
    int x, i = -1;
    ValuatorClassPtr val;

    val = dev->valuator;

    if (!val)
        return FALSE;
#if 0 //NR
    if (IsMaster(dev) && scheme != PtrAccelNoOp)
        return FALSE;
#endif //NR
    for (x = 0; pointerAccelerationScheme[x].number >= 0; x++) {
        if (pointerAccelerationScheme[x].number == scheme) {
            i = x;
            break;
        }
    }

    if (-1 == i)
        return FALSE;

    if (val->accelScheme.AccelCleanupProc)
        val->accelScheme.AccelCleanupProc(dev);

    if (pointerAccelerationScheme[i].AccelInitProc) {
        if (!pointerAccelerationScheme[i].AccelInitProc(dev,
                                            &pointerAccelerationScheme[i])) {
            return FALSE;
        }
    }
    else {
        val->accelScheme = pointerAccelerationScheme[i];
    }
    return TRUE;
}
