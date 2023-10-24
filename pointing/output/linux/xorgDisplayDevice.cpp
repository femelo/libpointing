/* -*- mode: c++ -*-
 *
 * pointing/output/linux/xorgDisplayDevice.h --
 *
 * Initial software
 * Authors: Nicolas Roussel
 * Copyright Inria
 *
 * http://libpointing.org/
 *
 * This software may be used and distributed according to the terms of
 * the GNU General Public License version 2 or any later version.
 *
 */

#include <pointing/output/linux/xorgDisplayDevice.h>

#include <X11/extensions/Xrandr.h>

#include <iostream>
#include <sstream>
#include <stdexcept>

namespace pointing {

xorgDisplayDevice::xorgDisplayDevice(void) {
    displayID = -1;
    initialize();
}

xorgDisplayDevice::xorgDisplayDevice(URI uri) {
    displayID = -1;
    
    if (!uri.path.empty()) {
        std::istringstream to_int ;
        to_int.str(uri.path.erase(0,1)) ; // Skip the leading '/'
        to_int >> displayID ;
    }
    
    initialize();
}

XRRModeInfo* xorgDisplayDevice::get_mode_info(RRMode mode) {
    for (int m = 0; m < res->nmode; ++m) {
        XRRModeInfo *mode_info = &res->modes[m];
        
        if (mode_info->id == mode) {
            return mode_info;
        }
    }
    
    return NULL;
}

bool xorgDisplayDevice::get_bounds(int* width, int* height, int* x, int* y) {
    // Need to get the crtc info
    RRCrtc crtc = output_info->crtc;
    XRRCrtcInfo* crtc_info = XRRGetCrtcInfo(dpy, res, crtc);
    if (crtc_info != NULL) {
        *width = crtc_info->width;
        *height = crtc_info->height;
        *x = crtc_info->x;
        *y = crtc_info->y;
        return true;
    } else {
        // This should not happen; but if it does, we just assume that the device use the preferred mode
        RRMode mode_id = output_info->modes[output_info->npreferred];
        XRRModeInfo* mode_info = get_mode_info(mode_id);
        
        if (mode_info != NULL) {
            *width = mode_info->width;
            *height = mode_info->height;
            return true;
        }
    }
    
    XRRFreeCrtcInfo(crtc_info);
    
    return false;
}

void xorgDisplayDevice::initialize() {
    dpy = XOpenDisplay(0);
    
    if (dpy == NULL) {
        throw std::runtime_error("xorgDisplayDevice: can't open display");
    }
    
    screen = DefaultScreen(dpy);
    root = RootWindow(dpy, screen);
    
    res = XRRGetScreenResources(dpy, root);
    
    if (res == NULL) {
        throw std::runtime_error("xorgDisplayDevice: can't retrieve screen resources") ;
    }
    
    if (displayID == -1) {
        // This happen when asking for "any" display
        displayID = get_any_display_id();
    }
    
    if (res->noutput <= displayID) {
        throw std::runtime_error("xorgDisplayDevice: invalid display ID") ;
    }
    
    output_info = XRRGetOutputInfo(dpy, res, res->outputs[displayID]);
    
    if (output_info == NULL) {
        throw std::runtime_error("xorgDisplayDevice: can't retrieve output info") ;
    }
    cached = NOTHING;
}


int xorgDisplayDevice::get_any_display_id() {
     for (int o = 0; o < res->noutput; o++) {
         XRROutputInfo *output_info = XRRGetOutputInfo(dpy, res, res->outputs[o]);
         bool connected = output_info->connection == RR_Connected;
         XRRFreeOutputInfo(output_info);
         
         // Any display but among the connected ones
         if (connected) {
             return o;
         }
     }
     
     return -1;
}

DisplayDevice::Bounds
xorgDisplayDevice::getBounds(Bounds */*defval*/) {
    if (cached & BOUNDS) return cached_bounds ;
    
    int width = 0;
    int height = 0;
    int x = 0;
    int y = 0;
    
    get_bounds(&width, &height, &x, &y);
    cached_bounds = DisplayDevice::Bounds(x, y, width, height);
    
    cached = cached | BOUNDS ;
    
    return cached_bounds ;
}

DisplayDevice::Size
xorgDisplayDevice::getSize(Size */*defval*/) {
    if (cached & SIZE) return cached_size ;
    
    cached_size = DisplayDevice::Size(output_info->mm_width, output_info->mm_height) ;
    cached = cached | SIZE ;
    
    return cached_size ;
}

double
xorgDisplayDevice::getRefreshRate(double */*defval*/) {
    if (cached & REFRESHRATE) return cached_refreshrate ;
    
    XRRScreenConfiguration *sc = XRRGetScreenInfo(dpy, root) ;
    cached_refreshrate = XRRConfigCurrentRate(sc) ;
    cached = cached | REFRESHRATE ;
    
    return cached_refreshrate ;
}

URI
xorgDisplayDevice::getURI(bool /*expanded*/) const {
    URI uri ;
    uri.scheme = "xorgdisplay" ;
    std::stringstream path ;
    path << "/" << displayID ;
    uri.path = path.str() ;
    return uri ;
}

xorgDisplayDevice::~xorgDisplayDevice(void) {
    XRRFreeScreenResources(res);
    XRRFreeOutputInfo(output_info);
    XCloseDisplay(dpy);
}

}
