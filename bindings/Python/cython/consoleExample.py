# -*- coding: utf-8 -*-

from libpointing import PointingDevice, DisplayDevice, TransferFunction
from libpointing import PointingDeviceManager, PointingDeviceDescriptor

import sys

def cb_man(desc, wasAdded):
	print(desc)
	print("was added" if wasAdded else "was removed")


pm = PointingDeviceManager()
PointingDevice.idle(100)
pm.addDeviceUpdateCallback(cb_man)

for desc in pm:
	print(desc)

"""
for desc in pm:
	print desc.devURI
	print desc.vendor, desc.product
	pdev = PointingDevice(desc.uri)
"""

try:
    pdev = PointingDevice(f"input:{sys.argv[1]}".encode('utf-8'))
except:
    pdev = PointingDevice(b"any:")

ddev = DisplayDevice.create("any:")
tfct = TransferFunction(b"system:", pdev, ddev)

def cb_fct(timestamp, dx, dy, button):
    rx,ry=tfct.applyd(dx, dy, timestamp)
    print("%s: %d %d %d -> %.2f %.2f"%(str(timestamp), dx, dy, button, rx, ry ))
    sys.stdout.flush()

pdev.setCallback(cb_fct)
print("Move the mouse of Press CTRL+C to exit")
for i in range(0, 10000):
    PointingDevice.idle(1)
