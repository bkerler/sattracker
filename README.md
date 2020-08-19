# Sattracker
Standalone Satellite Tracker Project
(c) B. Kerler 2020
Licensed under MIT

## What is this about:
Creating a device to display satellite names and positions by pointing at them at the sky, and for helping to point satellite dishes/antennas to the right spot using a standalone offline device, based on daily tle database data.


## Hardware
- Esp32 WRoom DevC
- BNO055 Sensor (LSM303D can be used instead as well, just uncomment in source and comment out the BNO055)
- Beitian BN-280 GPS
- SSD1306 OLED 128x64 Display

## Connections
- BNO055 GND=>GND VCC=>3V3 SCL=>G21 SDA=>G22 
- I2C    GND=>GND VCC=>3V3 SCL=>G25 SDA=>G26
- BN-280 GND=>GND VCC=>3V3 TX=>G16  RX=>G17
- 5V Power via battery pack connected using ESP32 microusb or external power regulator

## Installation

- I recommend flashing my freshly built esp32-idf3 firmware, or use 1.12 instead

- Flash using:

``bash
esptool/esptool.py --chip esp32 --port /dev/ttyUSB0 --baud 460800 write_flash -z 0x1000 bin/firmware.bin
``

- Use thonny or any other ide to upload the files (frozen mpy for better ram usage and stability)

``
  src/main.py
  bin/bno055.mpy
  bin/micropyGPS.mpy
  bin/plan13.mpy
  bin/ssd1306.mpy
  src/geo.txt
``

or for development :

``
  src/main.py
  src/bno055.py
  src/micropyGPS.py
  src/plan13.py
  src/ssd1306.py
  src/geo.txt
``

You can replace the geo.txt with your own daily tle file downloaded from:
https://www.celestrak.com/NORAD/elements/

## Issues
Due to RAM restrictions on the ESP32, make sure not to load more than 100 satellites at a time (or restrict them by sorting out the elevation of the satellite to be higher than 30.0 as I did in main.py).

## Thanks go to
- Plan13 in Basic by G6LVB (which I ported and modified to micropython)
See http://www.g6lvb.com/Articles/LVBTracker2/index.htm
- Plan13 modifications in C++ by Mark VandeWettering K6HX and dl9sec
- Adafruit for the ssd1306 and lsm303d Library which I modified
- MicropyGPS by Michael Calvin McCoy (used as is)

## ToDo for the community
Enjoy and please contribute your improvements/corrections or post your own hardware mods based on it via git issues/commits !
