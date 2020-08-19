# Satellite tracker, standalone ESP32
# (c) B. Kerler 2020

# The MIT License (MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


from machine import Pin, I2C
import ssd1306
from machine import Timer
from time import sleep
from math import atan2
from bno055 import *
from plan13 import *
import time
from micropyGPS import MicropyGPS
from machine import UART
import utime, gc, _thread

sensori2c  = I2C(-1, Pin(21), Pin(22),timeout=1000)
#from lsm303d import LSM303D
#lsm303 = LSM303D(sensori2c)
imu=BNO055(sensori2c)
calibrated=False

oledi2c = I2C(-1, scl=Pin(25), sda=Pin(26))
oled_width = 128
oled_height = 64
oled = ssd1306.SSD1306_I2C(oled_width, oled_height, oledi2c)

my_gps=MicropyGPS()
uart = UART(2, rx=16, tx=17, baudrate=9600)
my_gps = MicropyGPS(local_offset=0, location_formatting='ddm')

'''
def showheading_lsm303():
    #accel, mag = lsm303.read()
    mag_x,mag_y,mag_z=lsm303.magnetic
    raw_mag_x,raw_mag_y,raw_mag_z=lsm303._raw_magnetic
    heading=lsm303.heading
    pitch=lsm303.pitch
    
    #print("Acceleration (m/s^2): X=%0.3f Y=%0.3f Z=%0.3f"%accel)
    print("Magnetometer (micro-Teslas)): X=%0.3f Y=%0.3f Z=%0.3f"%(mag_x,mag_y,mag_z))
    print("Magnetometer (micro-Teslas)): X=%0.3f Y=%0.3f Z=%0.3f"%(raw_mag_x,raw_mag_y,raw_mag_z))
    print("Heading %0.3f"%heading)
    print("Pitch %0.3f"%pitch)
    return heading, pitch
'''

def showheading_bno066(debug):
    global calibrated
    if not calibrated:
        calibrated=imu.calibrated
        print('Calibration required: sys {} gyro {} accel {} mag {}'.format(*imu.cal_status()))
        return -1,-1
    else:
        if debug:
            print('Temperature {}°C'.format(imu.temperature()))
            print('Mag       x {:5.0f}    y {:5.0f}     z {:5.0f}'.format(*imu.mag()))
            print('Gyro      x {:5.0f}    y {:5.0f}     z {:5.0f}'.format(*imu.gyro()))
            print('Accel     x {:5.1f}    y {:5.1f}     z {:5.1f}'.format(*imu.accel()))
            print('Lin acc.  x {:5.1f}    y {:5.1f}     z {:5.1f}'.format(*imu.lin_acc()))
            print('Gravity   x {:5.1f}    y {:5.1f}     z {:5.1f}'.format(*imu.gravity()))
            print('Azimuth     {:4.0f} Roll {:4.0f} Elevation {:4.0f}'.format(heading,roll,pitch))
        heading,roll,pitch=imu.euler()
        pitch=90-pitch
        if pitch>90:
            pitch=180-pitch
            if heading>180:
                heading-=180
            else:
                heading+=180
        
        return heading, pitch


#timer = Timer(-1)
#timer.init(period=500, mode=Timer.PERIODIC, callback=lambda t:showheading())
#timer.stop()

def convertgpsposition():
    lon=my_gps.longitude
    lat=my_gps.latitude
    alt=my_gps.altitude
    if len(lat)>2:
       latitude = lat[0] + (lat[1]/60)
       if lat[2] == 'S':
           latitude = -latitude
    if len(lon)>2:
       longitude = lon[0] + (lon[1]/60)
       if lon[2] == 'W':
           longitude = -longitude

    hours,minutes,seconds=my_gps.timestamp
    seconds=int(seconds)
    day,month,year=my_gps.date
    year=2000+year
    return longitude,latitude,alt,day,month,year,hours,minutes,seconds

def getsatdata(tledata):
    longitude,latitude,alt,day,month,year,hours,minutes,seconds=convertgpsposition()
    
    freqRX  = 145.800      # Nominal downlink frequency
    freqTX  = 437.800      # Nominal uplink frequency

    #https://celestrak.com/NORAD/elements/active.txt
    satp=satpredict(latitude,longitude,alt,MAP_MAXX,MAP_MAXY)
    curx,cury=satp.curpos_to_xy()
    satp.settime(year,month,day,hours,minutes,seconds)
    satlat,satlon,sataz,satel, satx, saty=satp.sat_predict(tledata)
    curtime=satp.gettime()
    rx,tx=satp.getdoppler(freqRX,freqTX)
    return satp.sat.name, curtime, satlat, satlon, sataz, satel, curx, cury, satx, saty, rx,tx

def getgpsdata():
    n = 0
    mem_free = gc.mem_free()
    tm_last = 0
    rxlength = uart.any()
    if rxlength>0:
            b = uart.read(rxlength)
            for x in b:
                if 10 <= x <= 126:
                    try:
                        stat = my_gps.update(chr(x))
                    except:
                        stat=False
                    if stat:
                        tm = my_gps.timestamp
                        tm_now = (tm[0] * 3600) + (tm[1] * 60) + int(tm[2])
                        if (tm_now - tm_last) >= 10:
                            n += 1
                            tm_last = tm_now
                            if (n % 10) == 0:
                                print("Mem free:", gc.mem_free(), mem_free - gc.mem_free())
                                gc.collect()

#satellitedb=parse_tle_file("geo.txt")
tledata = ["ASTRA 1F",
           "1 23842U 96021A   20231.14647696  .00000131  00000-0  00000+0 0  9992",
           "2 23842   0.0411 342.3378 0002489 172.6506 269.0303  1.00275321 88984"
           ]
debug=False

while True:
    getgpsdata()
    longitude,latitude,alt,day,month,year,hours,minutes,seconds=convertgpsposition()
    if longitude==0.0 and latitude==0.0:
        oled.fill(0)
        oled.text("GPS isn't locked !", 0, 0)
        oled.text("Longitude: %02.2f"%longitude, 0, 10)
        oled.text("Latitude: %02.2f"%latitude, 0, 20)
        oled.text("SIU: %d"%my_gps.satellites_in_use,0,30)
        oled.text("SUD: %s"%",".join(str(my_gps.satellites_used)),0,30)
        oled.show()
    else:
        break
    time.sleep_us(500000)

satdb=[]
satm={}
entry=0

for tledata in tle_fields("","geo.txt").nextsatellite():
    satname, curtime, satlat, satlon, sataz, satel, curx, cury, satx, saty, rx,tx=getsatdata(tledata)
    if satel>30.0:
       if int(sataz) not in satm:
           satdb.append([tledata,sataz,satel])
           print("%s : Az: %02.2f El:%02.2f" % (satname,sataz,satel))
           satm[int(sataz)]=1
           entry+=1
           oled.fill(0)
           oled.text('P: %d visible sats !'%entry, 0, 0)
           oled.show()
       #if entry>100:
       #    break

while True:
    azimuth,elevation=showheading_bno066(debug)
    longitude=0.0
    if azimuth==-1 and elevation==-1:
        oled.fill(0)
        oled.text('Calibration required !', 0, 0)
        oled.show()
    else:
        oled.fill(0)
        oled.text('Azimuth: '+str(azimuth), 0, 0)
        oled.text('Elevation: '+str(elevation), 0, 10)
        factor=.5
        idx=0
        for sat in satdb:
            tledata=sat[0]
            sataz=sat[1]
            satel=sat[2]
            if sataz<azimuth+factor and sataz>azimuth-factor:
                satname, curtime, satlat, satlon, sataz, satel, curx, cury, satx, saty, rx,tx=getsatdata(tledata)
                print("Prediction for %s (MAP %dx%d: x = %d,y = %d):\n" % (satname, MAP_MAXX, MAP_MAXY, curx, cury))
                print("%s -> Lat: %.4f Lon: %.4f (MAP %dx%d: x = %d,y = %d) Az: %.2f El: %.2f" % (curtime, satlat, satlon, MAP_MAXX, MAP_MAXY, satx, saty, sataz, satel))
                print("RX: %.6f MHz, TX: %.6f MHz" % (rx,tx))
                oled.text('S: '+satname,0,20)
                #oled.text('SLon: %02f SLat: %02f '%(satlon,satlat),0,30)
                oled.text('SAz: %02.2f'%sataz,0,30)
                oled.text('SEl: %02.2f'%satel,0,40)
            idx+=1
        oled.show()
        
    time.sleep_us(500000)