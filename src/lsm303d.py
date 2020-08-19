# Code based on Adafruit
# Heavily modified by B. Kerler 2020

import struct
from math import *
PI=3.14159

class LSM303D(object):
    ADDRESS_ACCEL = (0x32 >> 1)  # 0011001x
    ADDRESS_MAG   = (0x3C >> 1)  # 0011110x
    REGISTER_ACCEL_CTRL_REG1_A = 0x20 # 00000111   rw
    REGISTER_ACCEL_CTRL_REG4_A = 0x23 # 00000000   rw
    REGISTER_ACCEL_OUT_X_L_A   = 0x28
    REGISTER_MAG_CRB_REG_M     = 0x01
    REGISTER_MAG_MR_REG_M      = 0x02
    REGISTER_MAG_OUT_X_H_M     = 0x03

    MAGGAIN_1_3 = 0x20 # +/- 1.3
    MAGGAIN_1_9 = 0x40 # +/- 1.9
    MAGGAIN_2_5 = 0x60 # +/- 2.5
    MAGGAIN_4_0 = 0x80 # +/- 4.0
    MAGGAIN_4_7 = 0xA0 # +/- 4.7
    MAGGAIN_5_6 = 0xC0 # +/- 5.6
    MAGGAIN_8_1 = 0xE0 # +/- 8.1

    MAGRATE_0_7 = const(0x00)  # 0.75 Hz
    MAGRATE_1_5 = const(0x01)  # 1.5 Hz
    MAGRATE_3_0 = const(0x02)  # 3.0 Hz
    MAGRATE_7_5 = const(0x03)  # 7.5 Hz
    MAGRATE_15 = const(0x04)  # 15 Hz
    MAGRATE_30 = const(0x05)  # 30 Hz
    MAGRATE_75 = const(0x06)  # 75 Hz
    MAGRATE_220 = const(0x07)  # 220 Hz

    _LSM303ACCEL_MG_LSB        = 16704.0
    _GRAVITY_STANDARD          = 9.80665      # Earth's gravity in m/s^2
    _GAUSS_TO_MICROTESLA       = 100.0        # Gauss to micro-Tesla multiplier

    def __init__(self, i2c, hires=True):
        self.i2c = i2c
        self._data = bytearray(6)
        i2c.writeto_mem(self.ADDRESS_ACCEL, self.REGISTER_ACCEL_CTRL_REG1_A, bytearray([0x27]))
        i2c.writeto_mem(self.ADDRESS_ACCEL, self.REGISTER_ACCEL_CTRL_REG4_A, bytearray([0b00001000 if hires else 0]))
        i2c.writeto_mem(self.ADDRESS_MAG, self.REGISTER_MAG_MR_REG_M, bytearray([0x00]))
        self._lsm303mag_gauss_lsb_xy = 1100.0
        self._lsm303mag_gauss_lsb_z = 980.0
        self._mag_gain = self.MAGGAIN_1_3
        self._mag_rate = self.MAGRATE_0_7

    @property
    def _raw_magnetic(self):
        """The raw magnetometer sensor values.
        A 3-tuple of X, Y, Z axis values that are 16-bit signed integers.
        """
        self.i2c.readfrom_mem_into(self.ADDRESS_MAG, self.REGISTER_MAG_OUT_X_H_M, self._data)
        raw_values = struct.unpack_from(">hhh", self._data[0:6])
        return (raw_values[0], raw_values[2], raw_values[1])

    @property
    def magnetic(self):
        """The processed magnetometer sensor values.
        A 3-tuple of X, Y, Z axis values in microteslas that are signed floats.
        """
        mag_x, mag_y, mag_z = self._raw_magnetic
        return (
            mag_x / self._lsm303mag_gauss_lsb_xy * self._GAUSS_TO_MICROTESLA,
            mag_y / self._lsm303mag_gauss_lsb_xy * self._GAUSS_TO_MICROTESLA,
            mag_z / self._lsm303mag_gauss_lsb_z * self._GAUSS_TO_MICROTESLA,
        )
    
    @property
    def mag_gain(self):
        """The magnetometer's gain."""
        return self._mag_gain
    
    
    def _read_u8(self, device, address):
        self.i2c.readfrom_mem_into(self.ADDRESS_MAG, address, self._data)
        return self._data[0]

    def _write_u8(self, device, address, val):
        with device as i2c:
            i2c.writeto_mem(self.ADDRESS_MAG, address&0xFF, bytearray([val&0xFF]))
    
    @property        
    def acceleration(self):
        # Read the accelerometer as signed 16-bit little endian values.
        self.i2c.readfrom_mem_into(self.ADDRESS_ACCEL, self.REGISTER_ACCEL_OUT_X_L_A | 0x80, self._data)
        accel = struct.unpack('<hhh', self._data)
        # Convert to 12-bit values by shifting unused bits.
        accel = (accel[0] >> 4, accel[1] >> 4, accel[2] >> 4)
        x,y,z=accel
        m=[0,0,0]
        for i in range(0,3):
            m[i]=copysign(min(fabs(accel[i]),1.0),accel[i])
        pitch = asin(-1*m[0]);
        pitch_print = pitch * 180 / PI
        roll = asin(m[1] / cos(pitch)) if abs(cos(pitch))>=abs(m[1]) else 0
        roll_print = roll * 180 / PI;
        return x,y,z,pitch,roll
        
    @property
    def pitch(self):
        acc_x,acc_y,acc_z,pitch,roll=self.acceleration
        tilt=atan2(fabs(acc_z),acc_x)*180/PI
        if tilt>90 and tilt<180:
            tilt=180-tilt
        return tilt
    
    def tilt_comp(self,mag_x,mag_y,mag_z,pitch,roll):
        x=(mag_x*cos(pitch))+(mag_z*sin(pitch))
        y=(mag_x*sin(roll)*sin(pitch))+(mag_y*cos(roll))-(mag_z*sin(roll)*cos(pitch))
        z=(mag_x*cos(roll)*sin(pitch))+(mag_y*sin(roll))+(mag_z*cos(roll)*cos(pitch))
        return x,y,z
    
    @property
    def heading(self):
        mag_x,mag_y,mag_z=self.magnetic
        acc_x,acc_y,acc_z,pitch,roll=self.acceleration
        mag_x,mag_y,mag_z=self.tilt_comp(mag_x,mag_y,mag_z,pitch,roll)
        heading=(atan2(mag_y,mag_x)*180)/PI
        tilt=(atan2(fabs(acc_z),acc_x)*180/PI)
        if heading<0:
            heading=360+heading
        if 90-tilt>0:
            heading=360-tilt
        return heading

    @mag_gain.setter
    def mag_gain(self, gain=0x20):
        """Set the magnetometer gain.  Gain should be one of the following
        constants:
         - LSM303_MAGGAIN_1_3 = +/- 1.3 (default)
         - LSM303_MAGGAIN_1_9 = +/- 1.9
         - LSM303_MAGGAIN_2_5 = +/- 2.5
         - LSM303_MAGGAIN_4_0 = +/- 4.0
         - LSM303_MAGGAIN_4_7 = +/- 4.7
         - LSM303_MAGGAIN_5_6 = +/- 5.6
         - LSM303_MAGGAIN_8_1 = +/- 8.1
        """
        self._mag_gain = gain
        self.i2c.writeto_mem(self.ADDRESS_MAG, self.REGISTER_MAG_CRB_REG_M, self._mag_gain)
        if self._mag_gain == MAGGAIN_1_3:
            self._lsm303mag_gauss_lsb_xy = 1100.0
            self._lsm303mag_gauss_lsb_z = 980.0
        elif self._mag_gain == MAGGAIN_1_9:
            self._lsm303mag_gauss_lsb_xy = 855.0
            self._lsm303mag_gauss_lsb_z = 760.0
        elif self._mag_gain == MAGGAIN_2_5:
            self._lsm303mag_gauss_lsb_xy = 670.0
            self._lsm303mag_gauss_lsb_z = 600.0
        elif self._mag_gain == MAGGAIN_4_0:
            self._lsm303mag_gauss_lsb_xy = 450.0
            self._lsm303mag_gauss_lsb_z = 400.0
        elif self._mag_gain == MAGGAIN_4_7:
            self._lsm303mag_gauss_lsb_xy = 400.0
            self._lsm303mag_gauss_lsb_z = 355.0
        elif self._mag_gain == MAGGAIN_5_6:
            self._lsm303mag_gauss_lsb_xy = 330.0
            self._lsm303mag_gauss_lsb_z = 295.0
        elif self._mag_gain == MAGGAIN_8_1:
            self._lsm303mag_gauss_lsb_xy = 230.0
            self._lsm303mag_gauss_lsb_z = 205.0
            

    @property
    def mag_rate(self):
        """The magnetometer update rate."""
        return self._mag_rate

    @mag_rate.setter
    def mag_rate(self, value):
        assert value in (
            MAGRATE_0_7,
            MAGRATE_1_5,
            MAGRATE_3_0,
            MAGRATE_7_5,
            MAGRATE_15,
            MAGRATE_30,
            MAGRATE_75,
            MAGRATE_220,
        )

        self._mag_rate = value
        reg_m = ((value & 0x07) << 2) & 0xFF
        self.i2c.writeto_mem(self.ADDRESS_MAG, self.REGISTER_MAG_CRA_REG_M, self._mag_gain)        
