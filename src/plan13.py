#Algorithm based on Plan13 by G6LVB in Python 3.x
#See http://www.g6lvb.com/Articles/LVBTracker2/index.htm
#Some parts from Mark VandeWettering K6HX and dl9sec
#(c) Bjoern Kerler DO2BJK 2020

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


from math import *

#Change max x and max y if you use a world map
MAP_MAXX=1150
MAP_MAXY=609

RE = 6378.137              # WGS-84 Earth ellipsoid
FL = 1.0 / 298.257224      # -"-
RP = RE * (1.0 - FL)       # -

GM = 3.986E5               # Earth's Gravitational constant km^3/s^2
J2 = 1.08263E-3            # 2nd Zonal coeff, Earth's Gravity Field

YM = 365.25                # Mean Year,     days
YT = 365.2421874           # Tropical year, days
WW = 2.0 * pi / YT         # Earth's rotation rate, rads/whole day
WE = 2.0 * pi + WW         # Earth's rotation rate, radians/day
W0 = WE / 86400.0          # Earth's rotation rate, radians/sec

#Sidereal and Solar data. Rarely needs changing. Valid to year ~2030
#https://www.ngdc.noaa.gov/geomag/WMM/soft.shtml World Magnetic Model
YG = 2014.0                # GHAA, Year YG, Jan 0.0
G0 = 99.5828               # -"-
MAS0 = 356.4105            # MA Sun and rate, deg, deg/day MEAN_ANOMALY
MASD = 0.98560028          # -"-
INS = radians(23.4375)     # Sun's inclination
CNS = cos(INS)             # -"-
SNS = sin(INS)             # -"-
EQC1 = 0.03340             # Sun's Equation of centre terms
EQC2 = 0.00035             # -"-
AU = 149.597870700E6       # 1 AU, mean range in km to the sun (astronomical unit)


#Convert date to day-number
def fnday(y,m,d):
    if m<3:
        m+=12
        y-=1
    a=int(y*YM)
    b=int((m+1)*30.6)
    c=int(d-428)
    return a+b+c

#Convert day-number to date; valid 1900 Mar 01 - 2100 Feb 28
def fndate(dt):
    dt+=428
    y=int((dt-122.1)/YM)
    dt-=int(y*YM)
    m=int(dt/30.61)
    dt-=int(m*30.6)
    m-=1
    if m>12:
        m-=12
        y+=1
    d=dt
    return y,m,d

#Converts latitude (Breitengrad) -90..90° / longitude (Laengengrad) -180..180°
#to x/y-coordinates of a map with maxamimum dimension MapMaxX * MapMaxY
def latlon2xy(lat,lon,MapMaxX,MapMaxY):
    x=int(((180.0+lon)/360.0)*MapMaxX)
    y=int(((90.0-lat)/180.0)*MapMaxY)
    return x,y

class P13DateTime():
    DN=0
    TN=0.0
    
    def __init__(self,iYear,iMonth,iDay,iHour,iMinute,iSecond):
        self.settime(iYear,iMonth,iDay,iHour,iMinute,iSecond)

    def settime(self,year,month,day,h,m,s):
        self.DN=fnday(year,month,day)
        self.TN=(h + m/60.0 + s/3600.0) / 24.0

    def gettime(self):
        year, month, day=fndate(self.DN)
        t = self.TN
        t *= 24.0
        h = int(t)
        t -= h
        t *= 60.0
        m = int(t)
        t -= m
        t *= 60.0
        s = int(t)
        return int(year),int(month),int(day),h,m,s

    def ascii(self):
        year,month,day,h,m,s=self.gettime()
        return "%4d-%02d-%02d %02d:%02d:%02d" % (year, month, day, h, m, s)

    def roundup(self,t):
        inc = t - fmod(self.TN, t)
        self.TN += inc
        self.DN += self.TN
        self.TN -= self.TN

class P13Observer():
    U=[0,0,0]
    E=[0,0,0]
    N=[0,0,0]
    O=[0,0,0]
    V=[0,0,0]
    
    def __init__(self,lat,lon,asl):
        LA=radians(lat)
        LO=radians(lon)
        HT=asl/1000.0
        CL=cos(LA)
        SL=sin(LA)
        CO=cos(LO)
        SO=sin(LO)
        self.U[0] = CL * CO
        self.U[1] = CL * SO
        self.U[2] = SL
        self.E[0] = -SO
        self.E[1] =  CO
        self.E[2] =  0.0
        self.N[0] = -SL * CO
        self.N[1] = -SL * SO
        self.N[2] =  CL
        D = sqrt(RE * RE * CL * CL + RP * RP * SL * SL)
        Rx = (RE * RE) / D + HT
        Rz = (RP * RP) / D + HT
        self.O[0] = Rx * self.U[0]
        self.O[1] = Rx * self.U[1]
        self.O[2] = Rz * self.U[2]
        self.V[0] = -self.O[1] * W0
        self.V[1] =  self.O[0] * W0
        self.V[2] =  0.0
 
class tle_fields():
    def __init__(self,tleline="",tlefile=""):
        self.tleline=tleline
        self.tlefile=tlefile

    def nextsatellite(self):
        if self.tlefile!="":
            for tledata in self.parse_tle_file(self.tlefile):
                yield self.parse(tledata)
        elif self.tleline!="":
            yield self.parse(self.tleline)

    def parse(self,ln):
        name=ln[0]
        l1=ln[1]
        l2=ln[2]
        return name,l1,l2

    def parse_tle_file(self,tlefile):
        with open(tlefile, "r") as rf:
            i = 0
            data = []
            for line in rf:
                data.append(line)
                i += 1
                if i == 3:
                    yield data
                    data = []
                    i=0

class tle():
    name=""
    N=0.0
    YE=0
    M2=0.0
    IN=0.0
    RA=0.0
    EC=0.0
    WP=0.0
    MA=0.0
    MM=0.0
    RV=0.0
    DE=0.0
    TE=0.0
    def __init__(self,tledata):
        self.name=tledata[0].rstrip()
        self.N = int(tledata[1][2:7])  # Satellite catalog number, L1 02-06
        self.YE = int(tledata[1][18:20])  # Epoch year year 11 18-19
        if self.YE < 58:
            self.YE += 2000
        else:
            self.YE += 1900
        self.TE = float(tledata[1][20:32])               # Epoch time days L1
        self.M2 = 2.0 * pi * float(tledata[1][33:43])    # Decay Rare rev/d/d L1
        self.IN = radians(float(tledata[2][8:16]))       # Inclination deg L2
        self.RA = radians(float(tledata[2][17:26]))      # R.A.A.N deg L2
        self.EC = float(tledata[2][26:33]) / 1.0E7       # Eccentricity L2
        self.WP = radians(float(tledata[2][34:42]))      # Arg of perigee deg L2
        self.MA = radians(float(tledata[2][43:51]))      # Mean anomaly deg L2 43-50
        self.MM = 2.0 * pi * float(tledata[2][52:63])    # Mean motion rev/d L2 52-62
        self.RV = int(tledata[2][63:68])                 # Orbit number
        self.DE = int(fnday(self.YE, 1, 0) + (self.TE))
        self.TE -= int(self.TE)

class P13Satellite():
    SAT=[0,0,0]
    VEL=[0,0,0]
    S=[0,0,0]
    V=[0,0,0]

    '''
    ALON=0.0   #Sat attitude deg
    ALAT=0.0   #Sat attitude deg
    DE=0       #Epoch Fraction of day
    
    #These values are stored, but could be calculated on the fly during calls to predict() 
    #Classic space/time tradeoff
    
    N0=0.0
    A_0=0.0
    PC=0.0
    QD=0.0
    WD=0.0
    DC=0.0
    
    RS=0.0     #Radius of satellite orbit
    RR=0.0     #Range rate for doppler calculation
    '''

    def __init__(self,tledata):
        self.tle=tle(tledata)
        self.name=self.tle.name
        self.N0=self.tle.MM / 86400.0
        self.A_0 = pow(GM / (self.N0 * self.N0), 1.0/3.0)
        self.B_0 = self.A_0 * sqrt(1.0 - self.tle.EC * self.tle.EC)
        self.PC  = RE * self.A_0 / (self.B_0 * self.B_0)
        self.PC  = 1.5 * J2 * self.PC * self.PC * self.tle.MM
        self.CI = cos(self.tle.IN)
        self.QD = -self.PC * self.CI
        self.WD =  self.PC * (5.0 * self.CI * self.CI - 1.0) / 2.0
        self.DC = -2.0 * self.tle.M2 / (3.0 * self.tle.MM)

    def predict(self,dt):
        '''
        long   DN;
        double TN;
        double GHAE, GHAA;
        double T, DT, KD, KDP;
        double M, DR, RN, EA;
        double DNOM, C_EA, S_EA;
        double A, B, D;
        double AP, CW, SW;
        double RAAN;
        double CQ, SQ;
        double CI, SI;
        double CG, SG;
        '''
        
        CX=[0,0,0]
        CY=[0,0,0]
        CZ=[0,0,0]
    
        DN = dt.DN
        TN = dt.TN

        TEG=self.tle.DE - fnday(YG, 1, 0) + self.tle.TE
        GHAE = radians(G0) + TEG * WE   #GHA Aries, epoch

        T   = (DN - self.tle.DE) + (TN - self.tle.TE)    # Elapsed T since epoch, days
        DT  = self.DC * T / 2.0                          # Linear drag terms
        KD  = 1.0 + 4.0 * DT                             # -"-
        KDP = 1.0 - 7.0 * DT                             # -"-
  
        M   = self.tle.MA + self.tle.MM * T * (1.0 - 3.0 * DT)   # Mean anomaly at YR,TN
        DR  = int(M / (2.0 * pi))                       # Strip out whole no of revs
        M  -= DR * 2.0 * pi                             # M now in range 0..2PI
        self.RN  = self.tle.RV + DR                          # Current Orbit number
    
        #Solve M = EA - EC*SIN(EA) for EA given M, by Newton's Method
        EA  = M                                #Initial solution

        D=2.0E-5
        while (fabs(D) > 1.0E-5):
            C_EA = cos(EA)
            S_EA = sin(EA)
            DNOM = 1.0 - self.tle.EC * C_EA
            D = (EA - self.tle.EC * S_EA - M) / DNOM     # Change to EA for better solution
            EA -= D                             # by this amount

        # Distances
        A = self.A_0 * KD
        B = self.B_0 * KD
        self.RS = A * DNOM

        # Calc satellite position & velocity in plane of ellipse
        self.S[0] = A * (C_EA - self.tle.EC)
        self.S[1] = B * S_EA
    
        self.V[0] = -A * S_EA / DNOM * self.N0
        self.V[1] =  B * C_EA / DNOM * self.N0

        AP = self.tle.WP + self.WD * T * KDP
        CW = cos(AP)
        SW = sin(AP)
        RAAN = self.tle.RA + self.QD * T * KDP
        CQ = cos(RAAN)
        SQ = sin(RAAN)

        # CX, CY, and CZ form a 3x3 matrix that converts between orbit
        # coordinates, and celestial coordinates.
    
        # Plane -> celestial coordinate transformation, [C] = [RAAN]*[IN]*[AP]
        CI = cos(self.tle.IN)
        SI = sin(self.tle.IN)
       
        CX[0] =  CW * CQ - SW * CI * SQ
        CX[1] = -SW * CQ - CW * CI * SQ
        CX[2] =  SI * SQ

        CY[0] =  CW * SQ + SW * CI * CQ
        CY[1] = -SW * SQ + CW * CI * CQ
        CY[2] = -SI * CQ

        CZ[0] = SW * SI
        CZ[1] = CW * SI
        CZ[2] = CI

        # Compute SATellite's position vector and VELocity in
        # CELESTIAL coordinates. (Note: Sz=S[2]=0, Vz=V[2]=0)
        self.SAT[0] = self.S[0] * CX[0] + self.S[1] * CX[1]
        self.SAT[1] = self.S[0] * CY[0] + self.S[1] * CY[1]
        self.SAT[2] = self.S[0] * CZ[0] + self.S[1] * CZ[1]

        self.VEL[0] = self.V[0] * CX[0] + self.V[1] * CX[1]
        self.VEL[1] = self.V[0] * CY[0] + self.V[1] * CY[1]
        self.VEL[2] = self.V[0] * CZ[0] + self.V[1] * CZ[1]

        # Also express SAT and VEL in GEOCENTRIC coordinates:
        GHAA = (GHAE + WE * T)  # GHA Aries at elapsed time T
        CG   = cos(-GHAA)
        SG   = sin(-GHAA)

        self.S[0] = self.SAT[0] * CG - self.SAT[1] * SG
        self.S[1] = self.SAT[0] * SG + self.SAT[1] * CG
        self.S[2] = self.SAT[2]

        self.V[0] = self.VEL[0] * CG - self.VEL[1]* SG
        self.V[1] = self.VEL[0] * SG + self.VEL[1]* CG
        self.V[2] = self.VEL[2]
    
    def latlon(self):
        lat=degrees(asin(self.S[2]/self.RS))
        lon=degrees(atan2(self.S[1],self.S[0]))
        return lat,lon
    
    def elaz(self,obs):
        R=[0,0,0]
        
        # Rangevec = Satvec - Obsvec
        R[0] = self.S[0] - obs.O[0]
        R[1] = self.S[1] - obs.O[1]
        R[2] = self.S[2] - obs.O[2]
        
        # Range magnitude
        r = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2])
        
        # Normalise Range vector
        R[0] /= r
        R[1] /= r
        R[2] /= r
        
        # UP Component of unit range
        u = R[0] * obs.U[0] + R[1] * obs.U[1] + R[2] * obs.U[2]
        # EAST
        e = R[0] * obs.E[0] + R[1] * obs.E[1]
        # NORTH
        n = R[0] * obs.N[0] + R[1] * obs.N[1] + R[2] * obs.N[2]
        
        # Azimuth
        az = degrees(atan2(e, n))
        
        if (az < 0.0):
            az += 360.0
        
        # Elevation
        el = degrees(asin(u))
        
        # Calculate range rate needed for doppler calculation
        # Resolve Sat-Obs velocity vector along unit range vector. (VOz=obs.V[2]=0)
        self.RR  = (self.V[0] - obs.V[0]) * R[0] + (self.V[1] - obs.V[1]) * R[1] + self.V[2] * R[2]   # Range rate, km/s
        return el, az
    
    def footprint(self,numberofpoints,MapMaxX,MapMaxY,satlat,satlon):
        '''
        int i;
        
        double srad;
        double cla, sla, clo, slo;
        double sra, cra;
        
        double a, x, y, z, Xfp, Yfp, Zfp;
        '''
            
        points=[]
        srad = acos(RE / self.RS)   # Radius of footprint circle
        sra  = sin(srad)       # Sin/Cos these to save time
        cra  = cos(srad)

        cla  = cos(radians(satlat))
        sla  = sin(radians(satlat))
        clo  = cos(radians(satlon))
        slo  = sin(radians(satlon))
        
        for i in range(0,numberofpoints):  # "numberofpoints" points to the circle
            a = 2.0 * pi * i / numberofpoints   # Angle around the circle
            Xfp = cra                                          # Circle of points centred on Lat=0, Lon=0
            Yfp = sra * sin(a)                                 # assuming Earth's radius = 1
            Zfp = sra * cos(a)
            
            x = Xfp * cla - Zfp * sla                        # Rotate point "up" by latitude "satlat"
            y = Yfp                                          # -"-
            z = Xfp * sla + Zfp * cla                        # -"-
            
            Xfp = x * clo - y * slo                          # Rotate point "around" through longitude "satlon"
            Yfp = x * slo + y * clo                          # -"-
            Zfp = z                                          # -"-
            
            # Convert point to Lat/Lon and convert/scale to a pixel map
            sx,sy=latlon2xy(degrees(asin(Zfp)), degrees(atan2(Yfp,Xfp)), MapMaxX, MapMaxY)
            points.append([sx,sy])
        return points

    # Returns the RX (dir = 0 or P13_FRX) or TX (dir = 1 or P13_FTX) frequency with doppler shift.
    def doppler(self,freqMHz,direction):
        dopplershift = -freqMHz * self.RR / 299792.0     #Speed of light is 299792.0 km/s
        if direction: #TX
            freqMHz = freqMHz - dopplershift
        else:
            freqMHz = freqMHz + dopplershift
        return freqMHz
    
    
class P13Sun():
    SUN=[0,0,0]
    H=[0,0,0]
    
    def predict(self,dt):
        '''
        long   DN;
        double TN;
        double T, GHAE, MRSE, MASE, TAS;
        double C, S;
        '''
        DN = dt.DN
        TN = dt.TN

        T    = (DN - fnday(YG, 1, 0)) + TN       #Elapsed Time: Epoch - YG
        GHAE = radians(G0) + T * WE              # GHA Aries, epoch
            
        MRSE = radians(G0) + T * WW + pi         #Mean RA Sun at Sat epoch
        MASE = radians(MAS0 + T * MASD)          #Mean anomaly of the Sun
        TAS  = MRSE + EQC1 * sin(MASE) + EQC2 * sin(2.0 * MASE)

        # Sin/Cos Sun's true anomaly
        C = cos(TAS)
        S = sin(TAS)

        # Sun unit vector - CELESTIAL coords
        self.SUN[0] = C
        self.SUN[1] = S * CNS
        self.SUN[2] = S * SNS
        
        # Obtain SUN unit vector in GEOCENTRIC coordinates
        C = cos(-GHAE) 
        S = sin(-GHAE)
            
        self.H[0] = self.SUN[0] * C - self.SUN[1] * S
        self.H[1] = self.SUN[0] * S + self.SUN[1] * C
        self.H[2] = self.SUN[2]

    def latlon(self):
        lat=degrees(asin(self.H[2]))
        lon=degrees(atan2(self.H[1],self.H[0]))
        return lat,lon
    
    def elaz(self,obs):
        el=0.0
        az=0.0
        # ToDo:
        # Convert the celestial coordinates SUN[] to
        # azimuth and elevation of an observer.
        return el,az
    
    # Generates the sunlight footprint at satlat/satlon and calculates rectangular
    # x/y coordinates scaled to a map with size MapMaxX/MapMaxY. The coordinates are stored
    # in a two dimensional array. points[n][0] stores x and points[n][1] stores y.
    # The coordinates can be concatenated with lines to create a footprint outline.
    # This is a simplified aproach with no real calculation of the distance to the sun at a
    # specific time. It is assumed that the nearest and farest distance of the sun makes almost no
    # difference in footprint radius, it is always almost 0.5*PI. Therefore one astronomical
    # unit is used for the distance. The same algorithm is used as for the satellite footprint
    # except, that RS is replaced by AU.

    def footprint(self,numberofpoints,MapMaxX,MapMaxY,sunlat,sunlon):
        '''
        int i;
        double srad;
        double cla, sla, clo, slo;
        double sra, cra;
        double a, x, y, z, Xfp, Yfp, Zfp;
        '''
        points=[]    
        srad = acos(RE / AU)   # Radius of sunlight footprint circle
        sra  = sin(srad)       # Sin/Cos these to save time
        cra  = cos(srad)

        cla  = cos(radians(sunlat))
        sla  = sin(radians(sunlat))
        clo  = cos(radians(sunlon))
        slo  = sin(radians(sunlon))
        
        for i in range(0,numberofpoints):                      # "numberofpoints" points to the circle
            a = 2.0 * pi * i / numberofpoints                  # Angle around the circle
            Xfp = cra                                          # Circle of points centred on Lat=0, Lon=0
            Yfp = sra * sin(a)                                 # assuming Earth's radius = 1
            Zfp = sra * cos(a)
            
            x   = Xfp * cla - Zfp * sla                        # Rotate point "up" by latitude "sunlat"
            y   = Yfp                                          # -"-
            z   = Xfp * sla + Zfp * cla                        # -"-
            
            Xfp = x * clo - y * slo                            # Rotate point "around" through longitude "sunlon"
            Yfp = x * slo + y * clo                            # -"-
            Zfp = z                                            # -"-
            
            #Convert point to Lat/Lon and convert/scale to a pixel map
            sx, sy=latlon2xy(degrees(asin(Zfp)), degrees(atan2(Yfp,Xfp)), MapMaxX, MapMaxY)
            points.append([sx,sy])
        return points
            
class satpredict():
    def __init__(self,lat,lon,alt,MAP_MAXX,MAP_MAXY):
        self.curpos = P13Observer(lat, lon, alt)
        self.curlat=lat
        self.curlon=lon
        self.mapmaxx=MAP_MAXX
        self.mapmaxy=MAP_MAXY
        self.curtime=None
        self.sat=None

    def curpos_to_xy(self):
        return latlon2xy(self.curlat, self.curlon, self.mapmaxx, self.mapmaxy)

    def settime(self,year,month,day,hour,minute,second):
        self.curtime = P13DateTime(year, month, day, hour, minute, second)

    def sat_predict(self,tledata):
        if self.curtime!=None:
            self.sat = P13Satellite(tledata)
            self.sat.predict(self.curtime)
            self.satlat, self.satlon = self.sat.latlon()
            self.satel, self.sataz = self.sat.elaz(self.curpos)
            self.satx, self.saty=latlon2xy(self.satlat, self.satlon, self.mapmaxx, self.mapmaxy)
            return self.satlat,self.satlon,self.sataz,self.satel, self.satx, self.saty
        return None

    def getdoppler(self,freqRX=0.0,freqTX=0.0):
        P13_FRX = 0
        P13_FTX = 1
        if self.sat!=None:
            return self.sat.doppler(freqRX, P13_FRX), self.sat.doppler(freqTX, P13_FTX)
        return None

    def gettime(self):
        if self.curtime!=None:
            return self.curtime.ascii()

    def getsatfootprint(self):
        if self.sat!=None:
            return self.sat.footprint(32, MAP_MAXX, MAP_MAXY, self.satlat, self.satlon)
        return None

#Uncomment here for testing purposes on PC
'''
def main():
    MyLAT   =  44.444444   # Latitude (Breitengrad): N -> +, S -> -
    MyLON   =  10.666666   # Longitude (Längengrad): E -> +, W -> -
    MyALT   = 386.0        # Altitude ASL (m)
    freqRX  = 145.800      # Nominal downlink frequency
    freqTX  = 437.800      # Nominal uplink frequency

    Year    = 2020         # Set start year
    Month   = 8            # Set start month
    Day     = 17           # Set start day
    Hour    = 14           # Set start hour
    Minute  = 30            # Set start minute
    Second  = 00           # Set start second

    #https://celestrak.com/NORAD/elements/active.txt
    #tledata = "ISS (ZARYA)\n1 25544U 98067A   20229.52810887  .00000218  00000-0  12046-4 0  9992\n2 25544  51.6457  49.8366 0001418  43.4994 351.5363 15.49164518241370\n"
    #tledata = "ASTRA 1F    \n1 23842U 96021A   20231.14647696  .00000131  00000-0  00000+0 0  9992\n2 23842   0.0411 342.3378 0002489 172.6506 269.0303  1.00275321 88984\n"
    for tledata in tle_fields("","geo.txt").nextsatellite():
        satp=satpredict(MyLAT,MyLON,MyALT,MAP_MAXX,MAP_MAXY)
        ixQTH,iyQTH=satp.curpos_to_xy()
        satp.settime(Year,Month,Day,Hour,Minute,Second)
        satlat,satlon,sataz,satel, satx, saty=satp.sat_predict(tledata)
        curtime=satp.gettime()
        print("Prediction for %s (MAP %dx%d: x = %d,y = %d):" % (satp.sat.name, MAP_MAXX, MAP_MAXY, ixQTH, iyQTH))
        print("%s -> Lat: %.4f Lon: %.4f (MAP %dx%d: x = %d,y = %d) Az: %.2f El: %.2f" % (curtime, satlat, satlon, MAP_MAXX, MAP_MAXY, satx, saty, sataz, satel))
        print("RX: %.6f MHz, TX: %.6f MHz" % satp.getdoppler(dfreqRX,dfreqTX))
        print("\n")

    print("Satellite footprint map coordinates:")
    aiSatFP = satp.getsatfootprint()
    for i in range(0,32):
        print("%2d: x = %d, y = %d" % (i, aiSatFP[i][0], aiSatFP[i][1]))

    Sun=P13Sun()
    Sun.predict(satp.curtime)
    dSunLAT,dSunLON=Sun.latlon()
    dSunEL,dSunAZ=Sun.elaz(satp.curpos)
    ixSUN,iySUN=latlon2xy(dSunLAT, dSunLON, MAP_MAXX, MAP_MAXY)
    print("\nSun -> Lat: %.4f Lon: %.4f (MAP %dx%d: x = %d,y = %d) Az: %.2f El: %.2f" % (dSunLAT, dSunLON, MAP_MAXX, MAP_MAXY, ixSUN, iySUN, dSunAZ, dSunEL))
    print("Sunlight footprint map coordinates:")
 
    aiSunFP=Sun.footprint(32, MAP_MAXX, MAP_MAXY, dSunLAT, dSunLON)
    for i in range(0,32):
        print("%2d: x = %d, y = %d" % (i, aiSunFP[i][0], aiSunFP[i][1]))

if __name__=="__main__":
    main()
'''