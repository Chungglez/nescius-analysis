from numpy import pi
p_atm = 14.7    # psia - atmospheric pressue at sea level

# Unit conversions functions
#-----------------------------------------------------------------------------------------------------------------------
mm2m = lambda length: length*.001
in2m = lambda length: length*0.0254
microm2m = lambda length: length*.000001
m2mm = lambda length: length*1000
ci2m3 = lambda volume: volume * .0000163871
psia2Pa = lambda psia: psia * 6894.76
psig2Pa = lambda psig: psia2Pa(psig + p_atm)
Pa2psig = lambda Pa: Pa/6894.76 + p_atm
Pa2psia = lambda Pa: Pa/6894.76
F2K = lambda temp: (temp - 32)*5/9 + 273.15
K2C = lambda temp: (temp-273.15)
K2R = lambda temp: temp*1.8
g2kg = lambda mass: mass/1000
kg2g = lambda mass: mass*1000
kg2lb = lambda mass: mass*2.20462
lb2N = lambda force: force*4.44822
N2lb = lambda force: force/4.44822
m2in = lambda length: length*39.3701
m22in2 = lambda area: area*1550
m32in3 = lambda volume: volume*61023.7
Pa2MPa = lambda stress: stress/1000000


# Convert from SI units to Hutzel gross Imperial units
#-----------------------------------------------------------------------------------------------------------------------
Cp2Imp = lambda Cp: Cp/4186.82              # J/kg-K       ▶ Btu/lb-F
visc2Imp = lambda visc: visc*0.055997410    # Pa-s         ▶ lb/in-s
h2SI = lambda h: h*2943623                  # Btu/in^2-s-F ▶ J/m^2-s-K

# Geometric conversions
#-----------------------------------------------------------------------------------------------------------------------
area2d = lambda area: (4*area/pi)**0.5
area2r = lambda area: (area/pi)**0.5
r2area = lambda radius: pi*radius**2
d2area = lambda diameter: pi/4*diameter**2
deg2rad = lambda deg: deg * pi/180