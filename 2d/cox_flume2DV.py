import numpy as np
from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus import WaveTools as wt


opts= Context.Options([
    ('ns_model',1,"ns_model={0,1} for {rans2p,rans3p}"),
    ("final_time",178.0,"Final time for simulation"),
    ("sampleRate",0.1,"Time interval to output solution"),
    ("gauges", True, "Collect data for validation"),
    ("cfl",0.9,"Desired CFL restriction"),
    ("he",0.04,"Max mesh element diameter"),
    ("hv",1.03,"water depth at vegetation"),
    ("Hm",0.20,"significant wave height"),
    ("Tp",1.85,"peak wave period"),
    ("mangrove_porous",False,"turn on porosity for mangroves"),
    ("filename",'wg1_regular_baseline_trial2.csv',"turn on porosity for mangroves"),
    ("wave_type",'Time',"runs simulation with random waves")
    ])


# general options
# sim time
T = opts.final_time
# initial step
dt_init = 0.001
# CFL value
cfl = 0.5
# mesh size
he = opts.he
# rate at which values are recorded
sampleRate = opts.sampleRate

# physical options
# water density
rho_0 = 998.2
# water kinematic viscosity
nu_0 = 1.004e-6
# air density
rho_1 = 1.205
# air kinematic viscosity
nu_1 = 1.5e-5
# gravitational acceleration
g = np.array([0., -9.81, 0.])



# wave options
water_level = opts.hv+0.85
wave_period = opts.Tp
wave_height = opts.Hm
wave_direction = np.array([1., 0., 0.])

#Monochromatic or Random
if opts.wave_type=='Monochromatic':
    wave = wt.MonochromaticWaves(period=wave_period, waveHeight=wave_height, mwl=water_level, depth=water_level, g=g, waveDir=wave_direction,waveType='Fenton',Nf=8)

elif opts.wave_type=='Time':
#    wave = wt.TimeSeries(timeSeriesFile="wg1_regular_baseline_trial2.csv", # e.g.= "Timeseries.txt",
    wave = wt.TimeSeries(timeSeriesFile=opts.filename, # e.g.= "Timeseries.txt",
                         skiprows=0,
                         timeSeriesPosition=np.array([0.,0.,0.]),
                         depth=water_level,
                         N=5000,          #number of frequency bins
                         mwl=water_level,        #mean water level
                         waveDir=wave_direction, 
                         g=g,
                         cutoffTotal = 0.001,
                         rec_direct = True,
                         window_params = None)

else:
    wave = wt.RandomWaves(Tp=wave_period,Hs=wave_height,mwl=water_level,depth=water_level,g=g,waveDir=wave_direction,spectName='JONSWAP',N=300,bandFactor=2.5)

wavelength = 10.

#  ____                        _
# |  _ \  ___  _ __ ___   __ _(_)_ __
# | | | |/ _ \| '_ ` _ \ / _` | | '_ \
# | |_| | (_) | | | | | | (_| | | | | |
# |____/ \___/|_| |_| |_|\__,_|_|_| |_|
# Domain
# All geometrical options go here (but not mesh options)

domain = Domain.PlanarStraightLineGraphDomain()

# ----- SHAPES ----- #

# ----- TANK ----- #

boundaryOrientations = {'y-': np.array([0., -1.,0.]),
                        'x+': np.array([+1., 0.,0.]),
                        'y+': np.array([0., +1.,0.]),
                        'x-': np.array([-1., 0.,0.]),
                        'sponge': np.array([+1., 0.,0.]),
                           }

boundaryTags = {'y-' : 1,
                'x+' : 2,
                'y+' : 3,
                'x-' : 4,
                'sponge':5,
               }

slope1=0.85
slope2=3.02
ymax=slope1+3.7

top = 1.0

vertices=[[0.0, 0.0],#0                                                                          
          [18.0-13.98, 0.0],#1                                                                         
          [25.0-13.98, slope1], #2
          [36.0-13.98, slope1],#3 
          [54-13.98, slope1],#4
          [61-13.98, slope1],#5                                                                 
          [87-13.98, slope2],#6                                                              
          [87-13.98, ymax],#7
          [0.0, ymax],]#8       
          #[-wavelength,ymax],#9
          #[-wavelength,0.0],]#10

vertexFlags=np.array([1, 1, 1, 1,
                      1, 1, 1, 
                      3, 3,])
                      #4, 4,])
    
segments=[[0,1],#0                                                                               
          [1,2],#1                                                                               
          [2,3],#2                                                                               
          [3,4],#3
          [4,5],#4                                                                               
          [5,6],#5                                                                               
          [6,7],#6
          [7,8],#7
          [8,0],]#8
          #[8,9],#9
          #[9,10], #10
          #[10,0],] #11
    
segmentFlags=np.array([1,1,1,1,1,1,
                       2,
                       3,
                       4])
#                       5,
#                       3,4,1])

regions = [[ 0.1, 0.1],]
          #[-0.1, 0.1],]

regionFlags=np.array([1])
#regionFlags=np.array([1,2])

if opts.mangrove_porous:
    vertices.append([36.0-13.98, 3.9]) # 11
    vertices.append([54.0-13.98, 3.9]) # 12
    vertexFlags=np.append(vertexFlags,5)
    vertexFlags=np.append(vertexFlags,5)
    segments.append([3,11])
    segments.append([11,12])
    segments.append([12,4])
    segmentFlags=np.append(segmentFlags,5)
    segmentFlags=np.append(segmentFlags,5)
    segmentFlags=np.append(segmentFlags,5)
    regions.append([45-13.98,0.86])
    regionFlags=np.append(regionFlags,3)

    print(len(segments),len(segmentFlags))


# TANK

tank = st.CustomShape(domain, vertices=vertices, vertexFlags=vertexFlags,
                      segments=segments, segmentFlags=segmentFlags,
                      regions=regions, regionFlags=regionFlags,
                      boundaryTags=boundaryTags, boundaryOrientations=boundaryOrientations)


# SPONGE LAYERS
# generation zone: 1 wavelength
# absorption zone: 2 wavelengths
#tank.setSponge(x_n=wavelength, x_p=2*wavelength)


#  ____                        _                   ____                _ _ _   _
# | __ )  ___  _   _ _ __   __| | __ _ _ __ _   _ / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
# |  _ \ / _ \| | | | '_ \ / _` |/ _` | '__| | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
# | |_) | (_) | |_| | | | | (_| | (_| | |  | |_| | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |____/ \___/ \__,_|_| |_|\__,_|\__,_|_|   \__, |\____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
#                                           |___/
# Boundary Conditions

# TANK

# atmosphere on top
tank.BC['y+'].setAtmosphere()
# free slip on bottom
tank.BC['y-'].setFreeSlip()
# free slip on the right
tank.BC['x+'].setFreeSlip()
# non material boundaries for sponge interface
tank.BC['sponge'].setNonMaterial()

# WAVE AND RELAXATION ZONES

smoothing = he*1.5
dragAlpha = 5*2*np.pi/wave_period/(1.004e-6)
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave,
                                               smoothing=smoothing,
                                               vert_axis=1)

#tank.setGenerationZones(flags=2,
#                   epsFact_porous=wavelength*0.5,
#                   center=[-wavelength*0.5,ymax*0.5],
#                   orientation=[1,0,0],
#                   waves=wave,
#                   dragAlpha=dragAlpha,
#                   vert_axis=1,
#                   porosity=1.,
#                   smoothing=smoothing)


if opts.mangrove_porous:
    tank.setPorousZones(flags=3,
                        dragAlpha=0,
                        dragBeta=0, #HD: 
                        porosity=0.375) #HD:0.973 LD:0.987

#column_gauge_locations=[((13.98,0.,0.),(13.98,ymax,0.)),((25.36,0.85,0.),(25.36,ymax,0.)),((28.68,0.85,0.),(28.68,ymax,0.)),((29.6,0.85,0.),(29.6,ymax,0.)),((30.52,0.85,0.),(30.52,ymax,0.)),((31.44,0.85,0.),(31.44,ymax,0.)),((32.04,0.85,0.),(32.04,ymax,0.)),((54.6,0.85,0.),(54.6,ymax,0.)),((55.81,0.85,0.),(55.81,ymax,0.)),((56.42,0.85,0.),(56.42,ymax,0.)),((57.33,0.85,0.),(57.33,ymax,0.)),((57.95,0.85,0.),(57.95,ymax,0.)),((65.17,0.85+(65.17-61)/(87-61)*(3.02-0.85),0.),(65.17,ymax,0.))]
column_gauge_locations=[((0.01,0.,0.),(0.01,ymax,0.)),((11.38,0.85,0.),(11.38,ymax,0.)),((14.7,0.85,0.),(14.7,ymax,0.)),((15.62,0.85,0.),(15.62,ymax,0.)),((16.54,0.85,0.),(16.54,ymax,0.)),((17.46,0.85,0.),(17.46,ymax,0.)),((18.06,0.85,0.),(18.06,ymax,0.)),((40.62,0.85,0.),(40.62,ymax,0.)),((41.83,0.85,0.),(41.83,ymax,0.)),((42.44,0.85,0.),(42.44,ymax,0.)),((43.35,0.85,0.),(43.35,ymax,0.)),((43.97,0.85,0.),(43.97,ymax,0.)),((51.19,0.85+(65.17-61)/(87-61)*(3.02-0.85),0.),(51.19,ymax,0.))]

tank.attachLineIntegralGauges('vof',
                              gauges=((('vof',), column_gauge_locations),),                     
                              fileName='column_gauges.csv') 

pressure_gauge_locations= ((35.89-13.98, 1.64, 0), (39.55-13.98, 1.64, 0),(43.10-13.98,1.62,0),(46.87-13.98,1.63,0),(50.52-13.98,1.63,0),(54.19-13.98,1.63,0))
tank.attachPointGauges('twp', gauges=((('p',), pressure_gauge_locations),), fileName='pressure_gaugeArray.csv')

velocity_gauge_locations=((32.24-13.98,1.25, 0), (43.09-13.98, 1.40, 0),(43.09-13.98,1.55,0),(43.09-13.98,1.72,0),(43.09-13.98,1.86,0),(57.83-13.98,1.38,0))
tank.attachPointGauges('twp', gauges=((('u','v'), velocity_gauge_locations),), fileName='velocity_gaugeArray.csv')





    
#  ___       _ _   _       _    ____                _ _ _   _
# |_ _|_ __ (_) |_(_) __ _| |  / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
#  | || '_ \| | __| |/ _` | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
#  | || | | | | |_| | (_| | | | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |___|_| |_|_|\__|_|\__,_|_|  \____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
# Initial Conditions

from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
smoothing = 1.5*he
nd = 2

class zero(object):
    def uOfXT(self,x,t):
        return 0.0

class P_IC:
    def uOfXT(self, x, t):
        p_L = 0.0
        phi_L = ymax - water_level
        phi = x[nd-1] - water_level
        p = p_L -g[nd-1]*(rho_0*(phi_L - phi)
                          +(rho_1 -rho_0)*(smoothedHeaviside_integral(smoothing,phi_L)
                                                -smoothedHeaviside_integral(smoothing,phi)))
        return p
class U_IC:
    def uOfXT(self, x, t):
        return 0.0
class V_IC:
    def uOfXT(self, x, t):
        return 0.0
class W_IC:
    def uOfXT(self, x, t):
        return 0.0
class VF_IC:
    def uOfXT(self, x, t):
        return smoothedHeaviside(smoothing, x[nd-1]-water_level)
class PHI_IC:
    def uOfXT(self, x, t):
        return x[nd-1] - water_level

# instanciating the classes for *_p.py files
initialConditions = {'pressure': P_IC(),
                     'vel_u': U_IC(),
                     'vel_v': V_IC(),
                     'vel_w': W_IC(),
                     'vof': VF_IC(),
                     'ncls': PHI_IC(),
                     'rdls': PHI_IC()}

#  __  __           _        ___        _   _
# |  \/  | ___  ___| |__    / _ \ _ __ | |_(_) ___  _ __  ___
# | |\/| |/ _ \/ __| '_ \  | | | | '_ \| __| |/ _ \| '_ \/ __|
# | |  | |  __/\__ \ | | | | |_| | |_) | |_| | (_) | | | \__ \
# |_|  |_|\___||___/_| |_|  \___/| .__/ \__|_|\___/|_| |_|___/
#                                |_|


domain.MeshOptions.genMesh = True
domain.MeshOptions.he = he
mesh_fileprefix = 'mesh'
domain.MeshOptions.setOutputFiles(mesh_fileprefix)

st.assembleDomain(domain)




#  _   _                           _
# | \ | |_   _ _ __ ___   ___ _ __(_) ___ ___
# |  \| | | | | '_ ` _ \ / _ \ '__| |/ __/ __|
# | |\  | |_| | | | | | |  __/ |  | | (__\__ \
# |_| \_|\__,_|_| |_| |_|\___|_|  |_|\___|___/
# Numerics

# outputStepping = TpFlow.OutputStepping(
#     final_time=T,
#     dt_init=dt_init,
#     # cfl=opts.cfl,
#     dt_output=sampleRate,
#     nDTout=None,
#     dt_fixed=None,
# )

outputStepping = TpFlow.OutputStepping()
outputStepping.final_time=T
outputStepping.dt_init=dt_init
outputStepping.dt_output=sampleRate
outputStepping.nDTout=None

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem()
myTpFlowProblem.domain=domain
#myTpFlowProblem.SystemPhysics.initialConditions=initialConditions
myTpFlowProblem.outputStepping=outputStepping



myTpFlowProblem.SystemNumerics.cfl=0.33
myTpFlowProblem.SystemNumerics.useSuperlu=True

myTpFlowProblem.SystemPhysics.setDefaults()
myTpFlowProblem.SystemPhysics.useDefaultModels()
myTpFlowProblem.SystemPhysics.gravity = np.array([0.0,-9.81,0.0])

m = myTpFlowProblem.SystemPhysics.modelDict
# ADD RELAXATION ZONES TO AUXILIARY VARIABLES


m['flow'].p.initialConditions['p'] = P_IC()
m['flow'].p.initialConditions['u'] = zero()
m['flow'].p.initialConditions['v'] = zero()
m['vof'].p.initialConditions['vof'] = VF_IC()
m['ncls'].p.initialConditions['phi'] = PHI_IC()
m['rdls'].p.initialConditions['phid'] = PHI_IC()
m['mcorr'].p.initialConditions['phiCorr'] = zero()

m['flow'].auxiliaryVariables = domain.auxiliaryVariables['twp']
m['vof'].auxiliaryVariables = domain.auxiliaryVariables['vof']
