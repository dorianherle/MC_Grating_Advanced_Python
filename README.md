# MC_Grating_Advanced_Python

```
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

from mc_grating_python_ver1.mc_grating_python import mc_grating_python as mc
```

## Create a geometry

```
def geometry(parameters):

    if parameters == "": # TO TEST
        parameters = [137,98,77]
        
    if len(parameters) != 3 and parameters != "":
        raise ValueError("List of input parameters needs to be 3")
              
    period, pillar_d, pillar_h\
        = parameters

    # geometry 
    geo = mc.geo(period, period)
    
    # Cover Layer
    geo.cover("Air (Special Formula)")
    
    #---PILLAR---#
    geo.cylinder(position=[0,0,0], diameter=pillar_d, height=pillar_h, 
                 name = "pillar", material = "Silicon (Table)", order = 1)
                 

    # Substrate
    geo.substrate("Fused Silica (Sellmeier)")
        
    
    # get geometry
    geometry = geo.print_g() 
    
    return geometry
```

## Visualize the Geometry

```
color={
           'Amorphous_Si_CMI (Table GS)': (0.3,0.3,0.3),
           'Amorphous_Si_CMI_Lossless': (0.3,0.3,0.1),
           'Silicon (Table)' : (0.4,0.4,0.4),
           'Fused Silica (Sellmeier)':(0,0,1),
           "SiO2 (Table GS)": (0.7,0.1,0.7),
           "Al (Table GS)": (0.1,0.5,0.1),
           "aSi_k_half": (0,1,0),
           "Al2O3 (Table GS)": (0,0,0.5),
           "Air (Special Formula)": (1,1,1),
            }


geometry_dict = geometry("")
visual = mc.show(geometry_dict, color=color)
cs = visual.show_cross_section(normal="y")
```
