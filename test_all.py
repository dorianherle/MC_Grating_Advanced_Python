#import sys
#sys.path.append(r"D:\Organized_PC\Simulation\MC_Grating_Python_ver2\mc_grating_python_new\mc_grating_python") 
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

from mc_grating_python_ver1.mc_grating_python import mc_grating_python as mc

#%%
# 1 | DEFINE GEOMETRY
 
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

#%%
        
# 2 | VISUALIZE GEOMETRY
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

# %%

# 3 | SIMULATE GEOMETRY 

sim = mc.sim(geometry_dict, single_wavelengh=400, nbr_oders_x=11, nbr_orders_y=11)


# %%
# 3.1 | FAR FIELD

ff = sim.far_field(parameter_row = "wavelength", \
                   start_row=400, end_row=700, nbr_of_points_row=20)
    
res_ff = mc.run(ff)

# %%
# 3.2 | NEAR FIELD

nf = sim.near_field(start_row=-200, end_row=200, plane="zy",  number_of_points = [200,50],\
                    position_of_plane=0.5, output_format="amplitude_re_im_and_power_flow")
 
res_nf = mc.run(nf)

#%%

# 4 | PLOT RESULTS

# 4.1 | PLOT FAR FIELD

p_ff = mc.plot(res_ff)
rta = p_ff.rta() # reflection - transmission - absorption plot

#%%
# 4.2 | PLOT NEAR FIELD

p_nf = mc.plot(res_nf, geometry_dict) 
ey_re = p_nf.field(component="Sz", plot=True, title="")
abs_e = p_nf.field(component="abs(E)", plot=True, title="|E|")
ey_re = p_nf.field(component="Ex mod", plot=True, title="")
ey_re = p_nf.field(component="Im(eps) |E|^2", plot=True, title="", periodicity=1)


