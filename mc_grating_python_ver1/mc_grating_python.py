# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 10:52:53 2021

@author: Dorian Herle
"""
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

print("MC GRATING PYTHON\nAuthor: Dorian Herle")


## -------- GEOMETRY ---------- ##

from geometry_3D import geometry3D

## ------ VISUALIZATION ------- ##
from visualizations import visualization


## ------- SIMULATION --------- ##
from analyze import analyze
from dict_to_mc_grating_script import dict_to_mc_grating_script
from run_simulation import launch_simulation

## ---------- ANALYZE ---------- ##
from visualizations import plot_res

## ---------- SAVE OPTIONS ---------- ##
# TODO

## ---------- ADVANCED SIM ---------- ##
# from advanced_sim import advanced_sim



class mc_grating_python():
    
    def __init__():
        pass
    
    def geo(periodx, periody):
        return geometry3D(periodx, periody)

    def show(geometry, color=""):
        return visualization(geometry, color=color)
    
    def sim(geometry, nbr_oders_x=11, nbr_orders_y=11,\
            polarization_type="Ex", angle_n = 0, angle_p = 0,\
            single_wavelengh=450):
 
        return analyze(geometry, 
                       nbr_oders_x=nbr_oders_x, 
                       nbr_orders_y=nbr_orders_y,
                       polarization_type=polarization_type,
                       angle_n = angle_n,
                       angle_p = angle_p,
                       single_wavelengh = single_wavelengh
                       )
    
    def run(analyze_dict, path="mc_grating_simulation_file", fix = False):
        a = dict_to_mc_grating_script(analyze_dict)
        mc_grating_str = a.gui_input_structure()+a.analyze_section()
        res = launch_simulation(mc_grating_str, path=path, fix=fix)
        return res
    
    def plot(res, geometry="", color=""):
        return plot_res(res, geometry=geometry, color=color)
    

   