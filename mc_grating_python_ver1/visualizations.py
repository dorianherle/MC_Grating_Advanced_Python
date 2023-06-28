# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 11:21:37 2021

@author: Dorian
"""

import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

from geometry_3D_visualization import visualize_3D, create_geometry_mesh
from geometry_cross_section_visualization import create_cross_section
from plot_results import orders, rta, far_field_scattering, field, field_line, poynting_vector, poynting_vector_plot, geometry_cross_section, plot_2D


import trimesh

class visualization():
    
    def __init__(self, geometry, color=""):
        # get 3D meshes and associated materials
        self.scene, self.materials = create_geometry_mesh(geometry)
        self.color = color
    
    def show_3D(self, color="", show_plot = True,show_axis=False, show_legends=True, save =""):

        return visualize_3D(self.scene, \
                            self.materials,\
                            color=self.color,\
                            show_axis=show_axis,\
                            show_legends=show_legends,\
                            show_plot = show_plot,\
                            save=save)
    
    def show_cross_section(self, pos=0.5, normal="x", show_color = True, plot=True, save=""):
        return create_cross_section(self.scene,\
                                    self.materials,\
                                    pos=pos, normal=normal, color=self.color, show_color = show_color, plot=plot, save=save)
            
            
        
class plot_res():
    def __init__(self, res, geometry="", color=""):
        
        self.res = res
        self.geometry_b = False
        self.color_b = False
        self.scene = None
        self.materials = None
        
        if geometry != "":
            # get 3D meshes and associated materials
            self.scene, self.materials = create_geometry_mesh(geometry)
            self.geometry_b = True
            
        if color != "": 
            self.color = color
            self.color_b = True
        
    def geometry_cross_section(self, ax=None):
        return geometry_cross_section(self.res, self.scene, material=self.materials, ax=ax)
    
    def orders(self, direction="cover", total=False, plot=True):
        return orders(self.res, direction=direction, total=total, plot=True)
    
    def rta(self, plot="True", save="", exclude=None):
        return rta(self.res, plot=plot, save=save, exclude=exclude)
    
    
    def far_field_scattering(self, param, direction="cover", plot=True, title=""):
        return far_field_scattering(self.res, direction=direction, plot=plot, title=title)
    
    def plot_2D_scan(self, data="Total Power", cover = True, title="", vmin=0, vmax=1, plot=True):
        return plot_2D(self.res, data=data, cover = cover, title=title, vmin=vmin, vmax=vmax, plot=plot)
    
    def field(self, component="Sx", plot=True, title="", min_field="", max_field="", epsR=1, periodicity=1):
        if self.geometry_b: 
            scene = self.scene
            material = self.materials
        else:
            scene = ""
            material = ""
            
        return field(self.res, scene, material, component=component, plot=plot, title=title, min_field=min_field, max_field=max_field, epsR=epsR, periodicity=periodicity)
    
    
    def field_line(self,component="Sx", plot=True):
        return field_line(self.res,component=component, plot=plot)
    
    def poynting_vector_plot(self, field_sx, field_sz, points, skip_every=8,
                title = "", xlabel="", ylabel="", ax = None,  start_points = None):
        
        return  poynting_vector_plot(field_sx,field_sz,points, skip_every=skip_every,
                title = title, xlabel=xlabel, ylabel=ylabel, ax = ax,  start_points = start_points)
    
    def poynting_vector(self, plot=True, plot3D=True, title="", min_field="", max_field=""):
        if self.geometry_b: 
            scene = self.scene
            material = self.materials
        else:
            scene = ""
            material = ""
            
        return poynting_vector(self.res, scene, material, plot=plot, plot3D=plot3D,\
                               title=title, min_field=min_field,\
                               max_field=max_field)
        
    