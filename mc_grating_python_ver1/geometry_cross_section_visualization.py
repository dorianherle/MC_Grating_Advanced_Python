# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 14:53:26 2021

@author: Dorian
"""

# For the visualization of the mesh
from mayavi import mlab
from mayavi.mlab import *
import mayavi
#mayavi.engine.current_scene.scene.off_screen_rendering = True
from mayavi.modules.text3d import Text3D

#%% Cross-section
import matplotlib.pyplot as plt
import numpy as np
from descartes import PolygonPatch

import trimesh
import utils
import itertools
from scipy.spatial.transform import Rotation as R
#from geometry_3D_visualization import visualize_3D, create_geometry_mesh



def create_cross_section(scene, materials, pos=0.5, normal="X", cross_section="", color="", show_color = True, plot=True, save=""):
    
    # Quick fix for numerical error
    if pos==1: pos = 0.9999
    if pos==0: pos = 0.0001
    
    # if no color is provided, assign all materials a random color
    if color == "":
        color = {material:(np.random.rand(1)[0],np.random.rand(1)[0],np.random.rand(1)[0],1) for material in materials}
        
   
    combined = trimesh.util.concatenate(scene)
    xmin,ymin,zmin = combined.bounds[0]
    xmax,ymax,zmax = combined.bounds[1]
    
    if normal.upper() == "X":
        index = 0
        plane_pos = xmin+(xmax-xmin)*pos
        xlabel = "Y"
        ylabel= "Z"
    if normal.upper() == "Y":
        index = 1
        plane_pos = ymin+(ymax-ymin)*pos
        xlabel = "X"
        ylabel= "Z"
    if normal.upper() == "Z":
        index = 2
        plane_pos = pos #zmin+(zmax-zmin)*pos
        xlabel = "X"
        ylabel= "Y"
        
    
    plane_normal = [0,0,0]
    plane_origin = [0,0,0]
    
    plane_normal[index] = 1
    plane_origin[index] = plane_pos
    
    # TODO
    # if cross_section != "":
    #     normal = "XYZ"
    #     normal = normal.replace(cross_section[0].upper(), "")
    #     normal = normal.replace(cross_section[1].upper(), "")
        
    #     plane_normal = [0,0,0]
    #     index = "XYZ".index(normal)
    #     plane_normal[index] = 1
        
        
    # rotate such that cross-section is upright
    rot_scene= []
    for mesh in scene:
        new_vert = []
        for v in mesh.vertices:
            new_vert.append(rotate_vector(v, rotation_axis = np.array(plane_normal)))
    
        mesh = trimesh.Trimesh(vertices= new_vert,
                               faces=mesh.faces)
        
        rot_scene.append(mesh)
        
    
    cross_section = {}
    for i,mesh in enumerate(rot_scene):
        
      
        slice_ = mesh.section(plane_origin=plane_origin, 
                              plane_normal=plane_normal)
        
        if slice_ is not None: # mesh might be outside of plane bounds
    
            # transformation matrix for to_planar 
            # I don't know why
            to_2D = trimesh.geometry.align_vectors(plane_normal, [0,0,1])
            slice_2D,_ = slice_.to_planar(to_2D = to_2D)
            polygon = slice_2D.polygons_full
            cross_section[i] = {"polygon": polygon, "material": materials[i]}
            
            
       
    
    if plot:
        fig, ax = plt.subplots(dpi=300)
        
        keep_track_of_material = []
        for k,item in cross_section.items(): 
           
            geom = item["polygon"]
            material = item["material"]
            
            for i, polygom in enumerate(geom):  
                
                if show_color and i==0:
                    patch = PolygonPatch(polygom.buffer(0),\
                                         fc=color[material], ec=color[material], alpha = 1,\
                                         label=material)
                    
                if show_color and i > 0 or material in keep_track_of_material:
                    patch = PolygonPatch(polygom.buffer(0),\
                                         fc=color[material], ec=color[material], alpha = 1)
                    
                if not show_color:
                    patch = PolygonPatch(polygom.buffer(0),\
                                         fc="none", ec='k', alpha = 1)
                    
                    
                ax.add_patch(patch)
                
                keep_track_of_material.append(material)
                
        
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        if show_color == True:
            ax.legend(loc='center left', bbox_to_anchor=(0, -0.25))
        if normal == "Z":
            ax.set_ylim(ymin,ymax)
            ax.set_xlim(xmin,xmax)
        if normal == "Y":
            ax.set_ylim(zmin,zmax)
            ax.set_xlim(xmin,xmax)
            ax.invert_yaxis()
        if normal == "X":
            ax.set_ylim(zmin,zmax)
            ax.set_xlim(ymin,ymax)
            ax.invert_yaxis()
        ax.axis("equal")
        print("BBB")
        ax.invert_yaxis()
        
        if save!= "":
            plt.savefig(save+".svg", bbox_inches ='tight')
        plt.show()
        
    
    return cross_section, combined.bounds



def vis(scene, title=""):
    "visualize mesh"
    # Create figure
    fig = mlab.figure(title, fgcolor=(0, 0, 0), bgcolor=(0, 0, 0))
    # Iterate through the meshes
    
    a = np.linspace(0.1,1,len(scene))
    for i,mesh in enumerate(scene):
        x,y,z = zip(*mesh.vertices) 
        triangular_mesh(x,y,z,mesh.faces, color=(a[i],a[i],a[i]), opacity=1)

    
    mlab.view(azimuth=270, elevation=90, roll=180, figure=fig)
    # View plot
    mlab.show()
    
# %% TURN 
def rotate_vector(vec, rotation_degrees=-90, rotation_axis = np.array([0, 1, 0])):

    rotation_radians = np.radians(rotation_degrees)
    rotation_vector = rotation_radians * rotation_axis
    rotation = R.from_rotvec(rotation_vector)
    rotated_vec = rotation.apply(vec)
    
    return rotated_vec
    




    

