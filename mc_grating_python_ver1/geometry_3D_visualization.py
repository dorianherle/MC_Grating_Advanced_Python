# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 10:39:48 2021

@author: Dorian
"""
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

import numpy as np
import math
# import matplotlib.patches as mpatches

import shapely
from shapely.ops import split
from shapely.geometry import LineString

from shapely.geometry import Polygon

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection

import trimesh

# For the visualization of the mesh
from mayavi import mlab
from mayavi.mlab import *
import mayavi
#mayavi.engine.current_scene.scene.off_screen_rendering = True
from mayavi.modules.text3d import Text3D

from tqdm import tqdm

def visualize_3D(scene, materials, color="", show_plot = True, show_axis=False, show_legends=True, show_vertices=False, save=""):

    # if no color is provided, assign all materials a random color
    if color == "":
        color = {material:(np.random.rand(1),np.random.rand(1),np.random.rand(1)) for material in materials}
        
    ## VISUALIZE MESH
    opacity = 1  
    # Create figure
    fig = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(0, 0, 0))
    # Iterate through the meshes
    xx = []
    yy = []
    zz = []
    for i,mesh in enumerate(scene):
        x,y,z = zip(*mesh.vertices) 
        triangular_mesh(x,y,z,mesh.faces, color=color[materials[i]], opacity=opacity)
        xx.extend(list(x))
        yy.extend(list(y))
        zz.extend(list(z))
        
    
    #mlab.points3d(x, y, z, opacity=0,transparent = True,scale_factor=1)
    x_min = np.min(xx)
    x_max = np.max(xx)
    y_min = np.min(yy)
    y_max = np.max(yy)
    z_min = np.min(zz)
    z_max = np.max(zz)
        
        
    if (show_axis == True):

        x = [x_min,x_max]
        y = [y_min,y_max]
        z = [z_min,z_max]

        axes = mlab.axes(color=(1, 1, 1), 
                        xlabel='x',
                        ylabel='y',
                        zlabel='z',
                        nb_labels = 5,
                        ranges = [x_min, x_max, y_min, y_max, z_min, z_max]
                        )
        axes.title_text_property.color = (1.0, 1.0, 1.0)
        axes.label_text_property.color = (1.0, 1.0, 1.0)
        axes.axes.font_factor = 0.70192
        
    if (show_legends == True):
 
        for i,material in enumerate(np.unique(materials)):
            mlab.text3d(x_min-100, 0, z_min+100*i,
                         material, scale = 10, color = color[material] ) 
            
    if (show_vertices == True):
        for i,mesh in enumerate(scene):
            x,y,z = zip(*mesh.vertices) 
            mlab.points3d(x, y, z, color=color[materials[i]])
    
   
    mlab.view(azimuth=270, elevation=90, roll=180, figure=fig)
    
    # Save plot
    if save != "":
        
        show_plot = False
        print("SAVING 3D PLOT")
        # mlab.options.offscreen = True
        mlab.savefig(save+".obj")
        
        
        mlab.close()
       
    # View plot
    if show_plot:
        mlab.show()
        

    
    
    return scene, materials

# Define Function that returns the points of a circle (anti-clockwise)
def PointsInCircum(r,pos=[0,0], n=500):
    points = [(math.cos(2*math.pi/n*x)*r,math.sin(2*math.pi/n*x)*r) for x in range(0,n+1)]
    x,y = zip(*points)
    x = np.array(list(x)) + pos[0]
    y = np.array(list(y)) + pos[1]
    return [x,y]

def PointsInRect(w,h,pos=[0,0]):
    points = [[w/2, -h/2],
              [w/2, h/2],
              [-w/2, h/2],
              [-w/2, -h/2],
              [w/2, -h/2]]
    
    x,y = zip(*points)
    x = np.array(list(x)) + pos[0]
    y = np.array(list(y)) + pos[1]
    
    return [x,y] 
    
    
def get_points_for_type(type_,geometry):
    
    if type_ == "circle":
        pos = geometry[0]
        diameter = geometry[1]
        return PointsInCircum(diameter/2,pos=pos)
    
    if type_ == "rectangle":
        pos = geometry[0]
        w = geometry[1]
        h = geometry[2]
        return PointsInRect(w,h,pos=pos)
    
    if type_ == "polygon":
        return geometry
    


def slice_boundary_and_MC_grating_symmetry(polygon, period_x, period_y):
    """
    Boundaries:
       3
     |---|
    0|   |2
     |---|
       1
       
    Parts Nomenclature defined by direction of normal:
               |
    ngative    |--> n positive
               |
    
    """
    dict_splits = {}
    
    
    lines = [
        LineString([(-period_x/2, -5e10), (-period_x/2, 5e10)]), #xleft
        LineString([(-5e10,-period_y/2), (5e10,-period_y/2)]), #ylower
        LineString([(period_x/2, -5e10), (period_x/2, 5e10)]),#xright
        LineString([(-5e10,period_y/2), (5e10,period_y/2)]), #yupper
        ]
    
    sliced = [[(polygon,"-1","o")]]
    for i, line in enumerate(lines):
        
        temp=[]
        for s in sliced[-1]:
            m,prev_i,prev_side = s
            split_res = split(m,line)
            if len(split_res) > 1:
                pi,ni = split_res
                pi_c = list(pi.centroid.centroid.coords)[0]
                ni_c = list(ni.centroid.centroid.coords)[0]
                
                if i == 0: #xleft
                    if pi_c[0] > ni_c[0]:
                        p = pi
                        n = ni
                    else:
                        p = ni
                        n = pi
                        
                if i == 1: #ylower
                    if pi_c[1] > ni_c[1]:
                        p = pi
                        n = ni
                    else:
                        p = ni
                        n = pi
                        
                if i == 2: #xright
                    if pi_c[0] < ni_c[0]:
                        p = pi
                        n = ni
                    else:
                        p = ni
                        n = pi
                        
                if i == 3: #yupper
                    if pi_c[1] < ni_c[1]:
                        p = pi
                        n = ni
                    else:
                        p = ni
                        n = pi
                        
                
                temp.extend([(p,prev_i+str(i),prev_side + "p"),\
                             (n,prev_i+str(i),prev_side + "n")])
                    
        if len(temp) > 0:
            sliced.append(temp)
            
        
    temp_m = []
    for s in sliced[-1]:
        m,prev_i,prev_side = s
        dict_splits[prev_i+prev_side] = m
        temp_m.append(m)
    
    return dict_splits

# Plots a Polygon to pyplot `ax`
def plot_polygon(ax, poly, **kwargs):
    path = Path.make_compound_path(
        Path(np.asarray(poly.exterior.coords)[:, :2]),
        *[Path(np.asarray(ring.coords)[:, :2]) for ring in poly.interiors])

    patch = PathPatch(path, **kwargs)
    collection = PatchCollection([patch], **kwargs)
    
    ax.add_collection(collection, autolim=True)
    ax.autoscale_view()
    return collection
    
def layer_polygons(dict_):

    period_x, period_y = dict_["general"]["period"]
    layer_dict = {}
    
    # cover
    cover_material = dict_["cover"]
    print(cover_material)
    if not cover_material.startswith("Air"):
    
        x,y = get_points_for_type("rectangle",
                                  [[0,0], period_x, period_y])
        
        polygon = Polygon(zip(x,y))
        
        layer_dict["cover"] = ([[polygon, cover_material]]) 
        

    
    layers = [k for k in dict_.keys() if k.isnumeric()]
    for layer in layers:
        list_polygons = []
        
        surrounding_material = dict_[layer]["surrounding_material"]
        
        if not surrounding_material.startswith("Air"):
            
            x,y = get_points_for_type("rectangle",
                                       ([0,0], period_x, period_y))
            polygon = Polygon(zip(x,y))
            list_polygons.append((polygon, surrounding_material))
            
    
        pillars = [k for k in dict_[layer].keys() if k.isnumeric()]
        
        for pillar in pillars[::]:
            
            x,y = get_points_for_type(dict_[layer][pillar]["type"],
                                      dict_[layer][pillar]["geometry"])
            polygon = Polygon(zip(x,y))
            
            
            # slice mesh on boundary and polygon
            test = slice_boundary_and_MC_grating_symmetry(polygon, 
                                                          period_x, period_y)
            
            for i,k in enumerate(test.keys()):
                m = test[k]
               
                if "n" in k:
                    planes=[a for a in k[2:] if a.isdigit()]
                    sides =[a for a in k[2:] if not a.isdigit() and "o" not in a]
                    move_planes = np.array(planes)[np.array(sides)=="n"]
                    
                    
                    for p in move_planes:
                        
                        if p == "0":
                            m = shapely.affinity.translate(m, xoff=period_x, 
                                                           yoff=0.0, 
                                                           zoff=0.0)
                        if p == "1":
                            m = shapely.affinity.translate(m, xoff=0, 
                                                           yoff=period_y, 
                                                           zoff=0.0)
                        if p == "2":
                            m = shapely.affinity.translate(m, xoff=-period_x, 
                                                           yoff=0, 
                                                           zoff=0.0)
                        if p == "3":
                            m = shapely.affinity.translate(m, xoff=0, 
                                                           yoff=-period_y, 
                                                           zoff=0.0)
    
                    
                else:
                    pass
                    
            
                list_polygons.append((m, dict_[layer][pillar]["material"]))
       
        layer_dict[layer] = list_polygons
        
    
    # only keep polygons that are not air 
    for layer, list_polygons in layer_dict.items():

        new_list_polygons = []
        for i, target_polygon in enumerate(list_polygons):
            target_p, target_m = target_polygon
            for polygon in list_polygons[i+1:]:
                p,_ = polygon
                target_p = target_p.difference(p)
                
            if not target_m.startswith("Air"):
                new_list_polygons.append((target_p, target_m))
        
        layer_dict[layer] = new_list_polygons
        
    
    # Add substrate
    substrate_material = dict_["substrate"]
    if not substrate_material.startswith("Air"):

        x,y = get_points_for_type("rectangle",
                                  [[0,0], period_x, period_y])
        
        polygon = Polygon(zip(x,y))
        
        layer_dict["substrate"] = ([[polygon, substrate_material]]) 
  
    
  
    # if plot:
        # for layer in layer_dict.keys():
            
            # fig, ax = plt.subplots(dpi=300)
            # plt.title("Layer: " + str(layer))
                
            # materials = []
            # for pillar in layer_dict[layer]:
                # polygon, material = pillar
                # materials.append(material)
                # plot_polygon(ax, polygon, 
                             # facecolor=color[material], 
                             # edgecolor='none')
                
            # plt.axis("equal")
            
            # # add legend
            # legend_patches = []
            # for mat in np.unique(materials):
                # color_patch = mpatches.Patch(color=color[mat], label=mat)
                # legend_patches.append(color_patch)
            
            # ax.legend(handles=legend_patches, loc='center left', bbox_to_anchor=(1, 0.5))
      
    
        
    return layer_dict

# For the creation fo the  mesh

def join_common_meshes(l):
    # join common : https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements
    
    out = []
    while len(l)>0:
        first, *rest = l
        first = set(first)
    
        lf = -1
        while len(first)>lf:
            lf = len(first)
    
            rest2 = []
            for r in rest:
                if len(first.intersection(set(r)))>0:
                    first |= set(r)
                else:
                    rest2.append(r)     
            rest = rest2
    
        out.append(first)
        l = rest
        
    return out

import matplotlib.pyplot as plt

def create_geometry_mesh(dict_, union=True):
    layer_dict = layer_polygons(dict_)
    
    
    # CONVERT POLYGONS TO MESHES
    materials = {}
    materials_list = []
    # current layer position (in z)
    layer_pos = 0
    # intiialize mesh dictionary
    mesh_dict = {}
    for layer in layer_dict.keys():
        temp_materials = []
        temp = []
        if layer == "cover" or layer == "substrate": # cover
            tickness = 100
        else:
            tickness = dict_[layer]["thickness"]
            
        for pillar in layer_dict[layer]:
            polygon, material = pillar
            temp_materials.append(material)
            mesh = trimesh.creation.extrude_polygon(polygon, height = tickness)
            mesh.apply_translation((0,0,layer_pos))
            #mesh.visual.vertex_colors = color[material]
            temp.append(mesh)
            
        # update layer poistion (in z)
        layer_pos += tickness
        materials_list.extend(temp_materials)
        materials[layer] = np.unique(temp_materials) 
        mesh_dict[layer] = temp
        
    # IDENTIFY WHICH MESHES SHOULD BE JOINED
    nbr_of_layers = len(layer_dict.keys())
    all_mehes = []
    to_join=[]
    for layer in layer_dict.keys():
        
        if layer != "cover" and layer != "substrate": 
            temp_join = []
            
            for material in materials[layer]:
                # print("CURRENT LAYER: ", layer)
                # print("MATERIAL: ", material)
                # gather all polygons of this material
                ps_layer0 = [(i,pm[0]) for i,pm in enumerate(layer_dict[layer]) if pm[1] == material]
                all_mehes.append([layer+"_"+str(i) for i,_ in ps_layer0])
                # gather all poylgons of next layer that share the same material
                next_layer = int(layer)+1
                # print("NEXT LAYER: ", next_layer)
                if next_layer < nbr_of_layers+1:
                    try: #TODO Figure out issue with numbering
                        ps_layer1 = [(i,pm[0]) for i,pm in enumerate(layer_dict[str(next_layer)]) if pm[1] == material]
                        for i,pi in ps_layer0:
                            for j,pj in ps_layer1:
                                
                                # show intersection
                                # xi, yi = pi.exterior.xy
                                # xj, yj = pj.exterior.xy
                                # plt.plot(xi, yi, c="red", label=layer+"_"+str(i))
                                # plt.plot(xj, yj, c="k", label=str(next_layer)+"_"+str(j))
                                # plt.legend()
                                # plt.pause(0.01)
                                # plt.show()
                              
                                    
                                if pi.intersects(pj) and True:
                                    temp_join.append([layer+"_"+str(i),str(next_layer)+"_"+str(j)])
                    except:
                        pass

            to_join.extend(temp_join)         
    #print(to_join)
    # next step
    # join common : https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements
    joined = join_common_meshes(to_join)
    #print(joined)
    materials_list = []
    scene = []         
    for j in tqdm(joined):
        meshes = []
        for name in j:
            layer = name.split("_")[0]
            index = int(name.split("_")[1])
            meshes.append(mesh_dict[layer][index])
            
          
        if union != True:
            scene.append(trimesh.util.concatenate(meshes))
            materials_list.extend([layer_dict[layer][index][1]]) #*len(meshes)
        else:
            
            mesh = meshes[0]
            for m in meshes[1:]:
                mesh = mesh.union(m)
                
            scene.append(mesh)
            materials_list.append(layer_dict[layer][index][1])
    
    #remaining meshes
    all_m = [j for sub in all_mehes for j in sub]
    joined_m = [j for sub in joined for j in sub]
    
    rem_m = list(set(all_m)-set(joined_m))
    
    for rm in rem_m:

        layer = rm.split("_")[0]
        index = int(rm.split("_")[1])
        
        mesh = mesh_dict[layer][index]
        materials_list.append(layer_dict[layer][index][1])
        scene.append(mesh)
        
    # add cover and substrate
    try:
        mesh = mesh_dict["cover"][0]
        scene.append(mesh)
        cover_material = dict_["cover"]
        materials_list.append(cover_material)
    except:
        pass
    
    try:
        mesh = mesh_dict["substrate"][0]
        scene.append(mesh)
        substrate_material = dict_["substrate"]
        materials_list.append(substrate_material)
    except:
        pass

    return scene, materials_list


    

        
