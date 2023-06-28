# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 16:46:34 2021

@author: Dorian
"""
import os 
import sys
file_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(file_dir)

import pandas as pd
import numpy as np
from geometry import geometry_class as geometry
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.ops import cascaded_union

class geometry3D():
    
    def __init__(self, period_x, period_y):
        self.period_X = period_x
        self.period_Y = period_y
        
        # intialize
        self.objects = []
        self.geo = geometry(period_x=period_x, period_y=period_y) 

    def cover(self, material):
        self.geo.cover(material)
        
    def substrate(self, material):
        self.geometry_to_mc_grating_script()
        self.geo.substrate(material)
            
    # OBJECTS
    def cylinder(self, position, diameter, height, name = "cylinder", material = "Silicon (Table)", order = 0):
        """
        Creates a 3D cylinder object
    
        Parameters
        ----------
        position : list
            Position of bottom center of cylinder. 
            Example: position = [0,0,0].
        radius : float
            Radius of cylinder.
            Example: radius = 49.
        height : float
            Height of cylinder 
            Example: height = 77.
    
        Returns
        -------
        None.
    
        """
        
        # points where cross-section changes - z Axis
        points_cross_section_change = [position[2], position[2]+height]
        
        # Create cylinder object -> As dictionnary 
        dict_ = {}
        dict_["name"] = name
        dict_["type"] = 'cylinder'
        dict_["parameters"] = [position, diameter, height]
        dict_["material"] = material
        dict_["points_cross_section_change"] = points_cross_section_change
        dict_["order"] = order
        
        self.objects.append(dict_)
        


    def cuboid(self, position, length, width, height, name = "cuboid", material = "Silicon (Table)", order = 0):
        """
        Creates a 3D cuboid
        
        """
        
        # points where cross-section changes - z Axis
        points_cross_section_change = [position[2], position[2]+height]
        
        # Create cylinder object -> As dictionnary 
        dict_ = {}
        dict_["name"] = name
        dict_["type"] = 'cuboid'
        dict_["parameters"] = [position, length, width, height]
        dict_["material"] = material
        dict_["points_cross_section_change"] = points_cross_section_change
        dict_["order"] = order
        
        self.objects.append(dict_)
        
        
    def print_g(self):
        return self.geo.print_g()
        

    #------------------------------------------------------------------------#
    def object_to_layer(self, object_, layer):
        
        # Cylinder
        if object_["type"] == 'cylinder':
            [position, radius, height] = object_["parameters"]
            material = object_["material"]
            
            layer.circle(position[:-1], radius, material)
            
        
        # Cuboid
        if object_["type"] == 'cuboid':
            [position, length, width, height] = object_["parameters"]
            material = object_["material"]
            
            
            layer.rectangle(position[:-1], length, width, material)
            
    def geometry_to_mc_grating_script(self):
      
        objects = self.objects
        geometry_objects = {}
        for obj in objects:
            geometry_objects[obj["name"]]=obj
            
        
        names_list = []
        points_cross_section_change_list = []
        for object_ in objects:
            for point_ in object_["points_cross_section_change"]:
                names_list.append(object_["name"])
                points_cross_section_change_list.append(point_)
                
        
        
        df = pd.DataFrame([names_list, points_cross_section_change_list], index=['Name', 'Cross-Section-Change']).T
        df = df.sort_values('Cross-Section-Change')  

        
        # Combine all objects that belong to the same layer
        # Source: https://stackoverflow.com/questions/65740018/pandas-dataframe-regrouping/65740157#65740157
        df = (df.pivot(index='Cross-Section-Change', columns='Name', values='Name')
           .apply(lambda x: x[x.notna().cumsum().ne(0)].bfill())
           .apply(lambda x: list(x.dropna()), axis=1)
        )
        
        # Get layer thicknesses and objects
        # From top to bottom -> MC Grating Layer Definition
        layer_thicknesses = np.abs(np.diff(df.index.tolist()[::-1]))
        
        objects_at_each_cross_section_change = list(df.values)[::-1]
        # Get object in each layer -> Keep strings which are the same in subsequent lists
        objects_in_each_layer = self.string_diff(objects_at_each_cross_section_change)
        
        # order objects in each layer according to their order -> important to make holes etc.
        sorted_objects_in_each_layer = []
        for objects_in_layer in objects_in_each_layer:
        
            orders = [geometry_objects[object_name]["order"] for object_name in objects_in_layer]
        
            # sort
            sorted_objects_in_layer = [x for _,x in sorted(zip(orders,objects_in_layer))][::-1]
            
            sorted_objects_in_each_layer.append(sorted_objects_in_layer)
            
        
        # Build MC Grating Script
        for i,layer_thickness in enumerate(layer_thicknesses):
            
            # check if cuboid in span domain
            
            check = [self.check_if_cuboid_spans_doamin(geometry_objects[object_name]) for object_name in sorted_objects_in_each_layer[i]]
            
            # check rare case that all pillars are from the same material
            # all_pillars_span_domain_same_material = False
            # materials = [geometry_objects[object_name]["material"] for object_name in sorted_objects_in_each_layer[i]]
            # if len(np.unique(materials)) == 1:
            #     # and all togehter span the entire domain
            #     if self.check_if_takes_entire_area([geometry_objects[object_name] for object_name in sorted_objects_in_each_layer[i]]):
            #         # LAYER and NO PILLARS
            #         layer = self.geo.layer(thickness = layer_thickness, surrounding_material = materials[0])
            #         all_pillars_span_domain_same_material = True
                    
         
            # BACKGROUND MATERIAL (LAYER) 
            if True in check:
                name = np.array(sorted_objects_in_each_layer[i])[check]
            
                if len(name) > 1:
                    raise Exception('Perfectly overlapping cuboids found in layer '+ str(i)+". Object names: " + str(name) +"\nPlease fix.")
                else:
                    name_of_cuboid_that_spans_domain = name[0]

                    
                object_ = geometry_objects[name_of_cuboid_that_spans_domain]
                #[_, _, _, height] = object_["parameters"]
                material = object_["material"]
                # LAYER 
                layer = self.geo.layer(thickness = layer_thickness, surrounding_material = material)
                
                # ADD "PILLARS" to LAYER
                for object_name in sorted_objects_in_each_layer[i]:
                    # if object does not span the domain 
                    if object_name != name_of_cuboid_that_spans_domain:
                        self.object_to_layer(geometry_objects[object_name], layer)
                    
                
            else:
                # LAYER 
                layer = self.geo.layer(thickness = layer_thickness, surrounding_material = "Air (Special Formula)")
                
                # ADD "PILLARS" to LAYER
                for object_name in sorted_objects_in_each_layer[i]:
                    self.object_to_layer(geometry_objects[object_name], layer)
                    
    def check_if_cuboid_spans_doamin(self, object_):
        # first check if imput object is a cuboid ...
        if object_["type"] == 'cuboid':
            # Check if the Cuboid spans the entire domain. If it is the case, then
            # it is better to define it as a layer than a rectangle
            [position, length, width, height] = object_["parameters"]
            pos_X, pos_Y, pos_Z = position
            if pos_X == 0 and pos_Y == 0 and length == self.period_X and width == self.period_Y:
                return True
            else:
                return False
        else:
            return False
        
        
    def check_if_takes_entire_area(self, objects_):
        print(objects_)
        polygons = []
        for object_ in objects_:
            type_ = object_["type"]
            geometry = object_["parameters"]
            
            if type_ ==  "circle":
                center = geometry[0]
                diameter = geometry[1]
                polygons.append(Point(center[0], center[1]).buffer(diameter/2))
            
            if type_ == "rectangle":
                center = geometry[0]
                width = geometry[1]
                height = geometry[2]
                # Define the outer border
                border = [(center[0]-width/2, center[0]-height/2), 
                          (center[0]+width/2, center[0]-height/2), 
                          (center[0]+width/2, center[0]+height/2), 
                          (center[0]-width/2, center[0]+height/2)]
                polygons.append(Polygon(shell=border))
                
            if type_ == "polygon":
                polygons.append(Polygon([[p.x, p.y] for p in geometry]))
                
        # merge all polygons
        u = cascaded_union(polygons)
        # create unit_cell polygon
        center = [0,0]
        width = self.period_X
        height = self.period_Y
        border = [(center[0]-width/2, center[0]-height/2), 
                   (center[0]+width/2, center[0]-height/2), 
                   (center[0]+width/2, center[0]+height/2), 
                   (center[0]-width/2, center[0]+height/2)]
        unitCell = Polygon(shell=border)
        # check overlab
        overlab = unitCell.intersection(u) 
         
        # if overlab is more than 99%
        if overlab.area/unitCell.area*100 > 99:
             return True
        else:
             return False
        
    def string_diff(self, list_of_list_of_strings):
        # When you have a list of list of strings 
        # and you wan to create a list of list of strings where only the strings are
        # inserted that are the same in subsequent lists
        # Example: temp1 = [['One', 'Two'],['One', 'Two','Three', 'Four']]
        # Res: [['One', 'Two']]
        res = []
        start = list_of_list_of_strings[0]
        for list_ in list_of_list_of_strings[1:]:
            s = set(start)
            temp = [x for x in list_ if x in s]
            res.append(temp)
            # update
            start = list_
        return res 

    


        
if __name__ == "__main__":  
    periodx, periody = [440, 440]       
    geo = geometry3D(periodx, periody)
    geo.cover("Air")
    
    mp = 100
    #geo.cylinder([0,0,0], 10, 100, name="c1")
    geo.cylinder([0,0,0], 60, 100, name="c1", material="Al (Table GS)", order = 1)
    geo.cylinder([0,0,mp], 70, 10, name="air_p", material="Air (Special Formula)", order = 2)
    geo.cylinder([0,0,0], 30, 100, name="air_p1", material="Air (Special Formula)")
    geo.cuboid([0,0,mp], 100, 100, 10, order = 3)
    #geo.cylinder([0,0,0], 10, 100, name="c1")
    geo.substrate("Si")
    a = geo.print_g()
    
    #mp = 200
    # #geo.cylinder([0,0,0], 10, 100, name="c1")
    # # geo.cylinder([0,0,0], 80, 300, name="oxide_stand", material="SiO2 (Table GS)")
    # # geo.cylinder([0,0,300], 280, 100, name="disk", material="Silicon (Table)")
    # geo.cylinder([0,0,mp], 380, 100, name="membrane_hole", material="Air (Special Formula)", order = 1)
    # geo.cylinder([0,0,mp], 500, 100, name="membrane", material="Al (Table GS)", order = 2)
    
    # geo.cylinder([periodx,periody,mp], 380, 100, name="membrane_hole1", material="Air (Special Formula)", order = 1)
    # geo.cylinder([periodx,periody,mp], 500, 100, name="membrane1", material="Al (Table GS)", order = 2)

    
    
    #geo.cylinder([0,0,0], 10, 100, name="c1")
    geo.substrate("Si")
    a = geo.print_g()
    print(a)
    
    # from visualize_geometry_vpython import visualize_geometry
    
    # visualize_geometry(a)
    
# %%

        
    
            
        
                            