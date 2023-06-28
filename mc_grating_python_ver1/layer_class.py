# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 16:04:46 2021

@author: Dorian
"""
class layer_class():
    
    def __init__(self, thickness, material, geo_dict):
        
        # Initialize geometry dictionarry
        self.thickness = thickness
        self.material = material
        self.geo_dict = geo_dict
        
        # Find latest added layer
        self.layer_nbr = str(max([int(k) for k in self.geo_dict.keys() if k.isnumeric()]))
        
        # Initialize shapes counter (MC_GRATING PILLAR)
        self.shapes_counter = 0
        
        
    def circle(self, pos_center, diameter, material):
        
        # update shapes counter
        self.shapes_counter += 1
        
        # update dict
        self.geo_dict[self.layer_nbr][str(self.shapes_counter)] = \
            {"type": "circle",\
             "geometry": [pos_center, diameter],\
             "material": material
             }
    
    def rectangle(self, pos_center, w_x, w_y, material):
        
        # update shapes counter
        self.shapes_counter += 1
     
        # update dict
        self.geo_dict[self.layer_nbr][str(self.shapes_counter)] = \
            {"type": "rectangle",\
             "geometry": [pos_center, w_x, w_y],\
             "material": material
             }
             
    def polygon(self, xy_data, material):
        
        # update shapes counter
        self.shapes_counter += 1
        
        # update dict
        self.geo_dict[self.layer_nbr][str(self.shapes_counter)] = \
            {"type": "polygon",\
             "geometry": xy_data,\
             "material": material
             }
             
                  
        
        