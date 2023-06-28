# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 15:18:46 2021

@author: Dorian
"""

from layer_class import layer_class

class geometry_class():
    
    def __init__(self, period_x, period_y):
        
        # initialize geometry dictionarry
        self.geo_dict = {}
        self.geo_dict["general"] = {"period": [period_x, period_y]}
        
        # intialize layer counter
        self.layer_counter = 0

    def cover(self, material):
        # update Dict
        self.geo_dict["cover"] = material

    def substrate(self, material):
        # update Dict
        self.geo_dict["substrate"] = material

    def layer(self, thickness, surrounding_material):
        # update layer counter
        self.layer_counter += 1
        self.geo_dict[str(self.layer_counter)] = {"thickness": thickness,\
                                                  "surrounding_material": surrounding_material}
        
        return layer_class(thickness, surrounding_material, self.geo_dict)

    def print_g(self):
        
        
        return self.geo_dict

        
if __name__ == "__main__":           
    geo = geometry(10, 100)
    l1 =geo.layer(1, "Al")
    l1.circle([0,0], 10, "Ag")
    l1.circle([0,0], 190, "Si")
    l1.rectangle([0,0], 190, 200, "Si")
    l2 = geo.layer(1, "Al")
    l2.circle([0,0], 10, "Ag")
    l3 = geo.layer(1, "Al")
    l3.polygon([[0,0],[0,1],[1,1],[1,0]], "Ag")
    a = geo.print_g() 
    print(a)      
        
        
        
        