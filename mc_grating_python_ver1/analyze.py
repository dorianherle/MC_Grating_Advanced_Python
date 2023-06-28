# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 17:54:23 2021

@author: Dorian
"""
#from far_field import far_field
import copy

class analyze():
    
    def __init__(self, geometry_dict, 
                 nbr_oders_x=11, nbr_orders_y=11,\
                 polarization_type="Ex", angle_n = 0, angle_p = 0,\
                 single_wavelengh=450):
          
          
          self.g = geometry_dict
          
          # Add analyze dict 
          self.g["analyze"] = {"general": {
                               "nbr_orders_x": nbr_oders_x,\
                               "nbr_orders_y": nbr_orders_y,\
                               "polarization_type": polarization_type,\
                               "angle_n": angle_n,\
                               "angle_p": angle_p,\
                               "single_wavelengh": single_wavelengh\
                              }}
              
    
    def far_field(self,
                   parameter_row = "fixed_parameter", parameter_column="",\
                   start_row=0, end_row=100, nbr_of_points_row=100,\
                   start_column=0, end_column=100, nbr_of_points_column=100,\
                   xy_const=False,\
                   layer_number="",\
                   output_format="power",
                  ): #fixed_wavelength=450
        
        if parameter_row == "":
            raise Exception("Please provide row parameter!")
        
        # remove near-field if exists
        if "near_field" in self.g["analyze"].keys():
             self.g["analyze"].pop("near_field")
         
        self.g["analyze"]["far_field"] = {"parameter_row": parameter_row,\
                                           "parameter_column": parameter_column,\
                                           "start_row": start_row,\
                                           "end_row": end_row,\
                                           "nbr_of_points_row": nbr_of_points_row,
                                           "start_column": start_column,\
                                           "end_column": end_column,\
                                           "nbr_of_points_column": nbr_of_points_column,\
                                           "xy_const": xy_const,\
                                           "layer_number": layer_number,\
                                           "output_format": output_format,\
                                           
                                           } #"fixed_wavelength": fixed_wavelength
     
        return copy.deepcopy(self.g)
             
    
    def near_field(self,\
                         start_row = 0, end_row = 1,\
                         start_column = 0, end_column = 1,\
                         number_of_points = [100,100],\
                         plane = "xy", position_of_plane=0.5,
                         line_position="",
                         output_format = "amplitude_re_im_and_power_flow"):
        
        print("POS: ", position_of_plane)
        
        # if scanning_direction not z -> start/end has to be euqal to 0 to 1
        if "z" not in plane:
            start = 0
            end=1
            
        
        # remove far-field
        if "far_field" in self.g["analyze"].keys():
            self.g["analyze"].pop("far_field")
        
        
        self.g["analyze"]["near_field"] = {
            "scanning_direction_row": plane[0], #scanning_direction_row,
            "scanning_direction_column": plane[1], #scanning_direction_column,
            "start_row": start_row, "end_row": end_row, 
            "start_column": start_column, "end_column": end_column, 
            "number_of_points": number_of_points,
            "plane": plane, "position_of_plane": position_of_plane,
            "line_position": line_position,
            "output_format": output_format
            }
        
        return copy.deepcopy(self.g)
        
    
    # class select_near_field_computation():
        
    #            def __init__(self):
    #                pass
               
    #            def one_d(self,\
    #                      scanning_direction = "z",\
    #                      start = 0, end = 1,\
    #                      plane = "xy", position_of_plane=0.5,
    #                      line_position=""):
                   
                   
    #                print("ONE D")
                    
    #                 pass
                
    #            def two_d(self,\
    #                      scanning_direction = "z",\
    #                      start = 0, end = 1,\
    
    def return_dict(self):
        return self.g.copy()
        

        

        
        