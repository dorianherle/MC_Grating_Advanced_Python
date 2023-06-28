# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 18:48:42 2021

@author: Dorian
"""
import numpy as np

from mc_grating_script_template import mc_grating_script_template

class dict_to_mc_grating_script():
      
    def __init__(self, dict_):
          self.dict_ = dict_
          self.mc = mc_grating_script_template()
          self.polarization_type = 0
          
    
    def geometry_type_to_points(self,dict_geometry, periodx, periody):
        """
        Outputs the points of a geometry defined by a certain type, i.e.
        "circle", "rectangle", "polygon"
    
        Parameters
        ----------
        dict_geometry : dict
            layer dictionary 
            
        periodx : float
            period in x-direction 
    
        periody : float
             period in y-direction 
             
        Returns
        -------
        points : 2D
            Nx2 list of points (coordinates, normalized to period).
    
        """
        
        # Extract type
        type_ = dict_geometry["type"]
        # Extract geometric dimensions
        dim = dict_geometry["geometry"]
        # convert those to numpy
        dim = [np.array(a) for a in dim]
        
        
        # circle
        if type_ == "circle":
            center_p, diameter = dim
            # Put center into middle (MC Grating Convention)
            center_p = center_p + np.array([periodx/2, periody/2])
            # 2 Points for cricle
            #   -----2
            # |      |
            # |      |
            # 1 ----- 
            points = [[center_p[0]-diameter/2,center_p[1]+diameter/2],\
                      [center_p[0]+diameter/2,center_p[1]-diameter/2]]
            
        if type_ == "rectangle":
            center_p, wx, wy = dim
            # Put center into middle (MC Grating Convention)
            center_p = center_p + np.array([periodx/2, periody/2])
            # 4 Points for rectangle
            # Anti-clockwise
            # 4 ----- 3
            # |       |
            # |       |
            # 1 ----- 2
            points = [[center_p[0]-wx/2, center_p[1]-wy/2],\
                      [center_p[0]+wx/2, center_p[1]-wy/2],\
                      [center_p[0]+wx/2, center_p[1]+wy/2],\
                      [center_p[0]-wx/2, center_p[1]+wy/2],\
                    ]
            
        if type_ == "polygon":
            points = dim
        
        
        # Normalize Coordiantes to Periods 
        points_norm = [[p[0]/periodx, p[1]/periody] for p in points]
        
        return points_norm


    def gui_input_structure(self, extract_geometry=False):
        
        # Initialize
        dict_ = self.dict_
        mc = self.mc
        
        # EXTRACT NBR of Layer
        
        layers = [k for k in dict_.keys() if k.isnumeric()]
        nbr_of_layer = len(layers)
        
        
        # HEADER
    
        # Extract general information
        periodx, periody = dict_["general"]["period"]
        nbr_orders_x = dict_["analyze"]["general"]["nbr_orders_x"]
        nbr_orders_y = dict_["analyze"]["general"]["nbr_orders_y"]
        polarization_type = dict_["analyze"]["general"]["polarization_type"]
        angle_n = dict_["analyze"]["general"]["angle_n"]
        angle_p = dict_["analyze"]["general"]["angle_p"]
        single_wavelengh = dict_["analyze"]["general"]["single_wavelengh"]
        
        # try: 
        #     single_wavelengh = dict_["analyze"]["far_field"]["fixed_wavelength"]
        # except:
        #     pass
        
        cover_material = dict_["cover"]
    
        # Order and POLARIZATION_TYPE
        if polarization_type == 'Any' : polarization_type = '0'
        if polarization_type == 'Hy': polarization_type = '1'
        if polarization_type == 'Ey': polarization_type = '2'
        if polarization_type == 'Hx': polarization_type = '3'
        if polarization_type == 'Ex': polarization_type = '4'
        
        self.polarization_type = polarization_type
            
        
        # get header template
        header = mc.header()
    
        # fill in template
        header = header.replace("ORDERS_X", str(nbr_orders_x))
        header = header.replace("ORDERS_Y", str(nbr_orders_y))
        header = header.replace("PERIOD_X", str(periodx))
        header = header.replace("PERIOD_Y", str(periody))
        header = header.replace("SINGLE_W", str(single_wavelengh))
        header = header.replace("POLARIZATION_STATE", str(0)) # WILL BE OVERWRITTEN BY THE POLYRIZATION TYPE
        header = header.replace("POLARIZATION_PHASE", str(0))
        header = header.replace("ANGLE_N", str(angle_n))
        header = header.replace("ANGLE_P", str(angle_p))
        header = header.replace("NBR_OF_LAYERS", str(nbr_of_layer))
        header = header.replace("COVER_MATERIAL", str(cover_material))
        
        # Populate layers
        layers_text = ""
        layers_dict = {"cover_material": cover_material}
        for layer in layers:
            #  Add layer key to layer dict
            layers_dict[layer] = {}
            
            # Get layer header template
            layer_header = mc.layer_header()
            
            # Extract number of pillars
            pillars = [k for k in dict_[layer].keys() if k.isnumeric()]
            nbr_of_pillars = len(pillars)
            
            # Get layer thickness + surrounding material
            thickness = dict_[layer]["thickness"]
            surrounding_material = dict_[layer]["surrounding_material"]
            
            #  Add layer thickness, and surrounding_material
            layers_dict[layer] = {"thickness": thickness, 
                                  "surrounding_material": surrounding_material}

            # Fill in template
            layer_header = layer_header.replace("LAYER_NUMBER", str(layer))
            layer_header = layer_header.replace("NUMBER_OF_PILLARS", str(nbr_of_pillars))
            layer_header = layer_header.replace("THICKNESS_LAYER", str(thickness))
            layer_header = layer_header.replace("SOURROUNDING_MATERIAL", str(surrounding_material))
        
            # Populate pillar
            pillars_text = ""
            for pillar in pillars:
                
                #  Add pillar key to layer dict
                layers_dict[layer][pillar]={}
                
                # Get pillar header template
                pillar_header = mc.pillar_header()
                 
                # Extract pillar material
                pillar_material = dict_[layer][pillar]["material"]
                
                # Points (Coordinates, normlaized to periods) 
                points = self.geometry_type_to_points(dict_[layer][pillar], 
                                                 periodx, periody)
                
                #  Add materials and points
                layers_dict[layer][pillar]={"material": pillar_material, "points": points}
                
                # Extract number of points
                number_of_points = len(points)
                
                # Fill in template
                pillar_header = pillar_header.replace("PILLAR_NUMBER", str(pillar))
                pillar_header = pillar_header.replace("PILLAR_MATERIAL", str(pillar_material))
                pillar_header = pillar_header.replace("NUMBER_OF_POINTS", str(number_of_points))
                
                points_text = ""
                for i, point in enumerate(points):
                    
                    # Get points string template
                    points_s = mc.points()
                    
                    # Point number
                    point_number = i+1
    
                    # Fill in template
                    points_s = points_s.replace("POINT_X", str(point[0]))
                    points_s = points_s.replace("POINT_Y", str(point[1]))
                    points_s = points_s.replace("POINT_NUMBER", str(point_number))
                    
                    # Add to points_text
                    points_text += points_s
                
                # Put it together -> Pillar text
                pillar_text = pillar_header + points_text
                # Add to pillars text
                pillars_text += pillar_text
            
            # Put it together
            layers_text += layer_header + pillars_text
            
            # Populate Substrate
            # Get Substrate template
            substrate_text = mc.substrate()
            
            # Get substrate material
            substrate_material = dict_["substrate"]
            layers_dict["substrate_material"]=substrate_material
            # Fill in the template
            substrate_text = substrate_text.replace("SUBSTRATE_MATERIAL", substrate_material)
            
            
        # PUT IT ALL TOGETHER
        if extract_geometry:
            return header+layers_text+substrate_text, layers_dict
        else:
            return header+layers_text+substrate_text
            
    def analyze_section(self):
        
        # Initialize
        dict_ = self.dict_
        mc = self.mc
        
        # Get analyze header template
        analyze_header = mc.analyze_header()

        
        # Fill in template
        analyze_header = analyze_header.replace("ORDER_C_X", "0") #In any case all orders are exported
        analyze_header = analyze_header.replace("ORDER_C_Y", "0")
        analyze_header = analyze_header.replace("ORDER_S_X", "0")
        analyze_header = analyze_header.replace("ORDER_S_Y", "0")
        
        # ------------------- #    
        # -- GET SCAN TYPE -- #
        # ------------------- #   
        
        scan_type = list(dict_["analyze"].keys())[1]
        print("SCAN TYPE: ", scan_type)
        
        if scan_type == 'far_field':
            # get field scan settings -> unimportant for far field
            field_scan = "\n".join(["\t".join([str(e) for e in elem]) for elem in mc.field_scan()])
            return self.far_field(analyze_header)+"\n"+field_scan
        
        if scan_type == 'near_field':
            return self.near_field(analyze_header)
        
        
    def far_field(self, analyze_header):
  
        # Initialize
        dict_ = self.dict_["analyze"]["far_field"]
        mc = self.mc
        far_field_text = ""
        
        # Extract all parameters
        parameter_row = dict_["parameter_row"]
        parameter_column = dict_["parameter_column"]
        start_row = dict_["start_row"]
        end_row = dict_["end_row"]
        nbr_of_points_row = dict_["nbr_of_points_row"]
        start_column = dict_["start_column"]
        end_column = dict_["end_column"]
        nbr_of_points_column = dict_["nbr_of_points_column"]
        xy_const = dict_["xy_const"]
        layer_number = dict_["layer_number"]
        output_format = dict_["output_format"]
        # fixed_wavelength = dict_["fixed_wavelength"]
        
        
        
        
        # --> scanning output
        # Convert to MC Grating format
        scanning_output_format = mc.get_MC_Grating_Scaning_Output_Format(output_format)
        # Fill the template
        analyze_header = analyze_header.replace("SCANNING_OUTPUT_FORMAT", scanning_output_format)
        
        # --> scanning parameter
        number_param, strings_of_scanning_parameter =\
            mc.get_MC_grating_scanning_parameter_and_string(\
                                                            parameter_row,\
                                                            constant_xy_ratio=xy_const,\
                                                            layer_number=layer_number)
      
        # Fill the template
        analyze_header = analyze_header.replace("STRING_OF_ROW_SCANNING_PARAMETER", str(strings_of_scanning_parameter))
        analyze_header = analyze_header.replace("ROW_SCANNING_PARAMETER", str(number_param))
        
        # Fill the template
        analyze_header = analyze_header.replace("NUMBER_OF_ROW_POINTS", str(nbr_of_points_row))
        
        # Add to text
        far_field_text += analyze_header
        # Check if single parameter scan
            
        # Check if xy == constant ratio
        if xy_const and parameter_column == "":
            # Extract xy_constant template
            constant_xy_ratio = mc.constant_xy_ratio()
            # Fill in template
            constant_xy_ratio = constant_xy_ratio.replace("NUMBER_OF_COLUMN_POINTS", str(nbr_of_points_row))
            # Add to text
            far_field_text += constant_xy_ratio
            
        # Check if 2D Scan
        if parameter_column != "":

            # --> scanning parameter
            number_param, strings_of_scanning_parameter =\
                mc.get_MC_grating_scanning_parameter_and_string(\
                                                                parameter_column,\
                                                                constant_xy_ratio=xy_const,\
                                                                layer_number=layer_number)
            
            # Get number of points to sample
            nbr_of_points = dict_["nbr_of_points_column"]
          
            # Extract template
            twoD_scan = mc.twoD_scan()
            # Fill in template
            twoD_scan = twoD_scan.replace("COLUMN_SCANNING_PARAMETER", str(number_param))
            twoD_scan = twoD_scan.replace("NUMBER_OF_COL_POINTS", str(nbr_of_points))
            twoD_scan = twoD_scan.replace("1STRING_1OF_1COLUMN_1SCANNING_PARAMETER ", str(strings_of_scanning_parameter)) #TODO replace not perfect
            
            # Add to text
            far_field_text += twoD_scan
        
        # Add Scans Template
        scanning_ranges = mc.scanning_ranges()
 
        # Fill Template
        parameters = ["angle_n", "angle_p", "period_x", "period_y", "wavelength", "layer"]
        
        # Fix for "Fixed parameter" -> Once fixed parameter set, program does not care about the ranges
        if parameter_row == "fixed_parameter": 
            index_row = 0
        else:
            index_row = parameters.index(parameter_row)
        try:
            index_colum = parameters.index(parameter_column)
        except: 
            index_colum = -1
        
        # Assign Start of Row Parameter
        scanning_ranges[index_row][0]=str(start_row)
        # Assign End of Row Parameter
        scanning_ranges[index_row][1]=str(end_row)
        
        if index_colum >= 0:
            # Assign Start of Column Parameter
            scanning_ranges[index_colum][0]=str(start_column)
            # Assign End of Column Parameter
            scanning_ranges[index_colum][1]=str(end_column)
            
        # Convert table to MC Grating text
        scanning_ranges = "\n"+"\n".join(["\t".join([str(e) for e in elem]) for elem in scanning_ranges])
        # Add to text:
        far_field_text += scanning_ranges
        
        
        # Add layers
        # get layers:
        layers = {k:v for k,v in self.dict_.items() if k.isnumeric()}
        
        # ADD LAYER SCAN
        for layer_nbr,layer_content in layers.items():
            
            # get layer template
            layer_template = mc.scanning_ranges_thickness_layer()
            layer_template = layer_template.replace("LAYER_NUMBER", str(layer_nbr)) 
            
            # check if layer scan
            if layer_number == layer_nbr:
                
                if parameter_column == "layer":
                    layer_template = layer_template.replace("LAYER_START", str(start_column))
                    layer_template = layer_template.replace("LAYER_END", str(end_column))
                
                if parameter_row == "layer":
                    layer_template = layer_template.replace("LAYER_START", str(start_row))
                    layer_template = layer_template.replace("LAYER_END", str(end_row))
                
            else:
                # get thickness
                layer_thickness = layer_content["thickness"]
                layer_template = layer_template.replace("LAYER_START", str(layer_thickness))
                layer_template = layer_template.replace("LAYER_END", str(layer_thickness))
                
            far_field_text+= layer_template
            
        # Add layer index -> for layer scan
        far_field_text+= "\nLAYER_INDEX                   Layer Index".replace("LAYER_INDEX", str(layer_nbr))
       
   
        return far_field_text
    
    def get_index(self, string, field_scan_template):
        
        return [i for i in range(len(field_scan_template)) if  string in field_scan_template[i][-1]][0]
    
    def near_field(self, analyze_header):
        
        # Initialize
        dict_ = self.dict_["analyze"]["near_field"]
        mc = self.mc
        far_field_text = ""
        
        # Fill the template
        analyze_header = analyze_header.replace("SCANNING_OUTPUT_FORMAT", "P")
        
        # --> scanning parameter
        number_param, strings_of_scanning_parameter =\
            mc.get_MC_grating_scanning_parameter_and_string(\
                                                            "fixed_parameter",\
                                                            constant_xy_ratio=False,\
                                                            layer_number="")
        # Fill the template
        analyze_header = analyze_header.replace("STRING_OF_ROW_SCANNING_PARAMETER", str(strings_of_scanning_parameter))
        analyze_header = analyze_header.replace("ROW_SCANNING_PARAMETER", str(number_param))
        
        # Fill the template
        analyze_header = analyze_header.replace("NUMBER_OF_ROW_POINTS", str(0))
        
        
        # Extract the parameters
        scanning_direction_row = dict_["scanning_direction_row"]
        scanning_direction_column = dict_["scanning_direction_column"]
        start_row = dict_["start_row"]
        end_row = dict_["end_row"]
        start_column = dict_["start_column"]
        end_column = dict_["end_column"]
        nbr_of_points_row_dir, nbr_of_points_column_dir = dict_["number_of_points"]
        plane = dict_["plane"]
        position_of_plane = dict_["position_of_plane"]
        line_position = dict_["line_position"]
        output_format = dict_["output_format"]
        
        # get coordiante of other scan direction (plane)
        scanning_direction_column = plane.replace(scanning_direction_row, "")
        
        # Get field scan template
        field_scan_template = mc.field_scan()
        
        # Initialize field calculation
        index_field_calculation = self.get_index("Field Calculation", field_scan_template)
        field_scan_template[index_field_calculation][0] = "true" #"FIELD_CALCULATION"
        
        
        
        # Enter number of points to compute -> ROW
        index_scan_dir = self.get_index("NPF"+scanning_direction_row, field_scan_template)
        field_scan_template[index_scan_dir][0] = nbr_of_points_row_dir # Nubmber of Points -1 for NPF...
        
        # if a 2D scan -> line_position not defined ("")
        if line_position == "":
            index_scanning_direction_column = self.get_index("NPF"+scanning_direction_column, field_scan_template)
            field_scan_template[index_scanning_direction_column][0] = nbr_of_points_column_dir
        
        # if a 1D scan -> line position defined -> NOT 2D FIELD SCAN but 1D
        if line_position != "":
            index_scanning_direction_column = self.get_index("2D Field Component", field_scan_template)
            field_scan_template[index_scanning_direction_column][0] = -1
            
        # Get field output format
        index_output = self.get_index("Field Output Format", field_scan_template) 
        field_scan_template[index_output][0] = mc.getting_field_output_format(output_format)
        
        # Get scanning direction (ROW)
        index_scanning_dir = self.get_index("Scanning Direction", field_scan_template) 
        field_scan_template[index_scanning_dir][0] = mc.getting_cross_section_plane(plane)
        
        # Get Scanning Range
        index_scanning_range = self.get_index("Range of Scanning " + scanning_direction_row.upper() +" Direction", field_scan_template)
        field_scan_template[index_scanning_range][0] = start_row
        field_scan_template[index_scanning_range][1] = end_row
        
         # if a 2D scan -> line_position not defined ("")
        if line_position == "":
            # Get Scanning Range
            index_scanning_range = self.get_index("Range of Scanning " + scanning_direction_column.upper() +" Direction", field_scan_template)
            field_scan_template[index_scanning_range][0] = start_column
            field_scan_template[index_scanning_range][1] = end_column
            
        # Position of plane
        coord = ("xyz".replace(scanning_direction_row, "")).replace(scanning_direction_column, "")
        index_position_of_plane = self.get_index("Fixed " + coord.upper(), field_scan_template)
        field_scan_template[index_position_of_plane][0] = position_of_plane
        
        # Position of line on second coordinate, if want to extract values on a line
        if line_position != "":
            coord = ("xyz".replace(scanning_direction_row, "")).replace(coord, "")
            index_position_of_plane = self.get_index("Fixed " + coord.upper(), field_scan_template)
            field_scan_template[index_position_of_plane][0] = line_position
        
        
        ranges = "\n".join(["\t".join([str(e) for e in elem]) for elem in mc.scanning_ranges()])
        
        # Add layers
        # get layers:
        layers = {k:v for k,v in self.dict_.items() if k.isnumeric()}
        
        # ADD LAYER SCAN
        for layer_nbr,layer_content in layers.items():
            
            # get layer template
            layer_template = mc.scanning_ranges_thickness_layer()
            layer_template = layer_template.replace("LAYER_NUMBER", str(layer_nbr)) 
            
           
            # get thickness
            layer_thickness = layer_content["thickness"]
            layer_template = layer_template.replace("LAYER_START", str(layer_thickness))
            layer_template = layer_template.replace("LAYER_END", str(layer_thickness))
                
            ranges += layer_template
        
        # Add layer index -> Default 
        ranges+= "\n1                   Layer Index"
        
        # ADD Polarization Type
        index_field_calculation = self.get_index("Polarization Type (Any, Hy=0, Ey=0, Hx=0, Ex=0)", field_scan_template)
        field_scan_template[index_field_calculation][0] = self.polarization_type #POLARIZATION TYPE
        print("POL: ",  self.polarization_type)
        

            
        
        return analyze_header+"\n"+ranges+"\n"+"\n".join(["\t".join([str(e) for e in elem]) for elem in field_scan_template])
        
        
        
    
        
        
        
        
        
        
                
            
            

if __name__ == "__main__":
    dict_ = {'general': {'period': [137, 137]}, 'cover': 'Air (Special Formula)', '1': {'thickness': 10, 'surrounding_material': 'Silicon (Table)'}, 'substrate': 'Fused Silica (Sellmeier)', 'analyze': {'general': {'nbr_orders_x': 11, 'nbr_orders_y': 11, 'polarization_type': 'Ex', 'angle_n': 0, 'angle_p': 0, 'single_wavelengh': 450}, 'near_field': {'scanning_direction_row': 'y', 'scanning_direction_column': 'z', 'start_row': 1, 'end_row': 0, 'start_column': -100, 'end_column': 110, 'number_of_points': [10, 10], 'plane': 'yz', 'position_of_plane': 0.5, 'line_position': '', 'output_format': 'amplitude_re_im_and_power_flow'}}}
    # print(dict_.keys())
    # print(dict_["general"]["period"])
    # print(dict_["analyze"]["general"])
    # print(dict_["1"])
    a = dict_to_mc_grating_script(dict_)
    print(a.gui_input_structure())
    print(a.analyze_section())
