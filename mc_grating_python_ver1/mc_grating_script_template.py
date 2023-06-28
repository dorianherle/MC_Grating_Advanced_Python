# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 18:49:34 2021

@author: Dorian
"""
import numpy as np

class mc_grating_script_template():
    
    def __init__(self):
        pass
        
        
    # -------------- #    
    # -- GEOMETRY -- #
    # -------------- #
    def header(self):
        header_s = \
            "ORDERS_X" + "\t" + "Number of Modes X" + "\n" +\
            "ORDERS_Y" + "\t" + "Number of Modes Y" + "\n" + \
            "PERIOD_X" + "\t" + "X-Period, nm" + "\n"  +\
            "PERIOD_Y" + "\t" + "Y-Period, nm" + "\n" +\
            "SINGLE_W" + "\t" + "Wavelength, nm" + "\n"  +\
            "POLARIZATION_STATE" + "\t" + "Pol State, deg" + "\n" +\
            "POLARIZATION_PHASE" + "\t" + "Pol phase, def" + "\n" +\
            "ANGLE_N" + "\t" + "AngleN, deg (Cover, Plane Base)" + "\n" +\
            "ANGLE_P" + "\t" + "AngleP, deg (Cover, Plane Base)" + "\n" +\
            "NBR_OF_LAYERS" + "\t" + "Number of Layers" + "\n" +\
            "1" + "\t" + "0" + "\t" + "Eps.Re; Eps.Im; Cover; Material:" + " COVER_MATERIAL" + "\n"
        
        return header_s
  
    def layer_header(self):
        layer_header_s = \
            "THICKNESS_LAYER" + "\t" +  "NUMBER_OF_PILLARS" + "\t" + "Layer[" + "LAYER_NUMBER" +"]: Thickness, nm; Number of Pillars" + "\n" \
            "1" + "\t" + "0" + "\t" + "Base" + "\t" + "Eps.Re; Eps.Im; Material:" + " " + "SOURROUNDING_MATERIAL" + "\n" 
        return layer_header_s
    
    def pillar_header(self):
        # Pillar header string
        pillar_header_s = \
            '1' +   '\t' + '0' + '\t' + "NUMBER_OF_POINTS" + '\t' + 'Eps.Re; Eps.Im; Npoints of Pillar[' + "PILLAR_NUMBER" + ']; Material:' + ' ' + "PILLAR_MATERIAL" + "\n" 
        return pillar_header_s
    
    def points(self):
        # Points string
        points_s = \
            "POINT_X" + "\t" + "POINT_Y" + "\t" + "x[" + "POINT_NUMBER" + "]/(X-Period); y[" + "POINT_NUMBER" + "]/(Y-Period)" + "\n" 
        return points_s
    
    def substrate(self):
        # Substrate string
        substrate_s = "1" + "\t" + "0" + "\t" + "Eps.Re; Eps.Im; Substrate; Material:" + " " + "SUBSTRATE_MATERIAL" + "\n"
        return substrate_s
    
    # -------------- #    
    # -- SCANNING -- #
    # -------------- #        
    def analyze_header(self):
        
        analyze_header  = \
        """
 Settings For FMM Crossed
true                UsingNF
false               UsingSmooth
0                   SmoothStart
0                   SmoothWidth
true                Background Calculation
false               Memory Consumption
false               Calculation Direction From Substrate to Cover
Y                   Number of Modes Correction
3                   Calculation Priority
true                Cover Diffraction Orders
true                Substrate Diffraction Orders
ORDER_C_X                   ORDER_C_Y                   X,Y Cover Order
ORDER_S_X                   ORDER_S_Y                   X,Y Substrate Order    
-1000                   10000                   MinMax X Cover Order
-1000                   10000                   MinMax Y Cover Order
-1000                   10000                   MinMax X Substrate Order
-1000                   10000                   MinMax Y Substrate Order
SCANNING_OUTPUT_FORMAT                   Scanning Output Format: Power; Amplitude (Re,Im); Amplitude (Mod,Ph)
true                Diffraction Angle
%10.5g              Decimal Digits Format
STRING_OF_ROW_SCANNING_PARAMETER               String of Scanning Parameter
ROW_SCANNING_PARAMETER                   Row Scanning Parameter
NUMBER_OF_ROW_POINTS                  Nubmber of Column Points -1"""
        return analyze_header
    
    def get_MC_Grating_Scaning_Output_Format(self, output_format):
        
        if output_format == "power":
            return "P"
        
        if output_format == "amplitude_re_im":
            return "A"
        
        if output_format == "amplitude_mod_ph":
            return "M"
        
        if output_format == "amplitude_re_im_and_power":
            return "a"
        
        if output_format == "amplitude_mod_ph_and_power":
            return "m"

    
    def get_MC_grating_scanning_parameter_and_string(self,scan_type, 
                                                     constant_xy_ratio="False",
                                                     layer_number=""):
        
        if scan_type == "fixed_parameter":
            return [0,"Fixed Parameters"]
        
        if scan_type == "angle_n":
            return [1,"AngN"]
        
        if scan_type == "angle_p":
            return [2,"AngP"]
        
        if scan_type == "period_x":
            if constant_xy_ratio:
                return [3,"X-per(Y)"]
            else:
                return [3,"X-per"]
            
        if scan_type == "period_y":
            if constant_xy_ratio:
                return [4,"Y-per(X)"]
            else:
                return [4,"Y-per"]
       
        if scan_type == "wavelength":
            return [5,"WaveL"]
        
        if scan_type == "layer":
            if layer_number == "":
                raise Exception("Please provide the layer number.")
            else:
                return [6,"Thick["+str(layer_number)+"]"]
            
    
    def constant_xy_ratio(self):
        return "NUMBER_OF_COLUMN_POINTS" + "\t" + "X-per/Y-per=const"
    
    def scanning_ranges(self):
        
        scanning_ranges = np.array([
            [0,0,"Range of Scanning AngleN"],
            [0,0,"Range of Scanning AngleP"],
            [0,0,"Range of Scanning X-Period"],
            [0,0,"Range of Scanning Y-Period"],
            [0,0,"Range of Scanning Wavelength"],
            [0,0,"Range of Scanning RowVar"],
            [0,0,"Range of Scanning ColVar"]
            ], dtype=object)
        return scanning_ranges
        
    def scanning_ranges_thickness_layer(self):
        return  "\n"+"LAYER_START" + '\t' + "LAYER_END" + '\t' + "Range of Scanning Thickness of Layer[LAYER_NUMBER]"
    
    
    def twoD_scan(self):
        
        twoD_scan = """
COLUMN_SCANNING_PARAMETER                   Col Scanning Parameter
0                   2D fmt 
NUMBER_OF_COL_POINTS                  Number of Col Points -1 
1STRING_1OF_1COLUMN_1SCANNING_PARAMETER                String of Col Scanning Parameter"""
        
        return twoD_scan
    
    # -------------- #    
    # --   FIELD  -- #
    # -------------- # 
    
    def field_scan(self):
        
        settings_bottom = np.array([
            #["1", "Layer Index"],
            ["true","Advanced Options"],
            ["false","Field Calculation"],
            ["true","Automatically Starts Second Stage"],
            [100,"Nubmber of Points -1 for NPFx"],
            [100,"Nubmber of Points -1 for NPFy"],
            [100,"Nubmber of Points -1 for NPFz"],
            [0,"2D Field Component"],
            [0,"Field Output Format: Power Flow; Amplitude (Re,Im); Amplitude (Mod,Ph)"],
            ["%10.5g","Decimal Digits Format"],
            ["1","Scanning Direction"],
            ["0"," Units Along X"],
            ["0","1","Range of Scanning X Direction"],
            ["0","1","Range of Scanning Y Direction"],
            ["0","1","Range of Scanning Z Direction"],
            ["0","Fixed X"],
            ["0","Fixed Y"],
            ["0","Fixed Z"],
            ["0","Graph Index for X-axis"],
            ["1","Graph Index for Y-axis"],
            ["2","Graph Index for Z-axis"],
            ["0","Polarization Type (Any, Hy=0, Ey=0, Hx=0, Ex=0)"],
            ["false","Output Polarizer"],
            ["0","Polarizer Angle, deg"]
            ], dtype=object)
        
        return settings_bottom
    
    
    
    def input_polarization(self, polarization_type="Any"):
        
        # Seeting scan type variabel 
        if polarization_type == 'Any' : polarization_type = '0'
        if polarization_type == 'Hy': polarization_type = '1'
        if polarization_type == 'Ey': polarization_type = '2'
        if polarization_type == 'Hx': polarization_type = '3'
        if polarization_type == 'Ex': polarization_type = '4'
        
        return polarization_type
    
    
    def getting_cross_section_plane(self,cross_section_plane):
        
        # getting scan direction
        if cross_section_plane == "zx": scan_dir_ = "11"
        if cross_section_plane == "zy": scan_dir_ = "1"
        if cross_section_plane == "yz": scan_dir_ = "0"
        if cross_section_plane == "xz": scan_dir_ = "10"
        if cross_section_plane == "yx": scan_dir_ = "21"
        if cross_section_plane == "xy": scan_dir_ = "20"
        
        return scan_dir_
    
    
    def getting_field_output_format(self, output_format):
        
        if output_format == "power_flow":
            return "P"
        
        if output_format == "amplitude_re_im":
            return "A"
        
        if output_format == "amplitude_re_im_and_power_flow":
            return "a"
        
        if output_format == "amplitude_mod_ph":
            return "M"
        
        if output_format == "amplitude_mod_ph_and_power_flow":
            return "m"
        
        
        
        
        

#  """
# 1                   Layer Index
# ADVANCED_OPTICS               Advanced Options
# FIELD_CALCULATION               Field Calculation
# true                Automatically Starts Second Stage
# NBR_OF_POINTS_2D_X                 Nubmber of Points -1 for NPFx
# NBR_OF_POINTS_2D_Y                 Nubmber of Points -1 for NPFy
# NBR_OF_POINTS_2D_Z                 Nubmber of Points -1 for NPFz
# 2D_FIELD_COMPONENT                  2D Field Component
# m                   Field Output Format: Power Flow; Amplitude (Re,Im); Amplitude (Mod,Ph)
# %10.5g              Decimal Digits Format
# SCANNING_DIRECTION                   Scanning Direction
# 0                   Units Along X
# FIELD_X_MIN                   FIELD_X_MAX                   Range of Scanning X Direction
# FIELD_Y_MIN                   FIELD_Y_MAX                   Range of Scanning Y Direction
# FIELD_Z_MIN                   FIELD_Z_MAX                  Range of Scanning Z Direction
# FIELD_X_POS                   Fixed X
# FIELD_Y_POS                   Fixed Y
# FIELD_Z_POS                   Fixed Z
# 0                   Graph Index for X-axis
# 1                   Graph Index for Y-axis
# 2                   Graph Index for Z-axis
# POLARIZATION_TYPE                   Polarization Type (Any, Hy=0, Ey=0, Hx=0, Ex=0)
# false               Output Polarizer
# 0                   Polarizer Angle, deg
# """
        
    
    
    



