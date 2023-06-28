# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 14:40:08 2021

@author: Dorian
"""
import subprocess
import json
import numpy as np
import re

def is_number(string):
      try:
          float(string)
          return True
      except ValueError:
          return False
    
class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

    
def launch_simulation(document, path="mc_grating_simulation_file", fix=False):
    
    
    # Export Script
    export_doument(document,path)
    

    # Run the file as command line input
    subprocess.run(['C:\Program Files\MC Grating Software\Full\ModalCrossed.exe', path+'.mdl',  path+'_output.json'])
    
    # Extract result
    result = open_output(path, fix=fix)
    
    # get periods
    periodx = float((document.split("X-Period, nm")[0]).split("\n")[-1])
    periody = float((document.split("Y-Period, nm")[0]).split("\n")[-1])
        
    # add periods
    result["Period_X"] = periodx
    result["Period_Y"] = periody 
    
    if result["SimulationType"] == "Fields":
     
        def from_scan_dir_to_cross_section(scan_dir):
        
            # getting scan direction
            if scan_dir == "11": cross_section_plane = "zx"
            if scan_dir ==  "1": cross_section_plane = "zy"
            if scan_dir ==  "0": cross_section_plane = "yz"
            if scan_dir == "10": cross_section_plane = "xz"
            if scan_dir == "21": cross_section_plane = "yx"
            if scan_dir == "20": cross_section_plane = "xy"
            
            return cross_section_plane
        
        # ADD CROSS-SECTION NORMAL 
        # get normal
        scan_dir = str(int((document.split("Scanning Direction")[0]).split("\n")[-1]))
        cross_section_plane = from_scan_dir_to_cross_section(scan_dir)
    
        normal = "xyz"
        normal = normal.replace(cross_section_plane[0], "")
        normal = normal.replace(cross_section_plane[1], "")
        normal = normal.upper()
        result["Normal_of_CrossSection"] = normal
        
        # ADD POSITION OF CROSS-SECTION
        fixed_pos = (document.split("Fixed "+normal)[0]).split("\n")[-1]

        if normal == "X":
            result["Pos_of_CrossSection"] = float(fixed_pos)#*periodx-periodx/2
            result["Pos_of_CrossSection_Value"] = periodx*(float(fixed_pos)-1/2)
        if normal == "Y":
            result["Pos_of_CrossSection"] = float(fixed_pos)#*periody-periody/2
            result["Pos_of_CrossSection_Value"] = periody*(float(fixed_pos)-1/2)
        if normal == "Z":
            result["Pos_of_CrossSection"] = float(fixed_pos)
            result["Pos_of_CrossSection_Value"] = float(fixed_pos)
            

        print("result[Pos_of_CrossSection]", result["Pos_of_CrossSection"])
        
        # ADD WAVELENGTH
        result["WaveL"] = str(float((document.split("Wavelength, nm")[0]).split("\n")[-1]))
    
        
        # get SCANNING PARAMETER ROW
        
        scan_param_row =  result["ScanPar"]["Name"][0] # name
        scan_param_row_nbr = result["ScanPar"]["NPoints"] # nbr of points
        scan_param_row_min = result["ScanPar"]["IniValue"] # initial value
        scan_param_row_max = result["ScanPar"]["FinValue"] # final value
        
       
        if scan_param_row == "X":
            start = periodx*scan_param_row_min-periodx/2
            end = periodx*scan_param_row_max-periodx/2
            row_points = np.linspace(start, end, scan_param_row_nbr)
            
            
        if scan_param_row == "Y":
            start = periody*scan_param_row_min-periody/2
            end = periody*scan_param_row_max-periody/2
            row_points = np.linspace(start, end, scan_param_row_nbr)

            
        if scan_param_row != "Z":
            # Reasign keys
            param_row = [k for k in result.keys() if is_number(k)]
            for i, old_key in enumerate(param_row):
                result[str(row_points[i])] = result.pop(old_key)
                
            result["ScanPar"]["Name"] = scan_param_row    
            result["ScanPar"]["IniValue"] = start
            result["ScanPar"]["FinValue"] = end
            
        
        
        if result["ScanType"] == 2: # 1D Line 
            # get SCANNING PARAMETER COLUMN
            scan_param_c =  result["ScanParC"]["Name"][0] # name
            scan_param_c_nbr = result["ScanParC"]["NPoints"] # nbr of points
            scan_param_c_min = result["ScanParC"]["IniValue"] # initial value
            scan_param_c_max = result["ScanParC"]["FinValue"] # final value
            
            
            if scan_param_c == "X":
                start = periodx*scan_param_c_min-periodx/2
                end = periodx*scan_param_c_max-periodx/2
                c_points = np.linspace(start, end, scan_param_c_nbr)
                
            if scan_param_c == "Y":
                start = periody*scan_param_c_min-periody/2
                end = periody*scan_param_c_max-periody/2
                c_points = np.linspace(start, end, scan_param_c_nbr)
                
            # if scan_param_c == "Z":
            #     start = periody*scan_param_c_min-periody/2
            #     end = periody*scan_param_c_max-periody/2
            #     c_points = np.linspace(start, end, scan_param_c_nbr)
               
            if scan_param_c != "Z":
                # Reasign keys
                param_row = [k for k in result.keys() if is_number(k)]
                for param_r in param_row:
                    subkeys = list(result[param_r].keys()).copy()
                    for i,param_c in enumerate(subkeys):
                        result[param_r][str(c_points[i])] = result[param_r].pop(param_c)
                    
                result["ScanParC"]["Name"] = scan_param_c
                result["ScanParC"]["IniValue"] = start
                result["ScanParC"]["FinValue"] = end
                    
            
        
        if result["ScanType"] == 1: # 1D Line 
             y_label = result["ScanPar"]["Name"][0]
             x_label = cross_section_plane.replace(y_label.lower(), "").upper()
             
             result["Pos_of_Line"] = str(float((document.split("Fixed "+x_label)[0]).split("\n")[-1]))
           
            
        
    else:
        # Add total power
        result = add_total_power(result)
        
    # Overwriten original document
    dumped = json.dumps(result, cls=NumpyEncoder)
    with open(path+'_output.json', 'w', encoding='utf-8') as f:
        json.dump(dumped, f, ensure_ascii=False, indent=4)
        
    return result
    
def export_doument(document,path):
    
    with open(path+ ".mdl", 'w') as new_file: 
        new_file.write(document)
    new_file.close()
    

def open_output(path, fix = False):
    
    if fix: 
        # TEMPORARY FIX FOR ISSUE WITH FIELD LINE
        f = open(path+'_output.json', "r"); t = f.read()
        
        def replacer(s, newstring, index, nofail=False):
            # raise an error if index is outside of the string
            if not nofail and index not in range(len(s)):
                raise ValueError("index outside given string")
        
            # if not erroring, but the index is still not in the correct range..
            if index < 0:  # add it to the beginning
                return newstring + s
            if index > len(s):  # add it to the end
                return s + newstring
        
            # insert the new string between "slices" of the original
            return s[:index] + newstring + s[index + 1:]
    
        
        t_new = t
        offset = 0
        for i,m in enumerate([m.start() for m in re.finditer("\"0\"", t)]): 
            m = offset+m
            t_new = replacer(t_new, str(i), m+1)
            if len(str(i)) != 1:
                offset += len(str(i))-1
         
        with open(path+'_output.json', 'w') as f:
            json.dump(eval(t_new), f)
    
        # END TEMPORARY FIX 
    
    # Get MC Grating Result
    with open(path+'_output'+'.json') as json_file:
        output = json.load(json_file)
    
    return output
    
##### ADD POWER ############

def add_total_power(mc_grating_json):
    scan_type = mc_grating_json["ScanType"]
    
    # print(scan_type)
    if scan_type == 0:
        mc_grating_json = add_total_power_0D(mc_grating_json)
    if scan_type == 1:
        mc_grating_json = add_total_power_1D(mc_grating_json)
    if scan_type == 2:
        mc_grating_json = add_total_power_2D(mc_grating_json)
       
    return mc_grating_json
       

def add_total_power_0D(mc_grating_json):
    
    r_value = "Fixed"
    total_r, total_t = total_power(mc_grating_json[r_value])
    absorption = 1-mc_grating_json[r_value]["Balance"]
    # Add to JSON File
    mc_grating_json[r_value]["Cover orders"]["Total Power"] = total_r
    try: # if exists (substrate transparent)
        mc_grating_json[r_value]["Substrate orders"]["Total Power"] = total_t
    except:
        pass
    mc_grating_json[r_value]["Absorption"] = absorption
    
    mc_grating_json[r_value]["Cover orders"]["Power"] = list(\
            np.array(mc_grating_json[r_value]["Cover orders"]["Power (s)"]) + 
            np.array(mc_grating_json[r_value]["Cover orders"]["Power (p)"]))
        
    try: # if exists (substrate transparent)
        mc_grating_json[r_value]["Substrate orders"]["Power"] = list(\
                np.array(mc_grating_json[r_value]["Substrate orders"]["Power (s)"]) + 
                np.array(mc_grating_json[r_value]["Substrate orders"]["Power (p)"]))
    except:
        pass #mc_grating_json[r_value]["Substrate orders"]["Power"] = []
    
    return mc_grating_json
        

def add_total_power_1D(mc_grating_json):
    data = {k:v for k,v in mc_grating_json.items() if k.isnumeric()}
    
    for r_value in data.keys():
        total_r, total_t = total_power(mc_grating_json[r_value])
        absorption = 1-mc_grating_json[r_value]["Balance"]
        # Add to JSON File
        mc_grating_json[r_value]["Cover orders"]["Total Power"] = total_r
        try: # if exists (substrate transparent)
            mc_grating_json[r_value]["Substrate orders"]["Total Power"] = total_t
        except:
            pass #mc_grating_json[r_value]["Substrate orders"]["Total Power"] = 0
            
        mc_grating_json[r_value]["Absorption"] = absorption
        mc_grating_json[r_value]["Cover orders"]["Power"] = list(\
            np.array(mc_grating_json[r_value]["Cover orders"]["Power (s)"]) + 
            np.array(mc_grating_json[r_value]["Cover orders"]["Power (p)"]))
            
        try:
            mc_grating_json[r_value]["Substrate orders"]["Power"] = list(\
                np.array(mc_grating_json[r_value]["Substrate orders"]["Power (s)"]) + 
                np.array(mc_grating_json[r_value]["Substrate orders"]["Power (p)"]))
        except:
            pass #mc_grating_json[r_value]["Substrate orders"]["Power"] = []
            
    return mc_grating_json

def add_total_power_2D(mc_grating_json):
    data = {k:v for k,v in mc_grating_json.items() if k.lstrip('-').isnumeric()}
    
    for r_value in data.keys():
        for c_value in mc_grating_json[r_value].keys():
            # Get Total Power (Reflection, Transmission, Absoprtion)
            total_r, total_t = total_power(mc_grating_json[r_value][c_value])
            absorption = 1-mc_grating_json[r_value][c_value]["Balance"]
            # Add to JSON File
            mc_grating_json[r_value][c_value]["Cover orders"]["Total Power"] = total_r
            
            try: # if exists (substrate transparent)
                mc_grating_json[r_value][c_value]["Substrate orders"]["Total Power"] = total_t
            except:
                pass #mc_grating_json[r_value][c_value]["Substrate orders"]["Total Power"] = 0
            
            mc_grating_json[r_value][c_value]["Absorption"] = absorption
            
            mc_grating_json[r_value][c_value]["Cover orders"]["Power"] = list(\
                np.array(mc_grating_json[r_value][c_value]["Cover orders"]["Power (s)"]) + 
                np.array(mc_grating_json[r_value][c_value]["Cover orders"]["Power (p)"]))
            
            try: # if exists (substrate transparent)
                mc_grating_json[r_value][c_value]["Substrate orders"]["Power"] = list(\
                    np.array(mc_grating_json[r_value][c_value]["Substrate orders"]["Power (s)"]) + 
                    np.array(mc_grating_json[r_value][c_value]["Substrate orders"]["Power (p)"]))
            except:
                pass #mc_grating_json[r_value][c_value]["Substrate orders"]["Power"] = []
            
        
    return mc_grating_json


def total_power(fixed_wavelength_data):
    """
    Extract the total cover reflection or substrate transmission
    """
    try:
        sum_s_cov = np.sum(fixed_wavelength_data["Cover orders"]["Power (s)"])
        sum_p_cov = np.sum(fixed_wavelength_data["Cover orders"]["Power (p)"])
        total_cov = sum_s_cov+sum_p_cov
    except:
        #print("no cover orders")
        total_cov = 0
        
    try:
        sum_s_sub = np.sum(fixed_wavelength_data["Substrate orders"]["Power (s)"])
        sum_p_sub = np.sum(fixed_wavelength_data["Substrate orders"]["Power (p)"])
        total_sub = sum_s_sub+sum_p_sub
    except:
        #print("no substrate orders")
        total_sub = 0
        
    return [total_cov, total_sub]
