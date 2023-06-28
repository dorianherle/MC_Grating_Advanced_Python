# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 14:46:17 2021

@author: Dorian
"""
import matplotlib.pyplot as plt 
import json
import numpy as np
import pandas as pd
from json_numpy_encoder_decoder import NumpyArrayDecoder
import trimesh
from geometry_cross_section_visualization import create_cross_section
from descartes import PolygonPatch
import cmath
import shapely
import copy
import re
import csv
# %% Visualization of the poyniting vector in 3D
from mayavi import mlab
from mayavi.mlab import *

# %%
## FAR- FIELD

def total_power(fixed_wavelength_data, direction="cover"):
        
        if direction == "cover":
            try:
                total_cov = fixed_wavelength_data["Cover orders"]["Total Power"]
            except:
                print("no cover orders")
                total_cov = 0
            
            return total_cov
        
        if direction == "substrate":
            try:
                total_sub = fixed_wavelength_data["Substrate orders"]["Total Power"]
            except:
                print("no substrate orders")
                total_sub = 0
            
            return total_sub
        
        if direction == "abs":
            try:
                abs_ = fixed_wavelength_data["Absorption"]
            except:
                print("no absorption data found -> please check simulation input.")
                abs_ = 0
            
            return abs_
        
 
def orders(data, direction="cover", total=False, plot=True):
   
    scan_param = data["ScanPar"]["Name"]

    wavelengths = [k for k in data.keys() if k.isnumeric()]
    
    orders = data[wavelengths[0]][direction.capitalize()+" orders"]["Orders"]
    
    orders_dict = {order:{scan_param:[],"power":[], "angle_n": [], "angle_p":[]} for order in orders}
    
    if not total:
        for w in wavelengths:
            orders_w = data[w][direction.capitalize()+" orders"]["Orders"]
            power_s_w = list(data[w][direction.capitalize()+" orders"]["Power (s)"])
            power_p_w = list(data[w][direction.capitalize()+" orders"]["Power (p)"])
            angle_N = list(data[w][direction.capitalize()+" orders"]["AngleN, deg"])
            angle_P = list(data[w][direction.capitalize()+" orders"]["AngleP, deg"])
            
            for order in orders_w:
                # iterate through each order
                order_idx = list(orders_w).index(order)
                power_s = float(power_s_w[order_idx])
             
                power_p = float(power_p_w[order_idx])
                power = power_s + power_p
                
                
                orders_dict[order]["power"].extend([power])
                orders_dict[order][scan_param].extend([float(w)])
                
                orders_dict[order]["angle_n"].extend([float(angle_N[order_idx])])
                orders_dict[order]["angle_p"].extend([float(angle_P[order_idx])])                
         
    if total:
        orders_dict = {}
        for w in wavelengths:
            orders_dict[w] = total_power(data[w], direction=direction)

            
    # Plot
    if plot and not total:
        plt.subplots(dpi=300)
        wavelengths_n = [float(w) for w in wavelengths]
        for order, w_and_power in orders_dict.items():
            w = w_and_power[scan_param]
            p = w_and_power["power"]
            plt.plot(w, p, label = order)

      
        plt.title(direction.capitalize()+" orders")
        plt.xlabel(scan_param)
        if direction == "cover":
            plt.ylabel("R")
        if direction == "substrate":
            plt.ylabel("T")
        plt.xlim(np.min(wavelengths_n), np.max(wavelengths_n))
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show()
    
    if plot and total:
        plt.subplots(dpi=300)
        wavelengths_n = [float(w) for w in wavelengths]
        total_power_all = [v for k,v in orders_dict.items()]
        plt.plot(wavelengths_n, total_power_all)
          
        plt.title(direction.capitalize()+"- Sum of All Orders")
        plt.xlabel(scan_param)
        plt.xlim(np.min(wavelengths_n), np.max(wavelengths_n))
        
        if direction == "cover":
            plt.ylabel("R")
        if direction == "substrate":
            plt.ylabel("T")
            
        plt.show()
    
    return orders_dict
    
def rta(data, plot=True, save="", exclude = None):
    scan_param = data["ScanPar"]["Name"]
    wavelengths = [k for k in data.keys() if k.isnumeric()]


    orders_dict = {"R":[], "T":[], "A":[], scan_param: wavelengths}
    for w in wavelengths:
        orders_dict["R"].append(float(total_power(data[w], direction="cover")))
        orders_dict["T"].append(float(total_power(data[w], direction="substrate")))
        orders_dict["A"].append(float(total_power(data[w], direction="abs")))
        
    orders_dict[scan_param] = [float(k) for k in data.keys() if k.isnumeric()]
    
    if plot:
        plt.figure(dpi=300)
        wavelengths_n = [float(w) for w in wavelengths]
    
        if exclude != "R":  
            plt.plot(wavelengths_n, orders_dict["R"], label="R")
        if exclude != "T":
            plt.plot(wavelengths_n, orders_dict["T"], label="T")
        if exclude != "A":  
            plt.plot(wavelengths_n, orders_dict["A"], label="A")
        
        if exclude != None:
            plt.title("R-T-A".replace(exclude, ""))
        else:
            plt.title("R-T-A")
        plt.xlabel(scan_param)
        plt.xlim(np.min(wavelengths_n), np.max(wavelengths_n))
        
        plt.ylabel("R/T/A")
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
        if save != "":
            plt.savefig(save, bbox_to_anchor=True)
            
        plt.ylim(0,1)
        plt.show()
        
        
    
    return orders_dict
    

def far_field_scattering(data, direction="cover", plot=True, title=""):
    scan_param = data["ScanPar"]["Name"]
    wavelengths = [k for k in data.keys() if k.isnumeric()]
    
    for w in wavelengths:
        
        df_single_w = data[w]
        result = df_single_w[direction.capitalize()+" orders"]
        
        ## Get Data
        angle_n = np.array(result['AngleN, deg'])
        angle_p = np.array(result['AngleP, deg'])
        power = np.array(result['Power'])
        
        ##  Far Field Scattering Plot
        angle_n_list = []
        angle_p_list = []
        
        for i,n in enumerate(angle_n):
            if n < 0:  
               #The choice of AngleP and AngleN 
               # is ambiguous for a specific order: 
               # the beam with AngleP2 = AngleP - 180 
               # and AngleN2 = - AngleN is directed exactly the same way. 
               angle_n_list.append(-n)
               angle_p_list.append(angle_p[i]-180)
            else:
                angle_n_list.append(n)
                angle_p_list.append(angle_p[i])
        
        r = np.array(angle_n_list)
        theta = np.radians(np.array(angle_p_list))
        #colors = np.divide(power,10**(-4))
        # np.nan is a float so you need to convert array to float before doing the boolean masking.
        #power = power.astype(float)
        #power[power<0.0001] = np.NaN
        
      
        
        palette = plt.get_cmap('jet')
        #palette.set_bad("white")
        
        fig = plt.figure(figsize=(8, 6), dpi=300)
        ax = fig.add_subplot(111, projection='polar')
        ax.set_rorigin(0)
        #c = 
        
        for i,a in enumerate(zip(-theta, r)):
            t, r = a
            c = ax.scatter(-t, r, c=power[i], s=100, cmap=palette, vmin=min(power), vmax=max(power), label = "{:.1f}, {:.1f}°".format(int(t*180/np.pi), r))
            #ax.annotate("{:.1f}, {:.1f}°".format(int(t*180/np.pi), r), xy=[t, r-1], fontsize=10)
            
        fig.colorbar(c)
        #plt.legend()
        # move title upwards, then adjust top spacing
        ax.set_title(direction.capitalize()+ " " + w, va='bottom', y=1.1, fontsize=15)
        plt.subplots_adjust(top=0.8)
        plt.show()
        plt.pause(0.001) 
        
# %% NEAR-FIELD 
     

path_material = r"C:\Users\BIEST\AppData\Roaming\MC Grating Software\Materials.dat"

# IMPORT MATERIAL 
def index_to_dielectricConst(n,k):
    # Source: https://www.researchgate.net/post/Kindly_elaborate_the_relation_between_Refractive_Index_and_Dielectric_Constant 
    return n**2-k**2, 2*n*k

def get_first_word_inside_brackets(s):
    match = re.search(r'\((.*?)\)', s)
    if match:
        return match.group(1).split()[0]
    else:
        return None

def schott_model(wavelength_um, c):
    refractive_index = np.sqrt(
        c[0] +
        c[1]*wavelength_um**2 +
        c[2]/wavelength_um**2 +
        c[3]/wavelength_um**4 +
        c[4]/wavelength_um**6 +
        c[5]/wavelength_um**8
    )
    return refractive_index

def sellmeier_model(wavelength_um, c):
    refractive_index = np.sqrt(
        1 +
        (c[0]*wavelength_um**2) / (wavelength_um**2 - c[3]) +
        (c[1]*wavelength_um**2) / (wavelength_um**2 - c[4]) +
        (c[2]*wavelength_um**2) / (wavelength_um**2 - c[5])
    )
    return refractive_index


def drude_model(wavelength, c):
    refractive_index = np.sqrt(1+(-c[1]+1j*wavelength*c[0]*c[1])/(c[0]**2+wavelength**-2))
    return refractive_index.real, refractive_index.imag


def herzberger_model(wavelength, c):
    L = 1 / (wavelength**2 - 0.028)
    return c[0] + c[1]*wavelength**2 + c[2]*wavelength**4 + c[3]*L + c[4]*L**2 + c[5]*L**3


datContent = []
with open(path_material, newline = '') as games:                                                                                          
    reader = csv.reader(games, delimiter='\t')
    for game in reader:
        datContent.append(game)


material_dict = {}
 # in um
wavelength_range = np.linspace(200, 3000, 301)/1000

for mat in datContent:
    name = mat[0]

    mat = list(filter(None, mat))  # remove all empty strings from the list
    model_type = get_first_word_inside_brackets(name)  # Parse the model type from the name
    coefficients = np.array([float(i) for i in mat[2:-2]])  # Parse coefficients 
    
    if model_type == 'Schott':
        n = schott_model(wavelength_range, coefficients)  # Calculate refractive index using Schott model
        k = np.zeros_like(n)  # assuming k = 0 for Schott model
    if model_type == 'Sellmeier':
        n = sellmeier_model(wavelength_range, coefficients)
        k = np.zeros_like(n)  # assuming k = 0 for Sellmeier model
    if model_type == 'Herzberger ':
        n = herzberger_model(wavelength_range, coefficients)
        k = np.zeros_like(n)  # assuming k = 0 for Herzberger model
    if model_type == 'Drude':
        coefficients = np.array([float(i) for i in mat[1:3]])  # Parse coefficients for Drude model
        n, k = drude_model(wavelength_range, coefficients)  # Calculate refractive index using Drude model

    if model_type == 'Table':
        w = np.array([float(i) for i in mat[3::3]])
        n = np.array([float(i) for i in mat[4::3]])
        k = np.array([float(i) for i in mat[5::3]])
        epsR, epsI = index_to_dielectricConst(n,k)
        material_dict[name] = {"w": w, "n": n, "k": k, "epsR": epsR, "epsI": epsI}
    else:
        epsR, epsI = index_to_dielectricConst(n,k)
        material_dict[name] = {"w": wavelength_range, "n": n, "k": k, "epsR": epsR, "epsI": epsI}
        
        
def is_number(string):
      try:
          float(string)
          return True
      except ValueError:
          return False
    
    
def extract_field(mc_grating_json, component="Sx"):
    # Extract field component from JSON output data
    # Returns field and associated points
    
    ## X 
    sc_name = mc_grating_json["ScanParC"]["Name"]
    

    ## Y
    s_name = mc_grating_json["ScanPar"]["Name"]
    
    # get index of component
    component_index = mc_grating_json["Components"].index(component)
    

    def cheat(df_list, s,sc, value):
        
        
        if float(sc) > 0:
            df_list.append([float(s),float(sc)-137/2,float(value)])
        if float(sc) == 0:
            df_list.append([float(s),float(sc)+137/2,float(value)])
            df_list.append([float(s),float(sc)-137/2,float(value)])
        if float(sc) < 0 and float(sc) != -137/2:
            df_list.append([float(s),float(sc)+137/2,float(value)])

            
       
        
        
        
    field = []
    df_list = []
    
    s_positions = [k for k in mc_grating_json.keys() if is_number(k)]
    
    for s in s_positions:
        sub_data = mc_grating_json[str(s)]
        temp = []
        for sc in sub_data.keys():
            value = sub_data[sc][component_index]
            temp.extend([value])
            #cheat(df_list, s,sc, value)
            df_list.append([float(s),float(sc),float(value)])
        field.append(temp)
        

    df = pd.DataFrame(df_list,
               columns =[s_name, sc_name, component])
    
    # Sort 
    df.sort_values([s_name, sc_name], inplace=True, ignore_index=True)
        
    return np.array(field), df


def complexe_modulo(z):
     a = z.real
     b = z.imag
     return np.sqrt(a**2+b**2)
 

def geometry_cross_section(mc_grating_json, scene, material=None, ax=None):
    
    if ax is None:
            ax = plt.gca()
            
    ## GET NORMAL
    normal = mc_grating_json['Normal_of_CrossSection']
    
    ## GET POSITION OF CROSS_SECTION
    pos = mc_grating_json['Pos_of_CrossSection']
         
    
    cross_section, bounds = create_cross_section(scene, material, pos=pos, normal=normal, plot=False)
    
    # If cross-section geometry is provided
    if scene != "": 
        # Add cross-section
        for item in cross_section.items(): 
            name, geom = item
            for i, polygom in enumerate(geom["polygon"]):
                
                
                patch = PolygonPatch(polygom.buffer(0),\
                                     fc="none", ec='k', alpha = 1)
                # else:
                #     patch = PolygonPatch(polygom.buffer(0),\
                #                          fc=material[name], ec='k', alpha = 1)
                ax.add_patch(patch)
                
    return ax, cross_section


def field(mc_grating_json, scene, material, component="Sx", plot=True, title="",\
               min_field="", max_field="", epsR =1, periodicity=1):
    

    
    # periodicity is used to show copies of the plot n times to visualize a periodic structure
    
   
    
    ## GET WAVELENGTH
    single_w = mc_grating_json["WaveL"]
    
    ## GET NORMAL
    normal = mc_grating_json['Normal_of_CrossSection']
    
    ## GET POSITION OF CROSS_SECTION
    pos = mc_grating_json['Pos_of_CrossSection']
         
    
     # IF GEOMETRY INCLUDED
    if scene != "" and material != "":
        cross_section, bounds = create_cross_section(scene, material, pos=pos, normal=normal, plot=False)
        
    ## X
    xlabel = mc_grating_json["ScanParC"]["Name"]
    initial_value_x = float(mc_grating_json["ScanParC"]["IniValue"])
    final_value_x = float(mc_grating_json["ScanParC"]["FinValue"])
    

    ## Y
    ylabel = mc_grating_json["ScanPar"]["Name"]
    initial_value_y = float(mc_grating_json["ScanPar"]["IniValue"])
    final_value_y = float(mc_grating_json["ScanPar"]["FinValue"])
    
    # IF ABSOLUTE VALUE 
    if "mod" in component:
        # get field type
        field_type = component.split(" mod")[0] 
        Re,re_df_field = extract_field(mc_grating_json, component=field_type+" Re")
        Im,_ = extract_field(mc_grating_json, component=field_type+" Im")
        field = complexe_modulo(Re + 1j *Im)
        
        # Rename component
        component = "|"+field_type+"|"

        # Create dataset
        s,sc,_ = re_df_field.to_numpy().T
        column_names = re_df_field.columns.to_list()[:-1] # avoid the field name
        column_names.extend([component]) # add modulus name

        df_field = pd.DataFrame(zip(s,sc,field.ravel()),columns=(column_names))
        
    elif "abs(" in component:
        type_ = component.split("abs(")[1]
        type_ = type_.split(")")[0]
        field_ex_im, df_field = extract_field(mc_grating_json, component= type_+"x Im")
        field_ex_re, _ = extract_field(mc_grating_json, component= type_+"x Re")
        
        field_ey_re, _ = extract_field(mc_grating_json, component=type_+"y Im")
        field_ey_im, _ = extract_field(mc_grating_json, component=type_+"y Re")
        
        field_ez_im, _ = extract_field(mc_grating_json, component=type_+"z Im")
        field_ez_re, _ = extract_field(mc_grating_json, component=type_+"z Re")
        
        field_e_complex_vector = np.array([-field_ex_re- field_ex_im*1j,
                                           field_ey_re+ field_ey_im*1j,
                                           field_ez_re+ field_ez_im*1j])
        
        field = np.sqrt(sum([complexe_modulo(component)**2 for component in  field_e_complex_vector]))
        
        s,sc,_ = df_field.to_numpy().T
        column_names = df_field.columns.to_list()[:-1] # avoid the field name
        column_names.extend([component]) # add modulus name

        df_field = pd.DataFrame(zip(s,sc,field.ravel()),columns=(column_names)) 
        
    elif component[0] == "|":
        # extract comp
        comp = component.split("|")[1]
        
        field_e_im, df_field = extract_field(mc_grating_json, component=comp+" Im")
        field_e_re, _ = extract_field(mc_grating_json, component=comp+" Re")
        
        field_e_complex_vector = np.array([field_e_re+ field_e_im*1j])
        
        field = np.sqrt(sum([complexe_modulo(component)**2 for component in  field_e_complex_vector]))
        
        s,sc,_ = df_field.to_numpy().T
        column_names = df_field.columns.to_list()[:-1] # avoid the field name
        column_names.extend([component]) # add modulus name

        df_field = pd.DataFrame(zip(s,sc,field.ravel()),columns=(column_names))
        
        
    elif "Im(eps) |E|^2" in component:
        
        component= "abs(E)"
        from shapely.geometry import Point
        
        field_ex_im, df_field = extract_field(mc_grating_json, component="Ex Im")
        field_ex_re, _ = extract_field(mc_grating_json, component="Ex Re")
        
        field_ey_re, _ = extract_field(mc_grating_json, component="Ey Im")
        field_ey_im, _ = extract_field(mc_grating_json, component="Ey Re")
        
        field_ez_im, _ = extract_field(mc_grating_json, component="Ez Im")
        field_ez_re, _ = extract_field(mc_grating_json, component="Ez Re")
        
        field_e_complex_vector = np.array([field_ex_re+ field_ex_im*1j,
                                           field_ey_re+ field_ey_im*1j,
                                           field_ez_re+ field_ez_im*1j])
        
        field = np.sqrt(sum([complexe_modulo(component)**2 for component in  field_e_complex_vector]))
        
        s,sc,_ = df_field.to_numpy().T
        
        column_names = df_field.columns.to_list()[:-1] # avoid the field name
        column_names.extend([component]) # add modulus name

        df = pd.DataFrame(zip(s,sc,field.ravel()),columns=(column_names))
        
        points = list(zip(s.tolist(), sc.tolist()))
        
        temp = np.array([0]*len(df[component]))
        for item in cross_section.items(): 
            name, geom = item
            material = geom["material"]
            
            w = material_dict[material]["w"]
            ns = material_dict[material]["n"]
            ks = material_dict[material]["k"]
            
           
            # plt.title(material)
            # plt.plot(w,ns, label="n")
            # plt.plot(w,ks, label="k")
            # plt.legend()
            # plt.pause(0.1)
            

            for i, polygon in enumerate(geom["polygon"]):
                to_remove = []
                for i, point in enumerate(points):
               
                    si,sci = point
                    if polygon.contains(Point(sci, si)):
                        
                        
                        is_x = df.iloc[:, 0] == si
                        is_z = df.iloc[:, 1] == sci
                        idx = np.where((is_z & is_x).values == True)

                        n = np.interp(float(single_w), w, ns)
                        k = np.interp(float(single_w), w, ks)
                        
                        # CORE
                        epsRe, epsIm = index_to_dielectricConst(n,k)
                        # eps_r =complex(epsR, epsIm)
                        # print(eps_r)
                        # c = 29979245
                        # epsilon_0 = 8.8541878128*10**-12
                        # omega = 2*np.pi*c/float(single_w)
                        # print(omega)
                        # sigma = -1j*omega*epsilon_0*(eps_r-1)
                        # print(sigma)
                        # TO SMALL
                        temp[idx] = epsIm # EPSILON IMAG
                        
                        #print(idx)
                        #print(temp[idx])
                        #df['epsR'][is_z & is_x] = epsR #material_dict[material]["epsR"]
                        
                        
                        to_remove.append(i)
                        
                points_c = []
                # Which points to keep?
                for i, p in enumerate(points):
                    if i not in to_remove:
                        points_c.append(p)
                # Reasign theses points to be checked 
                points = points_c
                
   
        df['epsImag'] = temp  
        df["Im(eps) |E|^2"] = df[component]*df['epsImag']
        

        
        field = df["Im(eps) |E|^2"].to_numpy().reshape(len(np.unique(s)),-1)
        
        df_field = df
        component= "Im(eps)* |E|^2"
        
        
        
    
    else:
        field, df_field = extract_field(mc_grating_json, component=component)
        # if component == "Ey Re" or component == "Ex Re" or component == "Hz Re": #fix due to coordinate missmatech between COMSOL and MC 
        #     field = -field
        
    
    from shapely.geometry import LineString
    # IF PLOT
    if plot: 
        
        # number of copies you want
        n_copies = periodicity
        
        # repeat the field
        repeated_field = np.tile(field, (1, n_copies))
        
        fig = plt.figure(dpi=300)
        ax = plt.gca()
        
        if min_field == "" or max_field == "":
            min_field = np.min(repeated_field)
            max_field = np.max(repeated_field)
     
        
        init_x = -(final_value_x - initial_value_x) * n_copies/2
        final_x = (final_value_x - initial_value_x) * n_copies/2
        # create field plot
        Q = ax.imshow(repeated_field,extent=[init_x,final_x, final_value_y,initial_value_y],\
                      interpolation="none", cmap='jet',\
                      aspect="equal", vmin=min_field,\
                      vmax = max_field)


        
        # If cross-section geometry is provided
        if scene != "": 
            aspect_ratio = 1
            # Add cross-section
            for item in cross_section.items(): 
                name, geom = item
                for i, polygon in enumerate(geom["polygon"]): 
                    # Scale the polygon
                    polygon_scaled = shapely.affinity.scale(polygon, xfact=aspect_ratio, origin='center')
                    for copy_index in range(n_copies):
                        polygon_copy = copy.deepcopy(polygon_scaled)
                        # translate the polygon copy to the center of each repeated field
                        polygon_copy = shapely.affinity.translate(polygon_copy, xoff=init_x+(final_value_x - initial_value_x)/2+(final_value_x - initial_value_x) * copy_index)
                        patch = PolygonPatch(polygon_copy, fc="none", ec="w", linestyle = "--", linewidth = 1, alpha = 1)
                        ax.add_patch(patch)


                    
        if title == "":
            plt.title(component+"\n@ "+single_w +" nm" + "; " + normal.upper() + " = " + str(pos))
        else:
            plt.title(title)
            
        ax.set_ylabel(ylabel)  
        ax.set_xlabel(xlabel)
        ax.set_xlim(init_x,final_x)
        ax.set_ylim(initial_value_y, final_value_y)
        ax.invert_yaxis()
    
        fig.colorbar(Q, ax=fig.get_axes())
        plt.show()
        plt.pause(0.001)
    
    return df_field, field, single_w




def extract_2D_data(mc_grating_json, data="Total Power", cover = True):
    
    r_name = mc_grating_json["ScanPar"]["Name"]
    c_name = mc_grating_json["ScanParC"]["Name"]
    
    data_ = {k:v for k,v in mc_grating_json.items() if k.lstrip('-').isnumeric()}
    
    dict_ = {}
    for r_value in data_.keys():
        
        temp_ = []
        for c_value in mc_grating_json[r_value].keys():
            if cover: 
                # get data
                d = mc_grating_json[r_value][c_value]['Cover orders'][data]
                
            else: 
                # get data
                d = mc_grating_json[r_value][c_value]['Substrate orders'][data]
      
            temp_.append(d)     
            
        dict_[r_value] = temp_
        
    df = pd.DataFrame(dict_)
    df.columns.name = r_name
   
    df.index = list(mc_grating_json[r_value].keys())
    df.index.name = c_name
    
    return df




def extract_line(mc_grating_json, component="Sx"):
    
    # get index of component
    component_index = mc_grating_json["Components"].index(component)
    
    # get values
    line = [v[component_index] for k,v in mc_grating_json.items() if is_number(k)]
    
    return np.array(line)

def field_line(mc_grating_json, component="Sx", plot=True):
    
    ## GET WAVELENGTH
    single_w = mc_grating_json["WaveL"]
    
    # get normal
    normal = mc_grating_json["Normal_of_CrossSection"]
    # get cross-section
    cross_section = "XYZ"
    cross_section = cross_section.replace(normal, "")
    # get position of line
    pos_line = mc_grating_json["Pos_of_Line"]
    
    # get boundaries
    xmin = mc_grating_json["ScanPar"]["IniValue"]
    xmax = mc_grating_json["ScanPar"]["FinValue"]
    n_points= mc_grating_json["ScanPar"]["NPoints"]
    points = np.linspace(xmin,xmax,n_points)
    
    if "abs(" in component:
        # get field type
        type_ = component.split("abs(")[1]
        field_type = type_.split(")")[0]
        
        if len(field_type) == 1:
            xRe = extract_line(mc_grating_json, component=field_type+"x Re")
            xIm = extract_line(mc_grating_json, component=field_type+"x Im")
            yRe = extract_line(mc_grating_json, component=field_type+"y Re")
            yIm = extract_line(mc_grating_json, component=field_type+"y Im")
            zRe = extract_line(mc_grating_json, component=field_type+"z Re")
            zIm = extract_line(mc_grating_json, component=field_type+"z Im")
            
            #line = complexe_modulo(ExRe + 1j *ExIm)
            
            field_complex_vector = np.array([-xRe- xIm*1j,
                                              yRe+ yIm*1j,
                                              zRe+ zIm*1j])
            
            line = np.sqrt(sum([complexe_modulo(component)**2 for component in  field_complex_vector]))
        
        else:
            ExRe = extract_line(mc_grating_json, component=field_type+" Re")
            ExIm = extract_line(mc_grating_json, component=field_type+" Im")
            line = complexe_modulo(ExRe + 1j *ExIm)
            
        # Rename component
        component = "|"+field_type+"|"
        
    else:
        line = extract_line(mc_grating_json, component=component)

        
    if plot:
        xlabel = mc_grating_json["ScanPar"]["Name"]
        line_pos_cood = cross_section.replace(xlabel[0], "")
        line_pos_cood = line_pos_cood+"/Per"+line_pos_cood
        
        
        plt.figure(dpi=300)
        plt.title(component+" @" + single_w +"nm\nCross-Section: "+cross_section+"; Line Position : " + pos_line + " [" +line_pos_cood+"]")
        plt.plot(points, line)
        plt.xlabel(xlabel)
        plt.xlim(xmin,xmax)
    
    return line, points, single_w, cross_section, pos_line
    


def poynting_vector_plot(field_sx,field_sz,points, skip_every=8,
                title = "", xlabel="", ylabel="", ax = None,  start_points = None):
    

     xx, zz = points
     
     if ax is None:
            ax = plt.gca()
     
     # if start_points is None:
     #     x_min = np.min(xx)
     #     x_max = np.max(xx)
     #     z_min = np.min(zz)
     #     z_max = np.max(zz)
         
     #     start_points = [[i,z_min-10] for i in np.linspace(x_min, x_max)]
            
            
     # Reduce arrow density in plot
     skip2D=(slice(None,None,skip_every),slice(None,None,skip_every))
 
     plt.figure(dpi=300)
     
     lenght = np.sqrt(field_sx**2+field_sz**2)
     ax.quiver(xx[skip2D], zz[skip2D],\
               field_sx[skip2D], -field_sz[skip2D] ) #cmap="jet", 
               # angles="xy", scale_units='xy', scale=1)
     # Note: inverting a data axis will correspondingly invert the arrows only with angles='xy'.
     # But if use, then scaling issue: https://stackoverflow.com/questions/69146016/python-quiver-angles-issue-what-is-the-difference-between-angles-uv-and-angle/69146753?noredirect=1#comment122225592_69146753
     # Therefore added negative sign manually
     
     
     ax.invert_yaxis()  
     ax.set_title(title)
     ax.set_xlabel(xlabel)
     ax.set_ylabel(ylabel)
    
     return ax


def poynting_vector(mc_grating_json, scene, material, plot=True, plot3D=True, cross_section=True,\
                    title="",\
                    min_field="", max_field=""):
    
    # get fields
    sx = field(mc_grating_json, scene, material, component="Sx", plot=False)
    sy = field(mc_grating_json, scene, material, component="Sy", plot=False)
    sz = field(mc_grating_json, scene, material, component="Sz", plot=False)
        
    # extract pos_param,row_param (x),column_param (y) names
    pos_param_name = mc_grating_json["Normal_of_CrossSection"]
    s_name = mc_grating_json["ScanPar"]["Name"] # Z
    sc_name = mc_grating_json["ScanParC"]["Name"] # Y
    
    # extract position value
    pos = float(mc_grating_json["Pos_of_CrossSection_Value"])
    
    # extract x,y,z,u,v,w
    s,sc,v0 = sx[0].to_numpy().T # Y,Z
    s,sc,v1 = sy[0].to_numpy().T
    s,sc,v2 = sz[0].to_numpy().T
    p = np.array([pos]*len(s))
    
    uvw_dict = {
        "X": v0,
        "Y": v1,
        "Z": v2
        }
    
    df_list = zip(p,s,sc,\
                  uvw_dict[pos_param_name],\
                  uvw_dict[s_name],\
                  uvw_dict[sc_name]) # X,Y,Z
    
    def vector_component_name(input_):
            if input_ == "X":
                return "U"
            if input_ == "Y":
                return "V"
            if input_ == "Z":
                return "W"
    
    df_3D = pd.DataFrame(df_list,
               columns =[pos_param_name,\
                         s_name,\
                         sc_name,
                         vector_component_name(pos_param_name),\
                         vector_component_name(s_name),\
                         vector_component_name(sc_name)])
    
    # Sort 
    df_3D.sort_values([pos_param_name, s_name, sc_name],\
                   inplace=True, ignore_index=True)
    
    
    if plot3D:
        # Extract coordinates and vector components
        x = df_3D["X"].to_numpy()
        y = df_3D["Y"].to_numpy()
        z = df_3D["Z"].to_numpy()
        u = df_3D["U"].to_numpy()
        v = df_3D["V"].to_numpy()
        w = df_3D["W"].to_numpy()
        
        # Plot
        mlab.figure('Original Poynting Vector Field', bgcolor=(0, 0, 0))
        mlab.quiver3d(x, y, z, u, v, w, scale_factor=10)
        mlab.view(azimuth=270, elevation=90, roll=180)
        mlab.show()
    

    if cross_section:
        # Extract vector components
        
        u = df_3D["U"].to_numpy()
        v = df_3D["V"].to_numpy()
        w = df_3D["W"].to_numpy()
        
        # get normal
        if pos_param_name == "X": n = np.array([1,0,0])
        if pos_param_name == "Y": n = np.array([0,1,0])
        if pos_param_name == "Z": n = np.array([0,0,1])
        
        # project vector onto cross-section plane
        
        def projection_onto_plane(n, u):
            # Source: https://www.geeksforgeeks.org/vector-projection-using-python/
            # n: normal vector of plane
            # u: 3d vector to be projected

            # finding norm of the vector n 
            n_norm = np.sqrt(sum(n**2))    
               
            # Apply the formula as mentioned above
            # for projecting a vector onto the orthogonal vector n
            # find dot product using np.dot()
            proj_of_u_on_n = (np.dot(u, n)/n_norm**2)*n
              
            # subtract proj_of_u_on_n from u: 
            # this is the projection of u on Plane P
            return u - proj_of_u_on_n
        
        # projected vectors onto plane -> still 3D vectors
        pr_plane = np.array([projection_onto_plane(n, np.array(vect)) for vect in list(zip(u,v,w))])

        # remove constant 3D component -> TODO: generalize for other projection planes
        column_index_to_remove = np.where(n==1)[0][0]
        pr_2D = np.delete(pr_plane, column_index_to_remove, 1)
        
        # create database 
        v1_2D,v2_2D = zip(*pr_2D)
        v1_2D = np.array(v1_2D)
        v2_2D = np.array(v2_2D)
        
        df_list_cross = zip(sc,s,v1_2D,v2_2D) 
        
        df_2D = pd.DataFrame(df_list_cross,
               columns =[s_name,
                         sc_name,\
                         vector_component_name(s_name),\
                         vector_component_name(sc_name)])
    
        # Sort 
        df_2D.sort_values([s_name, sc_name],\
                       inplace=True, ignore_index=True)
        
        if plot:
            
            # Plot
            plt.figure(dpi=300)
            ax = plt.gca()
            ax.invert_yaxis()
            
            
            
            max_arrows = 30
            skip = int(len(np.unique(sc))/max_arrows)
            lenght = np.sqrt(v1_2D**2+v2_2D**2)
            # Normalized by lenght and color is lenght
            # ax.quiver(sc[::skip],s[::skip], v1_2D[::skip]/lenght[::skip], v2_2D[::skip]/lenght[::skip], 
            #            lenght[::skip], cmap="jet",linewidths=1, angles='xy', units ="xy", scale =0.08)
            
            ax.quiver(sc[::skip],s[::skip], v1_2D[::skip]/max(lenght[::skip]), v2_2D[::skip]/max(lenght[::skip]), 
                       lenght[::skip], cmap="jet",linewidths=1, angles='xy', units ="xy", scale =0.04)
               
           
            
            # If cross-section geometry is provided
            
            # IF GEOMETRY INCLUDED
            if scene != "" and material != "":
                pos_cross = mc_grating_json["Pos_of_CrossSection"]
                cross_section, bounds = create_cross_section(scene, material,\
                                                             pos=pos_cross,\
                                                             normal=pos_param_name,\
                                                             plot=False)

            if scene != "": 
                # Add cross-section
                for item in cross_section.items(): 
                    name, geom = item
        
                    for polygom in geom['polygon']:
                        patch = PolygonPatch(polygom.buffer(0),\
                                              fc="none", ec='k', alpha = 1)
                        ax.add_patch(patch)
    

            
            if title == "": title = "Projected Poynting Vector @" + str(pos_cross) 
            plt.title(title)
            ax.set_xlabel(sc_name)
            ax.set_ylabel(s_name)
            
            
            plt.pause(0.01)
            plt.show()
            
        return df_3D, df_2D
    
    return df_3D


def poynting_vector_from_fields(mc_grating_json, scene, material, plot=True, plot3D=True, cross_section=True,\
                    title="",\
                    min_field="", max_field=""):
    
    # get fields
    sx = field(mc_grating_json, scene, material, component="Ex", plot=False)
    sy = field(mc_grating_json, scene, material, component="Sy", plot=False)
    sz = field(mc_grating_json, scene, material, component="Sz", plot=False)
        
    # extract pos_param,row_param (x),column_param (y) names
    pos_param_name = mc_grating_json["Normal_of_CrossSection"]
    s_name = mc_grating_json["ScanPar"]["Name"] # Z
    sc_name = mc_grating_json["ScanParC"]["Name"] # Y
    
    # extract position value
    pos = float(mc_grating_json["Pos_of_CrossSection_Value"])
    
    # extract x,y,z,u,v,w
    s,sc,v0 = sx[0].to_numpy().T # Y,Z
    s,sc,v1 = sy[0].to_numpy().T
    s,sc,v2 = sz[0].to_numpy().T
    p = np.array([pos]*len(s))
    
    uvw_dict = {
        "X": v0,
        "Y": v1,
        "Z": v2
        }
    
    df_list = zip(p,s,sc,\
                  uvw_dict[pos_param_name],\
                  uvw_dict[s_name],\
                  uvw_dict[sc_name]) # X,Y,Z
    
    def vector_component_name(input_):
            if input_ == "X":
                return "U"
            if input_ == "Y":
                return "V"
            if input_ == "Z":
                return "W"
    
    df_3D = pd.DataFrame(df_list,
               columns =[pos_param_name,\
                         s_name,\
                         sc_name,
                         vector_component_name(pos_param_name),\
                         vector_component_name(s_name),\
                         vector_component_name(sc_name)])
    
    # Sort 
    df_3D.sort_values([pos_param_name, s_name, sc_name],\
                   inplace=True, ignore_index=True)
    
    
    if plot3D:
        # Extract coordinates and vector components
        x = df_3D["X"].to_numpy()
        y = df_3D["Y"].to_numpy()
        z = df_3D["Z"].to_numpy()
        u = df_3D["U"].to_numpy()
        v = df_3D["V"].to_numpy()
        w = df_3D["W"].to_numpy()
        
        # Plot
        mlab.figure('Original Poynting Vector Field', bgcolor=(0, 0, 0))
        mlab.quiver3d(x, y, z, u, v, w, scale_factor=0.5)
        mlab.view(azimuth=270, elevation=90, roll=180)
        mlab.show()
    

    if cross_section:
        # Extract vector components
        
        u = df_3D["U"].to_numpy()
        v = df_3D["V"].to_numpy()
        w = df_3D["W"].to_numpy()
        
        # get normal
        if pos_param_name == "X": n = np.array([1,0,0])
        if pos_param_name == "Y": n = np.array([0,1,0])
        if pos_param_name == "Z": n = np.array([0,0,1])
        
        # project vector onto cross-section plane
        
        def projection_onto_plane(n, u):
            # Source: https://www.geeksforgeeks.org/vector-projection-using-python/
            # n: normal vector of plane
            # u: 3d vector to be projected

            # finding norm of the vector n 
            n_norm = np.sqrt(sum(n**2))    
               
            # Apply the formula as mentioned above
            # for projecting a vector onto the orthogonal vector n
            # find dot product using np.dot()
            proj_of_u_on_n = (np.dot(u, n)/n_norm**2)*n
              
            # subtract proj_of_u_on_n from u: 
            # this is the projection of u on Plane P
            return u - proj_of_u_on_n
        
        # projected vectors onto plane -> still 3D vectors
        pr_plane = np.array([projection_onto_plane(n, np.array(vect)) for vect in list(zip(u,v,w))])

        # remove constant 3D component -> TODO: generalize for other projection planes
        column_index_to_remove = np.where(n==1)[0][0]
        pr_2D = np.delete(pr_plane, column_index_to_remove, 1)
        
        # create database 
        v1_2D,v2_2D = zip(*pr_2D)
        v1_2D = np.array(v1_2D)
        v2_2D = np.array(v2_2D)
        
        df_list_cross = zip(sc,s,v1_2D,v2_2D) 
        
        df_2D = pd.DataFrame(df_list_cross,
               columns =[s_name,
                         sc_name,\
                         vector_component_name(s_name),\
                         vector_component_name(sc_name)])
    
        # Sort 
        df_2D.sort_values([s_name, sc_name],\
                       inplace=True, ignore_index=True)
        
        if plot:
            
            # Plot
            plt.figure(dpi=300)
            ax = plt.gca()
            ax.invert_yaxis()
            
            max_arrows = 10
            skip = int(len(np.unique(sc))/max_arrows)
            lenght = np.sqrt(v1_2D**2+v2_2D**2)
            # Normalized by lenght and color is lenght
            ax.quiver(sc[::skip],s[::skip], v1_2D[::skip]/lenght[::skip], v2_2D[::skip]/lenght[::skip], 
                       lenght[::skip], angles="xy", cmap="jet", headwidth=50, headlength=20)
               
           
            
            # If cross-section geometry is provided
            
            # IF GEOMETRY INCLUDED
            if scene != "" and material != "":
                pos_cross = mc_grating_json["Pos_of_CrossSection"]
                cross_section, bounds = create_cross_section(scene, material,\
                                                             pos=pos_cross,\
                                                             normal=pos_param_name,\
                                                             plot=False)

            if scene != "": 
                # Add cross-section
                for item in cross_section.items(): 
                    name, geom = item
                    for i, polygom in enumerate(geom): 
                        patch = PolygonPatch(polygom.buffer(0),\
                                              fc="none", ec='k', alpha = 1)
                        ax.add_patch(patch)


            
            if title == "": title = "Projected Poynting Vector @" + str(pos) 
            plt.title(title)
            ax.set_xlabel(sc_name)
            ax.set_ylabel(s_name)
            
            plt.pause(0.01)
            plt.show()
            
        return df_3D, df_2D
    
    return df_3D
            
        
        


        
        
        
    



#%%


            
def plot_2D(mc_grating_json, data="Total Power", cover = True, title="", vmin=0, vmax=0, plot=True):
    
    df = extract_2D_data(mc_grating_json, data=data, cover = cover)
    
    r_values = [float(elem) for elem in df.keys()]
    c_values = [float(elem) for elem in df.index]
    
    min_r = float(min(r_values))
    max_r = float(max(r_values))
    
    min_c = float(min(c_values))
    max_c = float(max(c_values))
    
    if vmin == vmax:
        vmin = np.min(df.values)
        vmax = np.max(df.values)
    
    if plot:
        plt.figure(dpi=300)
        plt.imshow(df, origin="lower",\
                    extent=[\
                        min_r,max_r,\
                        min_c, max_c],\
                    cmap="jet", aspect='auto', vmin=vmin, vmax=vmax)
        plt.colorbar()
        plt.xlabel(df.columns.name)
        plt.ylabel(df.index.name)
        
        if title!="":
            plt.title(title)
            
    return df
    




    
def vector_plot(field_sx,field_sz,points, skip_every=8,
                title = "", xlabel="", ylabel=""):
     
     xx, zz = points
     # Reduce arrow density in plot
     skip2D=(slice(None,None,skip_every),slice(None,None,skip_every))
 
     plt.figure(dpi=300)
     
     lenght = np.sqrt(field_sx**2+field_sz**2)
     Q = plt.quiver(xx[skip2D], zz[skip2D],\
                    field_sx[skip2D], -field_sz[skip2D], cmap="jet")
                    #angles="xy", scale_units='xy', scale=1)
     # Note: inverting a data axis will correspondingly invert the arrows only with angles='xy'.
     # But if use, then scaling issue: https://stackoverflow.com/questions/69146016/python-quiver-angles-issue-what-is-the-difference-between-angles-uv-and-angle/69146753?noredirect=1#comment122225592_69146753
     # Therefore added negative sign manually
     plt.gca().invert_yaxis()  
     plt.title(title)
     plt.xlabel(xlabel)
     plt.ylabel(ylabel)
     plt.pause(0.01)
     plt.show()

if __name__ == "__main__":
    
    
    def open_output(path):
        # Get MC Grating Result
        with open(path+'_output_scan_type_2'+'.json') as json_file:
            output = json.load(json_file)
        
        return output
    
    def open_output(path):
        # Get MC Grating Result
        with open(path+'.json') as json_file:
            output = json.load(json_file, cls=NumpyArrayDecoder)
        
        return output
    
    
    
    json_output = open_output("fields")
    # points = np.meshgrid(np.linspace(0, 1, 101), np.linspace(-200, 77, 101))
    
    # for position,field in json_output.items():
    #     # Since zx plane 
    #     points, field_sx = extract_field(field,  component="Sx", plot=False)
    #     _, field_sz = extract_field(field, component="Sz", plot=False)
    #     xx,zz = points
    #     vector_plot(field_sx,field_sz,points, title="Position: "+ position+ "\nVector Field Sx-Sz", xlabel="X/Per(x)", ylabel="Z")

    #
    
    field =  json_output['0.5555555555555556']
    points, field_sx = extract_field(field,  component="Sx", plot=False)
    _, field_sz = extract_field(field, component="Sz", plot=False)
    
    xx, zz = points
    xx = (xx-np.min(xx))/(np.max(xx)-np.min(xx))
    zz = (zz-np.min(zz))/(np.max(zz)-np.min(zz))
    # Plotting stream plot
    b = [0.]*20
    a = np.linspace(0,1,20)
    #b = [min(Z[skip1D])+20]*len(X[skip1D])+list(Z[skip1D][3:-2])+list(Z[skip1D][2:-2])
    skip_every = 1
    skip1D= slice(None,None,skip_every)
    seed_points = list(zip(a,b))
    plt.streamplot(xx, zz, field_sx, field_sz, 
                   cmap='jet', arrowstyle='-|>',
                   density=(len(xx[0][skip1D])/20,len(zz[0][skip1D])/20), start_points=seed_points)
    plt.gca().invert_yaxis()  
    # import plotly.graph_objs as go
    # import plotly.express as px
    # import pandas as pd
    # import numpy as np
    # import plotly.io as pio
    # import pickle
    # pio.renderers.default='browser'
    
    # # Axis
    # X = np.linspace(0,1,101)
    # Y = np.linspace(0,1,10)
    # Z = np.linspace(-500,200,101)
    
    # # Initialize List of all fiels
    # DX = []
    # DY = []
    # DZ = []
    
    # for position,field in json_output.items():
             
    #         _, dx = extract_field(field, component="Sx", plot=False)
    #         _, dy = extract_field(field, component="Sz", plot=False)
    #         _, dz = extract_field(field, component="Sz", plot=False)
            
    #         # Make them numpy imediatley
    #         dx = np.array(dx)
    #         dy = np.array(dy)
    #         dz = -np.array(dz)
            
    #         # Apppend
    #         DX.append(dx)
    #         DY.append(dy)
    #         DZ.append(dz)
    
    # #Convert to numpy
    
    # DX = np.array(DX)
    # DY = np.array(DY)
    # DZ = np.array(DZ)
        
        
    # # Create 3D Quiver Plot with color gradient
    # # Source: https://stackoverflow.com/questions/65254887/how-to-plot-with-matplotlib-a-3d-quiver-plot-with-color-gradient-for-length-giv
    # def plot_3d_quiver(x, y, z, u, v, w):
    #     # COMPUTE LENGTH OF VECTOR -> MAGNITUDE
    #     c = np.sqrt(np.abs(v) ** 2 + np.abs(u) ** 2 + np.abs(w) ** 2)
        
    #     c = (c.ravel() - c.min()) / c.ptp()
    #     # Repeat for each body line and two head lines
    #     c = np.concatenate((c, np.repeat(c, 2)))
    #     # Colormap
    #     c = plt.cm.jet(c)
    
    #     fig = plt.figure(dpi =300)
    #     ax = fig.gca(projection='3d')
    #     ax.quiver(x,y,z,u,v,w, colors=c, length=0.2, arrow_length_ratio=0.7,
    #               alpha=1)
    #     plt.gca().invert_zaxis()
    #     plt.show()
        
    
    
    # # Create Mesh !
    # xi, yi, zi = np.meshgrid(X, Y, Z, indexing='xy')
    # skip_every = 20
    # skip_slice = 10
    # skip3D=(slice(None,None,skip_slice),slice(None,None,skip_every),slice(None,None,skip_every))
    
    # # Source: https://stackoverflow.com/questions/68690442/python-plotting-3d-vector-field
    
    # plot_3d_quiver(xi[skip3D], yi[skip3D], zi[skip3D]/1000, DX[skip3D], DY[skip3D], DZ[skip3D])
    

    # # # Axis
    # # X = np.linspace(0,1,101)
    # Y = np.linspace(0,1,10)
    # Z = np.linspace(-0.2,0.077,101)
    
    # xpos,ypos = np.meshgrid(X[::5],Y, indexing="xy")
    # #xpos = xpos.reshape(1,-1)[0]
    # #ypos = ypos.reshape(1,-1)[0]
    # xpos = np.ravel(xpos)
    # ypos = np.ravel(ypos)
    # # Initialize List of all fields
    # DX = []
    # DY = []
    # DZ = []

    # for position,field in json_output.items():

    #     _, dx = extract_field(field, component="Sx", plot=False)
    #     _, dy = extract_field(field, component="Sz", plot=False)
    #     _, dz = extract_field(field, component="Sz", plot=False)
        
    #     # Make them numpy imediatley
    #     dx = np.array(dx)
    #     dy = np.array(dy)
    #     dz = -np.array(dz)

    #     # Apppend
    #     DX.append(dx)
    #     DY.append(dy)
    #     DZ.append(dz)

       
    # #Convert to numpy
    # move_i = [0, 1, 2]
    # move_e = [1, 2, 0]
    # DX = np.moveaxis(np.array(DX), move_i, move_e)
    # DY = np.moveaxis(np.array(DY), move_i, move_e)
    # DZ = np.moveaxis(np.array(DZ), move_i, move_e)
    
    # u1 = np.ravel(DX)
    # v1 = np.ravel(DY)
    # w1 = np.ravel(DZ)
    
    # np_df = np.array([u1, v1, w1])
    # vecf_norm = np.linalg.norm(np_df, 2, axis=0)
    # max_norm = np.max(vecf_norm)
    # min_norm = np.min(vecf_norm)
    
    # u2 = u1 * (vecf_norm - min_norm) / (max_norm - min_norm)
    # v2 = v1 * (vecf_norm - min_norm) / (max_norm - min_norm)
    # w2 = w1 * (vecf_norm - min_norm) / (max_norm - min_norm)
    
    # # Create Mesh !
    # xi, yi, zi = np.meshgrid(X, Y, Z, indexing="ij")


    # data_plot = [go.Streamtube(
    #     x = np.ravel(xi),
    #     y = np.ravel(yi),
    #     z = np.ravel(zi),
    #     u = u2,
    #     v = v2,
    #     w = w2,
    #     starts = dict(                           #Determines the streamtubes starting position.
    #         x=xpos,
    #         y=ypos,
    #         z=np.array([0.077]*len(xpos)
    #         )),
    #     #sizeref = 0.3,
    #     colorscale = 'jet',
    #     showscale = True,
    #     maxdisplayed = 300                      #Determines the maximum segments displayed in a streamtube.
    #     )]
    
    # fig = go.Figure(data=data_plot)
    # fig.update_scenes(zaxis_autorange="reversed")
    # fig.show()
        
    
    #df = extract_2D_data(json_output,cover = False)
    #far_field_scattering_plot(json_output["Fixed"], title="Far Field Scattering Plot")
    