# MC_Grating_Advanced_Python

This project contains advanced Python code for MC Grating. Follow the setup instructions below to get started.

## SETUP

### 1. Download & Install Anaconda

Anaconda is a popular Python distribution for data science and machine learning. Download and install it from the [official website](https://www.anaconda.com/).

### 2. Download & Unzip the MC Grating Python code

Download the MC Grating Python code from this GitHub repository and unzip it. 

![Download & Unzip Code](https://github.com/dorianherle/MC_Grating_Advanced_Python/blob/main/visualization_readme/Screenshot%202023-07-05%20113334.png)

### 3. Create a new Conda environment (Optional but Recommended) 

It's a good practice to create a separate environment for each project. To create a new environment named `mc_grating` with Python version 3.9, use the following command:

```
conda create -n mc_grating python=3.9
```

After running the command, Conda will ask for your permission to proceed with the installation. Type `y` and hit Enter to proceed. Once the environment is created, activate it using the following command:

```
conda activate mc_grating
```

### 4. Install Required Packages

First, navigate to the directory containing the extracted code:

```
cd "<path_to_folder>"
```

Replace <path_to_folder> with the actual path to the folder where you unzipped the MC Grating Python code, such as:

```
cd "C:\Users\MC_Grating\Downloads\MC_Grating_Advanced_Python-main\MC_Grating_Advanced_Python-main"
```

Then, install the necessary Python packages listed in the `requirements.txt` file:

```
pip install -r requirements.txt
```

### 5. Activate the Environment in Spyder

In order to use the mc_grating environment in Spyder, you need to activate it as shown in the images below:

   ![download](https://github.com/dorianherle/MC_Grating_Advanced_Python/blob/main/visualization_readme/Screenshot%202023-07-05%20114829.png
   )
   ![download](https://github.com/dorianherle/MC_Grating_Advanced_Python/blob/main/visualization_readme/envs.png)
   

### 6. Run the Test Script
Finally, open and run the test_all.py script to verify the setup.



# EXAMPLE:
   


```
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

from mc_grating_python_ver1.mc_grating_python import mc_grating_python as mc
```

## Create a geometry

```
def geometry(parameters):

    if parameters == "": # TO TEST
        parameters = [137,98,77]
        
    if len(parameters) != 3 and parameters != "":
        raise ValueError("List of input parameters needs to be 3")
              
    period, pillar_d, pillar_h\
        = parameters

    # geometry 
    geo = mc.geo(period, period)
    
    # Cover Layer
    geo.cover("Air (Special Formula)")
    
    #---PILLAR---#
    geo.cylinder(position=[0,0,0], diameter=pillar_d, height=pillar_h, 
                 name = "pillar", material = "Silicon (Table)", order = 1)
                 

    # Substrate
    geo.substrate("Fused Silica (Sellmeier)")
        
    
    # get geometry
    geometry = geo.print_g() 
    
    return geometry
```

## Visualize the Geometry 

```
color={
           'Amorphous_Si_CMI (Table GS)': (0.3,0.3,0.3),
           'Amorphous_Si_CMI_Lossless': (0.3,0.3,0.1),
           'Silicon (Table)' : (0.4,0.4,0.4),
           'Fused Silica (Sellmeier)':(0,0,1),
           "SiO2 (Table GS)": (0.7,0.1,0.7),
           "Al (Table GS)": (0.1,0.5,0.1),
           "aSi_k_half": (0,1,0),
           "Al2O3 (Table GS)": (0,0,0.5),
           "Air (Special Formula)": (1,1,1),
            }


geometry_dict = geometry("")
visual = mc.show(geometry_dict, color=color)
```
### 2D
```
cs = visual.show_cross_section(normal="y")
```

![2D](https://github.com/dorianherle/MC_Grating_Advanced_Python/blob/main/visualization_readme/2d.jpg)

### 3D

```
visual.show_3D()
```
![3D](https://github.com/dorianherle/MC_Grating_Advanced_Python/blob/main/visualization_readme/3d.jpg)

## Simulate 

```
sim = mc.sim(geometry_dict, single_wavelengh=400, nbr_oders_x=11, nbr_orders_y=11)
```
### Far-Field
```
ff = sim.far_field(parameter_row = "wavelength", \
                   start_row=400, end_row=700, nbr_of_points_row=200)
    
res_ff = mc.run(ff)
```
### Near-Field
```
nf = sim.near_field(start_row=-200, end_row=200, plane="zy",  number_of_points = [200,50],\
                    position_of_plane=0.5, output_format="amplitude_re_im_and_power_flow")
 
res_nf = mc.run(nf)
```

## Plot Results

### Far-Field
```
p_ff = mc.plot(res_ff)
rta = p_ff.rta() # reflection - transmission - absorption plot
```
![rta](https://github.com/dorianherle/MC_Grating_Advanced_Python/blob/main/visualization_readme/rta.png)

### Near-Field
```
p_nf = mc.plot(res_nf, geometry_dict) 
sz = p_nf.field(component="Sz", plot=True, title="")
abs_e = p_nf.field(component="abs(E)", plot=True, title="|E|")
ex_mod = p_nf.field(component="Ex mod", plot=True, title="")
im_e2 = p_nf.field(component="Im(eps) |E|^2", plot=True, title="", periodicity=1)
```
![sz](https://github.com/dorianherle/MC_Grating_Advanced_Python/blob/main/visualization_readme/sz.png)
![abs_e](https://github.com/dorianherle/MC_Grating_Advanced_Python/blob/main/visualization_readme/abs_e.png)
![abs_ex](https://github.com/dorianherle/MC_Grating_Advanced_Python/blob/main/visualization_readme/abs_ex.png)
![a](https://github.com/dorianherle/MC_Grating_Advanced_Python/blob/main/visualization_readme/absorption.png)
