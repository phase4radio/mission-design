# Satellite Solar Power Budget
Calculator for estimating PV panel output power based on earth orbit, PV panel parameters and satellite attitude  

### Input parameters
* Orbit TLE
* Initial pose
* Angular velocities
* PV panel area per panel per side
* Number of PV panels per side
* PV panel efficiency

### Outputs
* Mean power coefficient for each side
* Total Mean power coefficient
* Mean power for each side
* Total Mean coefficient
* Power plot for each side

### Limitations
* Orbit is simulated via SGP4 using TLEs  
* Satellite shape is assumed to be a Rectangular cuboid  
* Individual PV panels of each side are assumed to have the same area  
* Spacecraft parameters are common for all orbits  

## Install

```bash
git clone https://gitlab.com/drid/satellite-solar-power-budget.git
pip3 install -r requirements.txt
```

Copy **parameters-sample.yaml** to **parameters.yaml** and customize according to needs.  
Multiple TLEs can be entered simulating different orbits for the same satellite
```bash
cp parameters-sample.yaml parameters.yaml
vim parameters.yaml  # or nano or whatever editor
python3 power_budget.py
```

## FAQ

#### PVs for a single side are not of equal size, how do I enter their sizes?

Calculate the total area of all the PVs of that side and use a PV count of 1

#### Is satellite rotation taken into account?
Yes, enter a non zero value for angular velocities  
Angular velocity of 0 means that the satellite does not rotate in respact to the velocity vector

#### How are axis oriented?
X+ Velocity vector  
Z+ Nadir

#### What units are used?
Rotation: rad/min  
Pose: Degrees  
Power: mW  
Efficiency: 0.00 - 1.00

## Attribution
Initial code by 
Yannis Koveos (https://gitlab.com/ykoveos)
Agis Zisimatos (https://gitlab.com/zisi)
Manthos Papamatthaiou (https://gitlab.com/papamat)

