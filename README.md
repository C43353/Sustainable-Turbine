# Sustainable-Turbine

Files in order of use:
1. Profile_Plotting.py
  - Used to generate plots of the aerofoil shape profiles, and save them to the "Shape Profiles" file.

2. Interpolation_Tests.py
  - Just some initial tests on the type of interpolation to use, compares different methods against eachother.

3. CLD_Profile_Plotting.py
  - Used to generate plots of the aerofoil CLD profiles, and save them to the "CLD Profiles" file.
  - Also used to generate plots fo the aerofil CLD profiles compared to their interpolated version, and save them to the "Compare CLD Profiles" file.

4. CLD+Profile_Plotting.py
  - Used to generate plots of the aerofoil shape and CLD profiles side by side on one plot and save them to the "CLD + Shape Profiles" file.

5. Functions.py
  - This contains the functions used by the other calculation codes.
    - cld_func() creates functions that allow the interpolation of Cl and Cd values from a given CLD file from the Aerofoil-data.zip file.
    - nodal() takes the inputs for the wind turbine at a radial position and outputs the nodal values including nodal forces.
    - forces() integrates between the segments to find the torque and normal force for each segment of the blade.

6. Part_1.py
  - This is the proof of concept of the lecture notes, uses 1D analysis to find a suitable radius of the blade, then calculates the nodal outputs for a specific radial position on the blade.

7. BEM_500kW.py
  - Uses the lecture notes and variables stated in them to calculate the power output and other plots of a 500kW turbine at varying wind speeds and global pitch angles. (Should be compared to lecture slides to ensure correct calculations)

8. BEM_8MW.py
  - Uses the code from BEM_500kW.py, just changes the outer radius to 85m to allow larger power generation, the aerofoil profile, chord length and pitch angle will also vary in this compared to BEM_500kW.py.
