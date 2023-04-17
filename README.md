# Sustainable-Turbine

Files in order of use:
1. Shape_Plotting.py
  - Used to generate plots of the aerofoil shape profiles, and save them to the "Shape Profiles" file.

2. Interpolation_Tests.py
  - Just some initial tests on the type of interpolation to use, compares different methods against each other.

3. CLD_Profile_Plotting.py
  - Used to generate plots of the aerofoil CLD profiles, and save them to the "CLD Profiles" file.
  - Also used to generate plots fo the aerofil CLD profiles compared to their interpolated version, and save them to the "Compare CLD Profiles" file.

4. CLD+Shape_Plotting.py
  - Used to generate plots of the aerofoil shape and CLD profiles side by side on one plot and save them to the "CLD + Shape Profiles" file.

5. Functions.py
  - This contains the functions used by the other calculation codes.
    - cld_func() creates functions that allow the interpolation of Cl and Cd values from a given CLD file from the Aerofoil-data.zip file.
    - nodal() takes the inputs for the wind turbine at a radial position and outputs the nodal values including nodal forces.
    - forces() integrates between the segments to find the torque and normal force for each segment of the blade.
    - nodal_twist() for a given wind speed and nodal position will calculate the twist angle and chord length for optimum performance at that wind speed, rotation speed, and airfoil profile.

6. Part_1.py
  - This is the proof of concept of the lecture notes, uses 1D analysis to find a suitable radius of the blade, then calculates the nodal outputs for a specific radial position on the blade.

7. BEM_500kW.py
  - Uses the lecture notes and variables stated in them to calculate the power output and other plots of a 500kW turbine at varying wind speeds and global pitch angles. (Should be compared to lecture slides to ensure correct calculations)

8. Aerofoil_Ranking.py
  - Used to create a list of the airfoil profiles by ordering by lift to drag coefficient raio (higher better).
  - The profiles for the blade are decided by selecting every 3rd element from this list to create an even profile of airfoils down the length of the blade.

9. BEM_8MW.py, BEM_3MW.py and BEM_15MW.py
  - Uses Aerofoil_Ranking.py and nodal_twist() to create the optimum turbine blade for a given radius and rotation speed at a set wind speed (currently 10 m/s).
  - The created blade profile is used in BEM in the same way as BEM_500kW.py to calculate power output, forces etc.
