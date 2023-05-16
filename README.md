# Sustainable-Turbine

All scripts should be able to run straight after placing all files into same folder.
They call data from the zipfile and generate their own files to store plots and csv data.

**Files in order of use:**  
**1. Shape_Plotting.py**
  - Used to generate plots of the aerofoil shape profiles, and save them to the "Shape Profiles" file.

**2. Interpolation_Tests.py**
  - Just some initial tests on the type of interpolation to use, compares different methods against each other.

**3. CLD_Profile_Plotting.py**
  - Used to generate plots of the aerofoil CLD profiles, and save them to the "CLD Profiles" file.
  - Also used to generate plots fo the aerofil CLD profiles compared to their interpolated version, and save them to the "Compare CLD Profiles" file.

**4. CLD+Shape_Plotting.py**
  - Used to generate plots of the aerofoil shape and CLD profiles side by side on one plot and save them to the "CLD + Shape Profiles" file.

**5. Functions.py**
  - This contains the functions used by the other calculation codes.
    - cld_func() creates functions that allow the interpolation of Cl and Cd values from a given CLD file from the Aerofoil-data.zip file.
    - nodal() takes the inputs for the wind turbine at a radial position and outputs the nodal values including nodal forces.
    - forces() integrates between the segments to find the torque and normal force for each segment of the blade.
    - nodal_twist() for a given wind speed, nodal position and chord length will calculate the twist angle for optimu performance at that wind speed, rotation speed, and airfoil profile
    - nodal_chord() for a given wind speed and nodal position will calculate the twist angle and chord length for optimum performance at that wind speed, rotation speed, and airfoil profile.

**6. Part_1.py**
  - This is the proof of concept of the lecture notes, uses 1D analysis to find a suitable radius of the blade, then calculates the nodal outputs for a specific radial position on the blade.

**7. BEM_500kW.py**
  - Uses the lecture notes and variables stated in them to calculate the power output and other plots of a 500kW turbine at varying wind speeds and global pitch angles. (Should be compared to lecture slides to ensure correct calculations)

**8. Aerofoil_Ranking.py**
  - Used to create a list of the airfoil profiles by ordering by lift to drag coefficient raio (higher better).
  - The profiles for the blade are decided by selecting every 3rd element from this list to create an even profile of airfoils down the length of the blade.

**9. nodal_chords_testing.py**
  - Used to check that the chord length and twist angle calculation methods still allow the interation towards a constant value using BEM
  - (Shows that they do still iterate towards stable value)

**10. BEM_8MW.py, BEM_3MW.py and BEM_15MW.py**
  - Uses Aerofoil_Ranking.py and nodal_twist() to create the optimum turbine blade for a given radius and rotation speed at a set wind speed (currently 10 m/s).
  - The created blade profile is used in BEM in the same way as BEM_500kW.py to calculate power output, forces etc.

**11. BEM_xMW.py**
  - Is essentially the same as BEM_8MW, 3MW and 15MW but placed there with the intention of changing the desired power output or other conditions to generate whatever turbine is desired using the aerofoil geometry found for Iteraion_7_1(Final)

**Iteration Files:**  
Demonstrate the progression of initial designs  
  
**1. Iteration_1.py**
   - Uses cld data of lecture blade
   - Uses twist angle from lecture notes
   - Uses radius positions and chord length from lecture notes multiplied by 4

**2. Iteration_2.py**
  - Uses blade profile similar to lecture one (38, highest Cl/Cd airfoil provided)
  - Uses twist angle from lecture notes
  - Uses chord length from lecture notes multiplied by 4
  - Uses linespace for radial position from 4.5 - 84.5 m

**3. Iteration_3.py**
  - Uses blade profiles ranked from low Cl/Cd to high, using every 3rd profile
  - Uses twist angle from lecture notes
  - Uses chord length for lecture notes multiplied by 4
  - Uses lisnapce for radial positions from 4.5 - 84.5 m

**4. Iteration_4.py**
  - Uses blade profiles ranked from low Cl/Cd to high, using every 3rd profile
  - Calculates twist angle from optimum angle of attack for airfoil using nodal_twist function
  - Uses chord length for lecture notes multiplied by 4
  - Uses lisnapce for radial positions from 4.5 - 84.5 m

**5. Iteration_5.py**
  - Uses blade profiles ranked from low Cl/Cd to high, using every 3rd profile
  - Calculates twist angle and chord length using nodal_chord function (method 1 - https://www.ehow.co.uk/how_7697179_calculate-along-wind-turbine-blade.html)
  - Uses chord length for lecture notes multiplied by 4
  - Uses lisnapce for radial positions from 4.5 - 84.5 m  

**5.1. Iteration_5_.py**
  - Uses the same method but corrects chord length from going to infinity by mirroring about the third node
  - Produces higher forces than method 3

**6. Iteration_6.py**
  - Uses blade profiles ranked from low Cl/Cd to high, using every 3rd profile
  - Calculates twist angle and chord length using nodal_chord function (method 2 - https://www.mdpi.com/1996-1073/13/9/2320)
  - Uses chord length for lecture notes multiplied by 4
  - Uses lisnapce for radial positions from 4.5 - 84.5 m  

**6.1. Iteration_6_1.py**
  - Uses the same method but corrects chord length from going to infinity by mirroring about the third node
  - Produces an odd blade profile (starts thin, gets thicker, gets thin again)

**7. Iteration_7.py**
  - Uses blade profiles ranked from low Cl/Cd to high, using every 3rd profile
  - Calculates twist angle and chord length using nodal_chord function (method 3 - https://ieeexplore.ieee.org/abstract/document/7884538)
  - Uses chord length for lecture notes multiplied by 4
  - Uses lisnapce for radial positions from 4.5 - 84.5 m  

**7.1. Iteration_7_1(Final).py**
  - Uses the same method but corrects chord length from going to infinity by mirroring about the third node
  - Produces lower forces than method 1 and produces a sensible looking blade profile
  - Chosen "best" method

**7.2. Iteration_7_2.py**
  - Everything is the same as in Iteration_7_1.py however the profiles selected are (hopefully) better.
  - This design is not final as it was created too late into the design process to implement.

**8. Iteration_8.py**
  - Is a silly test, uses all but the cylindrical profile (50 nodes and airfoils) to generate a blade profile and power output
  - Uses method 3 - https://ieeexplore.ieee.org/abstract/document/7884538
