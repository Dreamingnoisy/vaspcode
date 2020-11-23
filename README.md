# vaspcode
Some scripts to postprocess the vasp data
IF YOU HAVE ANY QUESTIONS, FEEL FREE TO LEAVE YOUR COMMENTS!

## trajectory.py, movie.xyz, rdf_example.py and rdf.png
movie.xyz is the trajectory file from MD(Molecular Dynamics) calculations .  
trajectory.py is the python script to calculate and draw the pair correlation function of two selected elements.  
The rdf_example.py is to demonstrate the usage of the trajectory.py.  
Only support orthogonal cells.  
The definition of pair correlation function (g(r)) can be found in https://doi.org/10.1016/j.jcp.2011.01.048.  
John C. Crocker and Eric R. Weeks also provide useful imformation about g(r) at http://www.physics.emory.edu/faculty/weeks//idl/gofr2.html.  
In https://physicspython.wordpress.com/2019/07/31/radial-distribution-function/, Patrick Gono also wrote a Python program to handle the g(r) 
of O-O pairs at the interface.  
The trajectory.py gives a more convenient way to choose different elements pairs.  
The rdf.png is the image of g(r), it may seem different because the movie.xyz is the trajectory of MD of a crystal alloy.  
Maybe I will use Numpy package to write a more concise script, LOL!  

## XYZ2Data.cpp provided by DingChangJie
A Linux Shell program converting an .xyz file to Lammps-supported data file format.  
(1)To build, simply run build commands like " g++ XYZ2Data.cpp -o XYZ2Data.exe ". A pre-built executable "XYZ2Data.exe" (following the same suffix rule as in Windows) is also provided here.  
(2)Syntax: ./XYZ2Data.exe -i xyz_file_location -o data_file_location -d data_description  
All of the three arguments must not be omitted.  
An example input "Example.xyz" and the resulting data file "Example_converted" are attached. Run the following command to check:  
./XYZ2Data.exe -i Example.xyz -o Example_converted -d DemoFile  
