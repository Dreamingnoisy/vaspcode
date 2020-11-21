# vaspcode
Some scripts to postprocess the vasp data

## trajectory.py, movie.xyz, rdf_example.py and rdf.png
movie.xyz is the trajectory file from XDATCAR.  
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

IF YOU HAVE ANY QUESTIONS, FEEL FREE TO LEAVE YOUR COMMENTS!
