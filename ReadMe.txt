On my honor I completed the CIS survey for this course.

Built using freeglut and g++ and Eigen on Windows with the following commands: (run in the final folder)

Run the simulation
g++ jsoncpp.cpp Particle.cpp SimulationParameters.cpp FluidSimulator.cpp -lopengl32 -lglew32 -lfreeglut -lglu32 -IC:/Packages/Eigen3

Inputs for the simulation:
./a inputs/fluid.json    

Display the particles at each frame and generate pngs
g++ jsoncpp.cpp Particle.cpp SimulationParameters.cpp ParticleViewer.cpp -lopengl32 -lglew32 -lfreeglut -lglu32 -IC:/Packages/Eigen3

Inputs for the viewer:
./a fluid.%03d.part 

PIC_result.gif is a GIF of using 0, or purely PIC, as the FLIP ratio
result.gif is a GIF of using .95, blended PIC/FLIP, as the default result of the tutorial's given inputs

Basic and Incremental 1-6 contain code for each step of the tutorial as I worked through it from this source: 
https://unusualinsights.github.io/fluid_tutorial/
Final is incremental 7 in the tutorial

final/partsPIC contain the particle output for 0 ratio output
final/parts contain the particle output for .95 ratio output

(Optional)
For earlier incremental folders use this command to run the test:
g++ jsoncpp.cpp Particle.cpp SimulationParameters.cpp StaggeredGridTest.cpp -lopengl32 -lglew32 -lfreeglut -lglu32 -IC:/Packages/Eigen3
