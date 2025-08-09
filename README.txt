Jonas Tividor
CMSC435
Mon/Wed 1:00-2:15pm
Project 2 - Ray Tracing II

IMPORTANT NOTICES:
I own a windows machine and had a lot of trouble getting Eigen to set up. So my include path uses the entire enclude path on my machine.
THIS WILL NOT WORK ON OTHER MACHINES.
You might need to replace the #include <C:/Toolbox/Eigen/...> to just Eigen/Dense in order to run my code.

Extra Credit Smooth Shading was completed: ./shade -p {input} {output}

To Run program:
Type "make" and an executable should be created
To run this executable just type "./shade {input} {output}"

This program takes in an nff file and generates polygons and ray traces to calculate the colors of pixels on the screen.
Certain mathmatical equations are used to convert between camera, world, and image coordinates. The expected image, as extra credit
was not completed, should be a red stack of triangles on a blue background. For this project, lighting and reflections are considered.

Resources Used:
Bargteil Lectures
Eigen Library Look-up
"Ray Tracing in One Weekend C++ Tutortial" - Avery on Youtube

(I did not do this in a weekend)

For issues with compiling or running contact @
jonast1@umbc.edu