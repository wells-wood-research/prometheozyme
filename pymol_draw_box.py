"""
Author: https://gist.github.com/ros-luc/c92acd623b8d4fd19242b61d7828359c

This script adds a function in PyMol that allows it to draw boxes.
You can specify the box center and its size. Very useful to represent docking boxes.

Usage:
1) In PyMol, go to File > Run script > this script
2) In the console, type:
   draw_box center=(1, 2, 3), size=(20, 20, 20), spacing=1, linewidth=3, color=(255, 255, 255)

Center, size and color (RGB) must be written as tuples.
Spacing is set to 1 as default, as used by AutodockVina, but can be changed to e.g. 0.375.
Linewidth is not appreciated in the viewing window, but will appear correctly once the PNG is rendered.

Adapted from https://pymolwiki.org/index.php/DrawBoundingBox

"""
from pymol.cgo import *
from pymol import cmd
from random import randint


def draw_box(center=(0, 0, 0),
             size=(30, 30, 30),
             spacing=1,
             color=(255, 255, 255),
             linewidth=3
             ):
    """
       4-----------8
      /|          /|
     / |         / |
    3-----------7  |
    |  |    º   |  |
    Y  2--------|--6
    | /         | Z
    |/          |/
    1---X-------5
	This will draw a box using the numbered vertices and with the specified center
    """
    if not type(center) == tuple:
        center = eval(center)
    if not type(size) == tuple:
        size = eval(size)
    if not type(color) == tuple:
        color = eval(color)
	
	# Use scaling
    x = float(spacing) * float(size[0])
    y = float(spacing) * float(size[1])
    z = float(spacing) * float(size[2])

	# Using the center, locate the min and max for each axis
    minX = float(center[0]) - (x/2)
    minY = float(center[1]) - (y/2)
    minZ = float(center[2]) - (z/2)
    maxX = float(center[0]) + (x/2)
    maxY = float(center[1]) + (y/2)
    maxZ = float(center[2]) + (z/2)

	# Each point has three coordinates, refer to drawing above
    v = {1: (minX, minY, minZ),
         2: (minX, minY, maxZ),
         3: (minX, maxY, minZ),
         4: (minX, maxY, maxZ),
         5: (maxX, minY, minZ),
         6: (maxX, minY, maxZ),
         7: (maxX, maxY, minZ),
         8: (maxX, maxY, maxZ),
         }

	# These are the pairs of vertices that will make the box
    points = [v[1], v[2],
              v[1], v[3],
              v[1], v[5],
              v[2], v[4],
              v[2], v[6],
              v[3], v[4],
              v[3], v[7],
              v[4], v[8],
              v[5], v[6],
              v[5], v[7],
              v[6], v[8],
              v[7], v[8]
              ]
	
	# Set the color
    r = float(color[0] / 255)
    g = float(color[1] / 255)
    b = float(color[2] / 255)

	# Initialize box
    box = [
        LINEWIDTH, float(linewidth),
        BEGIN, LINES,
        COLOR, r, g, b
    ]
	
	# Add all vertices to the box
    for point in points:
        box.append(VERTEX)
        for coordinate in point:
            box.append(float(coordinate))
    box.append(END)

	# Finally, give the object a name and load it
    box_name = "box_" + str(randint(0, 10000))
    while box_name in cmd.get_names():
        box_name = "box_" + str(randint(0, 10000))

    cmd.load_cgo(box, box_name)


cmd.extend("draw_box", draw_box)