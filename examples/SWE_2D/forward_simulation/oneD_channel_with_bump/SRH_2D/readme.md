# Mesh Preparation and SRH-2D Simulation Guide

## 1. Creating Mesh with GMSH: The `rect_channel_with_bump_refined.geo` file is the geometry file.

## 2. Converting to SRH-2D Format using [pyHMT2D](https://github.com/psu-efd/pyHMT2D): run `gen_mesh.py`

## 3. Running SRH-2D Simulation using [pyHMT2D](https://github.com/psu-efd/pyHMT2D): run `run_SRH_2D.py`

## The results are in the generated VTK files, which can be visualized using Paraview. The SRH-2D cases files (`.srhgeom`, `.srhmat`, `.srhhydro`) are also generated and can be used by Hydrograd.
