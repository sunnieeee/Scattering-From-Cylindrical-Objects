# Scattering from cylindrical objects

## Description

`CEM 2020 Fall final, School of Physics, Huazhong University of Science and Technology`

Use 2D-FDTD to calculate total fields and scattered fields in area where a cylindrical object is scattering in free space with given source, under Absorbing Boundary Conditions (ABC).

## Group member

吕大为，曹玉清

## Setup &  Results

### ./SRC/TFSF+PML3.0.py

- Contributor: 吕大为
- IDE: Spyder
- About: 见 ./DOC/Scattering From Cylindrical Objects一些说明.pdf
- Results

### ./SRC/UPML_Guass_fdtd.m

- Contributor: 曹玉清

- IDE: `MATLAB R2017b`及以上（`F5`运行）

- Set:

  ABC: UPML

  Source: Guassian source 

  Cylindrical object:

  - Radius = 20 (m)

  - Epsilon_r = 6.25

- Results

![CylinScat_UPML_Guass](./PICS/CylinScat_UPML_Guass.gif)

