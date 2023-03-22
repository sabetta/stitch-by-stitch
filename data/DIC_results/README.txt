**************************************
README.txt

by Krishma Singal
Georgia Institute of Technology
School of Physics

version 1 
March 21, 2023
**************************************


DATASET INFORMATION

The contents of this dataset describe the acrylic garter fabric sample undergoing uniaxial stretching. The sample was prepped to undergo Digital Image Correlation (DIC) analysis. Data was collected in December 2022.
It is to be used with the (as of file generation date, unpublished) manuscript:
"Programming Mechanics in Knitted Materials, Stitch by Stitch"
by Krishma Singal, Michael S. Dimitriyev, Sarah E. Gonzalez, Sam Quinn, and Elisabetta A. Matsumoto 

Corresponding author: 
Elisabetta A. Matsumoto
Georgia Institute of Technology
School of Physics
e-mail: sabetta@gatech.edu


CONTENTS

The folder "edited_images" contain 37 TIF images of the garter fabric undergoing stretching along its y-axis. They were converted to black and white using the program FIJI. 

The folders "xDataE1" and "yDataE1" contain CSV files tracking the deformation of a garter fabric undergoing uniaxial stretching. Each file corresponds to a different frame of the fabric being stretched. 

MATERIAL INFORMATION

The acrylic yarn is Brava worsted yarn (28455-White) from KnitPicks(™), which is 100% acrylic. The fabric was coated in graphite powder to aid with the tracking. 

ANALYSIS

We used Ncorr, an open source MATLAB software, to perform the DIC analysis. 