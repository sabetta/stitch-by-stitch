# uniaxial

This directory contains files relating to processing videos of uniaxial stretching experiments into force-vs-extension data.
The following scripts were run in order for each sample to produce the data in `experimentData` from the source videos:

- `measureFabric.m` digitally measures a reference object in the video to determine the image scale. It also measures the width of the fabric at the clamp to aid in determining the stress from the force.
- `uniaxial.m` and `forceVerify.m` use computer vision methods, combined with automatic and manual filtering to extract force data from the recorded screen of the force gauge, and the positions of the pins as functions of time.
- `viewResults.m` displays and exports the data into the `.h5` format in `experimentData`. See `experimentData/readme.md` for details on the file format

Once the `.h5` files are produced, the Mathematica notebook `template.nb` can be used to read the data and to output it to usable formats (such as CSV).
