# File Format

The `.h5` files in this directory contain experimental data from uniaxial stretching experiments, and are readable by the Mathematica notebook `template.nb`. Alternatively, their contents are fully described below:

The samples used in each experiment are described in the file name, and in the `description` attribute of `/` in the `.h5` file. Each file contains data derived from up to 5 (five) videos. Not every frame of every video has data available, so the following information is stored for each frame that does have data associated with it:

1. The frame number, in the `/video#/frames` record.
2. The (x,y) position of every tracked pin, in the `/video#/traj` record.
3. The force applied to the sample, in the `/video#/forces` record.

The order in which the tracked pins are stored is by convention:

1. Right pin
2. Left pin
3. Bottom pin (if present)
4. Top pin (if present)
5. Bottom clamp 
6. Top clamp 

Additionally, two properties of the sample are stored in their own record, and are the same across all videos in a file:

1. The distance from one edge of the clamp to the other, measured in pixels, in `/clampSize`. This is used to set the camera scale.
2. The width of the fabric at the clamp, in pixels, in `/fabricSize`. This is used to convert force to stress.


Data Note:

1. All of the pearlized cotton and acrylic samples have 6 points tracked except for pearlized cotton Rib y which only has 4 (does not include the top and bottom pin). 
2. All of the glove prototypes samples have 4 points tracked (no top or bottom pin).


