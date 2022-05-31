# SMEARS
## Statistical Modeler for Extended Astrophysical Radio Sources

SMEARS is a modeling package. It takes observations processed through the FHD (Fast Holographic Deconvolution) software package (found here: https://github.com/EoRImaging/FHD) and converts these observations into models for sources. These models consider the potential of these sources to be diffuse and can either create extended models for diffuse sources or recognize and model sources as point-like. The code also handles situations where sources have both diffuse and point like structure, producing an extended model for the diffuse structure and a point model for point-like structure.

SMEARS is made up of several different functions, each designed to make the process of searching for sources, modeling them, and outputting catalog files as intuitive as possible. This is done through several steps:

1. Building a Pandas DataFrame contianing the GLEAM objects from an input GLEAM catlog and the sources identified in the FHD processed observation files. Each row in this DataFrame represents an object from the GLEAM file and it's data from each observation. Objects not observed are indicated by a 0 magnitude.
2. Identify sources that may benefit from a more complete model. This can be done in many ways, but the easiest may be to use the source_search function, which allows a user to specify things like minimum magnitude in Jy, a block of sky to search, and several other parameters that may be helpful to consider.
3. Models can then be generated for such sources. This process is a bit more complex and is detailed below.
4. Writing out a new file in the SkyModel format from the pyradiosky  libraray containing the original GLEAM catalog but with the new sources replacing the old GLEAM objects. This is accomplished with the gen_SkyModel function.

### The Modeling Technique

The process  by which these models are generated takes place through several functions. The model generation sequence can be initiated by calling the gen_SkyModel function and feeding it a list of sources to model or by callin the create_plots function on a single source. Either of these methods will initiate the same sequence of code to generate the source models.

## Installation

## Dependencies

## 
