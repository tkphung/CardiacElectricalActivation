# CardiacElectricalActivation
Simple biventricular cardiac electrical activation model using shortest path tree algorithm. This model was used in a study published by [Oomen+ 2021](https://doi.org/10.1007/s10237-021-01532-7).

# Model Creation
The electrical model is based on a canine, biventricular geometry segmented from CINE MRI data. A finite-element mesh representing the left and right ventricular geometries was generated using methods previously published methods ([Herz+ 2005](https://doi.org/10.1007/s10439-005-3312-7), [Phung+ 2020](https://doi.org/10.1115/1.4043876)). The model geometry is stored in the variable `<MODEL>`. This information includes the neighbor connectivity as well as the directions of the myocardial muscle fibers in each element.

# Electrical Activation
The conduction velocities in the model are assigned using one variable to represent the myocardial fiber conduction velocity. The cross-fiber and radial velocities are assumed to be 0.4x the fiber velocity. And the endocardial walls are set to be 6x the fiber velocity. These assumptions of conduction anisotropy are based on work by [Lee+ 2019](https://doi.org/10.1016/j.media.2019.06.017).

The model is solved by setting one (or multiple) initiation elements and propagating the activation to all other elements in the biventricular geometry. This diffusion problem is solved by representing the finite-element geometry as a weighted graph based. The nodes are the centroid of each element, and the edges are weighted by the time it takes to conduct from one node to another. The activation pattern is solved using a shortest path tree algorithm implemented by MATLAB.

The subsequent activation timing is used in a pseudo-ECG formulation to calculate the unipolar lead voltage at locations in space representing the 12-lead ECG on the subject-specific canine torso. The shapes of the pseudo-ECG can be directly compared to recorded data (stored in the files `<DATA_EP_bl.mat>` (baseline) and `<DATA_EP_lbbb.mat>` (left bundle branch block)).

# Scripts
## EModel_RunModel.m
This script runs a single electrical activation model and generates the subsequent pseudo-ECG (12 Lead ECG).
The two model parameters to vary are
* v: the myocardial conduction velocity in the main fiber directions
* STAR: the intiation element(s) and time of initiation (in ms) for electrical activat ion

## EModel_Optimize.m
This script optimizes a single electrical activation model to match data from 12-lead ECG at baseline and after inducing Left Bundle Branch Block.
The optimization first tunes the location of electrical activation initiation. Then the absolute conduction velocity is tuned to match the depolarization time measured from the ECG data.

## MRI_to_Electrical_Pipeline.m
This script shows all of the steps from segmented cardiac MRI (using Medviso's Segment) to simulate an electrical activation model. The script is setup as steps that can be run one-by-one to follow the creation and simulation of the model.
The MRI segmentation files called in Step 1 can be accessed in the data directory and loaded into Medviso's Segment to visualize segmentation.

## EP_Slab_Model.m
This script creates a simple 3D slab for electrical activation simulation. This customizable parameters of this model include its geometry, fiber orientations, and conduction velocities.

[![DOI](https://zenodo.org/badge/312748480.svg)](https://zenodo.org/badge/latestdoi/312748480)
