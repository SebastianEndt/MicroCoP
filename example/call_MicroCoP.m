%% Parametrize and call MicroCoP
addpath('..'); % this example only works with the example-directory as pwd

%% Outputs
show_plots = true;
save_phantom = false;

%% Tissue class segmentations
% This segmentation is a modified, low-resolution version of a BrainWeb phantom with severe MS lesions (https://brainweb.bic.mni.mcgill.ca/anatomic_ms3.html).
% Size is (nx=90, ny=90, nsl=7, ntc=4) - so there are ntc=4 tissue classes available. seg_label gives the corresponding labels "NAWM", "GM", "CSF" and "lesion".
% Credits for BrainWeb go to the McConnell Brain Imaging Centre (BIC):
%   C.A. Cocosco, V. Kollokian, R.K.-S. Kwan, A.C. Evans: "BrainWeb: Online Interface to a 3D MRI Simulated Brain Database" NeuroImage, vol.5, no.4, part 2/4, S425, 1997.
%   R.K.-S. Kwan, A.C. Evans, G.B. Pike: "MRI simulation-based evaluation of image-processing and classification methods" IEEE Transactions on Medical Imaging. 18(11):1085-97, Nov 1999.
%   R.K.-S. Kwan, A.C. Evans, G.B. Pike: "An Extensible MRI Simulator for Post-Processing Evaluation" Visualization in Biomedical Computing (VBC'96). Lecture Notes in Computer Science, vol. 1131. Springer-Verlag, 1996. 135-140.
%   D.L. Collins, A.P. Zijdenbos, V. Kollokian, J.G. Sled, N.J. Kabani, C.J. Holmes, A.C. Evans: "Design and Construction of a Realistic Digital Brain Phantom" IEEE Transactions on Medical Imaging, vol.17, No.3, p.463--468, June 1998.
seg = load('seg_brainweb_MS_mod.mat','seg','seg_label'); % struct "seg" needs the fields "seg" and "seg_label"
%seg.seg = seg.seg(:,:,4,:); % take only one or few slices to save memory and accelerate computation time.

%% Tissue sub-compartments for all tissue classes
% if you don't specify par2, a phantom with 1D spectra will be generated
% NAWM with sub-compartments MWF, IC, EC
tissue_comps{1}.subcomp_label = {'MWF', 'IC','EC'};
tissue_comps{1}.par1 = [120, 950, 700];
tissue_comps{1}.par2 = [20, 100, 75];
tissue_comps{1}.f = [20, 42.5, 37.5]/100; % should sum up to 1

% GM with sub-compartments MWF, IC, EC
tissue_comps{2}.subcomp_label = {'MWF', 'IC','EC'};
tissue_comps{2}.par1 = tissue_comps{1}.par1 * 1.20;
tissue_comps{2}.par2 = tissue_comps{1}.par2 * 1.05;
tissue_comps{2}.f = [5, 42.5, 52.5]/100; % should sum up to 1

% CSF with sub-compartment CSF
tissue_comps{3}.subcomp_label = {'CSF'};
tissue_comps{3}.par1 = 2750;
tissue_comps{3}.par2 = 290;
tissue_comps{3}.f = 1; % should sum up to 1

% Lesion with sub-compartments IC, EC
tissue_comps{4}.subcomp_label = {'IC','EC'};
tissue_comps{4}.par1 = [tissue_comps{1}.par1(2),tissue_comps{2}.par1(3)];
tissue_comps{4}.par2 = [tissue_comps{1}.par2(2),tissue_comps{2}.par2(3)];
tissue_comps{4}.f = [tissue_comps{1}.f(2), 1-tissue_comps{1}.f(2)]; % should sum up to 1

%% Parametrize grid
% Here: 2D 60x60 grid with linear spacing.
% If you don't specify par2, a phantom with 1D spectra will be generated
grid_params.spacing = 'lin'; % 'lin' (default), 'quad', 'log'
grid_params.q = 60;

grid_params.par1.name = 'T1';
grid_params.par1.range = [50, 3000];
grid_params.par1.unit = 'ms';

grid_params.par2.name = 'T2';
grid_params.par2.range = [5, 300];
grid_params.par2.unit = 'ms';

%% Adjust details
adjustments.compartmental_gauss_std = 1; % How high should the std of gaussians in the spectrum generally be (in percent of the parameter range)? Scalar > 0.
adjustments.compartmental_gauss_std_increase = 1; % How much higher should the std of Gaussians at maximum par1/par2 be compared to minimum par1/par2? (scalar >= 1)
adjustments.vxl_vxl_std_par = 0.25; % How high should the parameters vary from voxel to voxel? (specify standard deviation in multiples of min(grid_params.par1/2.range) ) 
adjustments.vxl_vxl_std_f = 0.005; % How high should the parameters vary from voxel to voxel? (specify standard deviation. 1 is the sum of all f in a voxel)

%% Generate phantom
[phantom, cvf] = MicroCoP(seg, tissue_comps, grid_params, adjustments, show_plots, save_phantom);
