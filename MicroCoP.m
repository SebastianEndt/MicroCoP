function [phantom, cvf] = MicroCoP(seg, tissue_comps, grid_params, adjustments, show_plots, save_phantom)
    % MicroCoP generates Microstructure Correlation Phantoms
    %
    % Usage:
    % [phantom, cvf] = synth_phantom_gen(seg, tissue_comps, grid_params, adjustments, plot_results, save_phantom)
    %
    % Inputs:
    %   seg             struct. Required fields:
    %                       seg (array of size [Nx,Ny,Nsl,Ntc], containing the 3D segmentations of Ntc tissue classes. Segmentations should sum up to 1.)
    %                       seg_label (cell array with Ntc elements. Each element contains a string to label the corresponding tissue class in seg. Only for plotting.)
    %   tissue_comps    cell array with Ntc elements. Each cell contains a struct with the fields:
    %                       subcomp_label (cell array with Nsc elements,  containing Nsc strings to label the corresponding sub-compartments. Only for plotting.)
    %                       par1 (array with Nsc elements containing the (first) tissue parameter for each of Nsc sub-compartments in the corresponding tissue class. E.g. T1 relaxation times.)
    %                       f (array with Nsc volume fractions for each sub-compartment. Has to sum up to 1.)
    %                       par2 (optional, equal to par1 with a second parameter, needed for 2D correlation spectra. If not specified, 1D spectra with par1 only are generated.)
    %   grid_params     struct. Required fields:
    %                       spacing (string specifying the grid spacing, 'lin' (default), 'quad' or 'log')
    %                       q (integer specifying the grid size (q*1 for 1D or q*q for 2D)
    %                       par1 (struct with fields:
    %                           range (two scalars defining min and max value or par1 in the spectra)
    %                           name (string specifying the parameter, e.g. 'T1', only for plotting.)
    %                           unit (string specifying the unit of par1, e.g. 'ms', only for plotting))
    %                       par2 (equal to par1, specifying a second parameter. Needed for 2D correlation spectra. If not specified, 1D spectra for only par1 are generated.)
    %   adjustments     (optional) struct. Required fields:
    %                       compartmental_gauss_std (scalar specifying the standard deviation of Gaussian peaks in the spectra. (in percent of the parameter range). Minimum 1. Default: 1.);
    %                       compartmental_gauss_std_increase (Scalar specifying how much higher the std of Gaussians at maximum par1/par2 should be be compared to minimum par1/par2. 1 means same std for all peaks. Default: 1)
    %                       vxl_vxl_std_par (scalar specifying the voxel-to-voxel standard deviation of par1/par2 in multiples of the minimum par1/par2 in the grid)
    %                       vxl_vxl_std_f (scalar specifying the voxel-to-voxel standard deviation of f. 1 is the sum of all f in a voxel)
    %   show_plots      (optional) boolean. true=show plots for sanity check, false=don't show (default).
    %   save_phantom    (optional) boolean. true=save phantom data in a unique .mat-file, false=don't save (default).
    %
    % Outputs:
    %   phantom         array of size [Nx,Ny,Nysl,q,q] for 2D or [Nx,Ny,Nsl,q] for 1D containing the full multicomponent spectra for each voxel.
    %   cvf             cell array with Ntc entries. Each cell containing an array of size [Nx,Ny,Nsl,Nsc] with 3D compartmental volume fractions maps for each of Nsc sub-compartments in the corresponding tissue class.
    %
    % (c) Sebastian Endt, Technische Hochschule Ingolstadt, Ingolstadt, Germany 2023
    %
    % Whenever you use MicroCoP, please cite:
    %   Sebastian Endt, Carolin M. Pirkl, Marco Palombo, Marion I. Menzel "MicroCoP: digital Microstructure Correlation Phantoms for benchmarking of multicomponent MRI methods", Proceedings of ISMRM 2023
    
    %% Preparations

    % check optional inputs
    if nargin<6
        save_phantom = false;
        if nargin<5
            show_plots = false;
            if nargin<4
                adjustments.compartmental_gauss_std = 1;
                adjustments.compartmental_gauss_std_increase = 1;
                adjustments.vxl_vxl_std_par = 0.25;
                adjustments.vxl_vxl_std_f = 0.005;
            end%if
        end%if
    end%if
    
    % determine dimensionality
    if isfield(tissue_comps{1},'par2') && isfield(grid_params,'par2')
        grid_params.mode = 2;
    else
        grid_params.mode = 1;
    end%if
    
    % check inputs
    grid_params.par1.range = sort(grid_params.par1.range);
    if grid_params.mode == 2
        grid_params.par2.range = sort(grid_params.par2.range);
    end%if
    check_microcop_inputs(seg, tissue_comps, grid_params, adjustments);

    % make sure segmentations sum up to 1 in every voxel
    seg.seg = seg.seg./repmat(sum(seg.seg,4),[1,1,1,size(seg.seg,4)]);
    seg.seg(isnan(seg.seg))=0;

    % parametrize grid
    grid_params.par1.range = sort(grid_params.par1.range);
    if grid_params.mode == 2
        grid_params.par2.range = sort(grid_params.par2.range);
    end%if
    switch lower(grid_params.spacing)
        case 'lin'
            grid_params.par1.values = linspace(grid_params.par1.range(1),grid_params.par1.range(2),grid_params.q);
            if grid_params.mode == 2
                grid_params.par2.values = linspace(grid_params.par2.range(1),grid_params.par2.range(2),grid_params.q);
            end%if
        case 'quad'
            grid_params.par1.values = quadspace(grid_params.par1.range(1),grid_params.par1.range(2),grid_params.q);
            if grid_params.mode == 2
                grid_params.par2.values = quadspace(grid_params.par2.range(1),grid_params.par2.range(2),grid_params.q);
            end%if
        case 'log'
            grid_params.par1.values = logspace_SE(grid_params.par1.range(1),grid_params.par1.range(2),grid_params.q);
            if grid_params.mode == 2
                grid_params.par2.values = logspace_SE(grid_params.par2.range(1),grid_params.par2.range(2),grid_params.q);
            end%if
    end%switch
    if grid_params.mode == 2
        [grid_params.par1.grid, grid_params.par2.grid] = meshgrid(grid_params.par1.values, grid_params.par2.values);
    end%if
    
    %% Create the gaussian distributions for each compartment in 1D (for 
    % par1 and par2). stack the par1 distributions along the par2 dimension
    % and vice versa. then multiply them element-wise to obtain 2D
    % distributions. normalize the 2D distributions to a sum of 1

    % centers of spectra components for all tissue classes and their 
    % sub-compartments. value differs from voxel to voxel
    comp_center_std_par1 = adjustments.vxl_vxl_std_par*min(grid_params.par1.range); % std of the mean par1 of each sub-compartment (voxel to voxel comparison)
    comp_means_par1 = cell(size(tissue_comps));
    if grid_params.mode == 2
        comp_center_std_par2 = adjustments.vxl_vxl_std_par*min(grid_params.par2.range); % std of the mean par2 of each sub-compartment (voxel to voxel comparison)
        comp_means_par2 = cell(size(tissue_comps));
    end%if
    for j_tissue_class = 1:numel(tissue_comps)
        comp_means_par1{j_tissue_class} = comp_center_std_par1*randn([size(seg.seg,1:3),numel(tissue_comps{j_tissue_class}.par1)]) + ...
            repmat(permute(tissue_comps{j_tissue_class}.par1,[1,3,4,2]),[size(seg.seg,1:3),1]);
        comp_means_par1{j_tissue_class} = max(min(comp_means_par1{j_tissue_class},...
            grid_params.par1.range(2)*ones(size(comp_means_par1{j_tissue_class}))), ...
            grid_params.par1.range(1)*ones(size(comp_means_par1{j_tissue_class})));
        if grid_params.mode == 2
            comp_means_par2{j_tissue_class} = comp_center_std_par2*randn([size(seg.seg,1:3),numel(tissue_comps{j_tissue_class}.par2)]) + ...
                repmat(permute(tissue_comps{j_tissue_class}.par2,[1,3,4,2]),[size(seg.seg,1:3),1]);
            comp_means_par2{j_tissue_class} = max(min(comp_means_par2{j_tissue_class},...
                grid_params.par2.range(2)*ones(size(comp_means_par2{j_tissue_class}))), ...
                grid_params.par2.range(1)*ones(size(comp_means_par2{j_tissue_class})));
        end%if
    end%for j_tissue_class
    clearvars comp_center_std_par1 comp_center_std_par2

    % volume fraction of each sub-compartment. value differs from voxel to 
    % voxel. sum of all volume fractions has to be 1
    comp_means_f = cell(size(tissue_comps));
    for j_tissue_class = 1:numel(tissue_comps)
        comp_means_f{j_tissue_class} = adjustments.vxl_vxl_std_f*randn([size(seg.seg,1:3),numel(tissue_comps{j_tissue_class}.f)]).*...
            repmat(permute(tissue_comps{j_tissue_class}.f/min(tissue_comps{j_tissue_class}.f),[1,3,4,2]),[size(seg.seg,1:3),1]) + ...
            repmat(permute(tissue_comps{j_tissue_class}.f,[1,3,4,2]),[size(seg.seg,1:3),1]);
        comp_means_f{j_tissue_class} = max(comp_means_f{j_tissue_class},zeros(size(comp_means_f{j_tissue_class})));
        comp_means_f{j_tissue_class} = comp_means_f{j_tissue_class}./repmat(sum(comp_means_f{j_tissue_class},4),[1,1,1,size(comp_means_f{j_tissue_class},4)]);
    end%for j_tissue_class

    % stds for all Gaussians
    comp_dist_std_par1 = range(grid_params.par1.values)/100*adjustments.compartmental_gauss_std; % std of the gaussians in the spectra
    if grid_params.mode == 2
        comp_dist_std_par2 = range(grid_params.par2.values)/100*adjustments.compartmental_gauss_std; % std of the gaussians in the spectra
    end%if
    comp_dist_stds_par1 = cell(size(tissue_comps));
    if grid_params.mode == 2
        comp_dist_stds_par2 = cell(size(tissue_comps));
    end%if
    for j_tissue_class=1:numel(tissue_comps)
        comp_dist_stds_par1{j_tissue_class} = comp_dist_std_par1*(tissue_comps{j_tissue_class}.par1/max(grid_params.par1.values)*(adjustments.compartmental_gauss_std_increase-1)+1);
        if grid_params.mode == 2
            comp_dist_stds_par2{j_tissue_class} = comp_dist_std_par2*(tissue_comps{j_tissue_class}.par2/max(grid_params.par2.values)*(adjustments.compartmental_gauss_std_increase-1)+1);
        end%if
    end%for j_tissue_class
    clearvars comp_dist_std_par1 comp_dist_std_par2

    % create actual Gaussians
    if grid_params.mode == 2
        gaussians_1D_par1 = cell(size(tissue_comps));
        gaussians_1D_par2 = cell(size(tissue_comps));
        gaussians_2D = cell(size(tissue_comps));
    else
        gaussians_1D = cell(size(tissue_comps));
    end%if
    for j_tissue_class=1:numel(tissue_comps)
        if grid_params.mode == 2
            gaussians_1D_par1{j_tissue_class} = NaN([size(seg.seg,1:3), numel(grid_params.par2.values), numel(grid_params.par1.values), size(comp_means_par1{j_tissue_class},4)]);
            gaussians_1D_par2{j_tissue_class} = NaN(size(gaussians_1D_par1{j_tissue_class}));
            gaussians_2D{j_tissue_class} = NaN(size(gaussians_1D_par1{j_tissue_class}));
        else
            gaussians_1D{j_tissue_class} = NaN([size(seg.seg,1:3), numel(grid_params.par1.values), size(comp_means_par1{j_tissue_class},4)]);
        end%if
        for j_subcomp = 1:size(comp_means_par1{j_tissue_class},4)
            if grid_params.mode == 2
                gaussians_1D_par1{j_tissue_class}(:,:,:,:,:,j_subcomp) = repmat(permute(stack_of_gaussians(grid_params.par1.values, comp_means_par1{j_tissue_class}(:,:,:,j_subcomp), comp_dist_stds_par1{j_tissue_class}(j_subcomp)),[1,2,3,5,4,6]),[1,1,1,numel(grid_params.par2.values)]);
                gaussians_1D_par2{j_tissue_class}(:,:,:,:,:,j_subcomp) = repmat(stack_of_gaussians(grid_params.par2.values, comp_means_par2{j_tissue_class}(:,:,:,j_subcomp), comp_dist_stds_par2{j_tissue_class}(j_subcomp)),[1,1,1,1,numel(grid_params.par1.values)]);
                gaussians_2D{j_tissue_class}(:,:,:,:,:,j_subcomp) = gaussians_1D_par1{j_tissue_class}(:,:,:,:,:,j_subcomp).*gaussians_1D_par2{j_tissue_class}(:,:,:,:,:,j_subcomp);
                gaussians_2D{j_tissue_class}(:,:,:,:,:,j_subcomp) = gaussians_2D{j_tissue_class}(:,:,:,:,:,j_subcomp)./repmat(sum(gaussians_2D{j_tissue_class}(:,:,:,:,:,j_subcomp),4:5),[1,1,1,numel(grid_params.par2.values),numel(grid_params.par1.values)]);
            else
                gaussians_1D{j_tissue_class}(:,:,:,:,j_subcomp) = permute(stack_of_gaussians(grid_params.par1.values, comp_means_par1{j_tissue_class}(:,:,:,j_subcomp), comp_dist_stds_par1{j_tissue_class}(j_subcomp)),[1,2,3,5,4,6]);
            end%if
        end%for j_subcomp
    end%for j_tissue_class
    clearvars gaussians_1D_par1 gaussians_1D_par2 comp_dist_stds_par1 comp_dist_stds_par2 comp_means_par1 comp_means_par2

    %% Weight the sub-compartments according to tissue_comps and sum them up
    if grid_params.mode == 2
        spectra = single(NaN([size(gaussians_2D{1},1:5),numel(tissue_comps)]));
        for j_tissue_class=1:numel(tissue_comps)
            spectra(:,:,:,:,:,j_tissue_class) = single(sum(gaussians_2D{j_tissue_class}.*repmat(permute(comp_means_f{j_tissue_class},[1,2,3,5,6,4]),[1,1,1,numel(grid_params.par1.values),numel(grid_params.par2.values),1]),6));
            spectra(:,:,:,:,:,j_tissue_class) = spectra(:,:,:,:,:,j_tissue_class)./repmat(sum(spectra(:,:,:,:,:,j_tissue_class),4:5),[1,1,1,numel(grid_params.par1.values),numel(grid_params.par2.values)]);
        end%for j_tissue_class
    else
        spectra = single(NaN([size(gaussians_1D{1},1:4),numel(tissue_comps)]));
        for j_tissue_class=1:numel(tissue_comps)
            spectra(:,:,:,:,j_tissue_class) = single(sum(gaussians_1D{j_tissue_class}.*repmat(permute(comp_means_f{j_tissue_class},[1,2,3,5,4]),[1,1,1,numel(grid_params.par1.values),1]),5));
            spectra(:,:,:,:,j_tissue_class) = spectra(:,:,:,:,j_tissue_class)./repmat(sum(spectra(:,:,:,:,j_tissue_class),4),[1,1,1,numel(grid_params.par1.values)]);
        end%for j_tissue_class
    end%if
    clearvars gaussians_1D gaussians_2D

    %% Generate compartmental volume fraction (CVF) maps
    cvf = cell(size(tissue_comps));
    for j_tissue_class=1:numel(tissue_comps)
        cvf{j_tissue_class} = single(NaN(size(comp_means_f{j_tissue_class})));
        for j_subcomp = 1:size(comp_means_f{j_tissue_class},4)
            cvf{j_tissue_class}(:,:,:,j_subcomp) = comp_means_f{j_tissue_class}(:,:,:,j_subcomp).*seg.seg(:,:,:,j_tissue_class);
        end%for j_subcomp
    end%for j_tissue_class
    clearvars j_tissue_class j_subcomp comp_means_f

    %% Generate the actual phantom by summing up the spectra according to
    % the tissue-class segmentations
    if grid_params.mode == 2
        phantom = sum(spectra.*repmat(permute(seg.seg,[1,2,3,5,6,4]),[1,1,1,numel(grid_params.par2.values),numel(grid_params.par1.values),1]),6);
        phantom = single(phantom)./repmat(sum(single(phantom),4:5),[1,1,1,numel(grid_params.par1.values),numel(grid_params.par2.values)]);
    else
        phantom = sum(spectra.*repmat(permute(seg.seg,[1,2,3,5,4]),[1,1,1,numel(grid_params.par1.values),1]),5);
        phantom = single(phantom)./repmat(sum(single(phantom),4),[1,1,1,numel(grid_params.par1.values)]);
    end%if
    phantom(isnan(phantom)) = 0;
    clearvars spectra

    %% Plots for sanity checks
    if show_plots
        plot_microcop(tissue_comps, cvf, seg, phantom, grid_params);
    end%if show_plots

    %% Save phantom
    if save_phantom
        save_struct.phantom = phantom;
        save_struct.cvf = cvf;
        save_struct.seg = seg;
        save_struct.tissue_comps = tissue_comps;
        save_struct.grid_params = grid_params;
        save_struct.adjustments = adjustments;
        save([datestr(now,'yymmdd_HHMMSS_'),'MicroCoP_phantom.mat'],'-struct','save_struct')
    end%if save_phantom
end%function