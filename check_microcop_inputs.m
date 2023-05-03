function [] = check_microcop_inputs(seg, tissue_comps, grid_params, adjustments)
    
    % some checks for common mistakes in MicroCoP usage

    if size(seg.seg,4) > numel(seg.seg_label)
        error('You have to provide a label for each tissue class in the segmentation.')
    elseif size(seg.seg,4) < numel(seg.seg_label)
        warning('There are more tissue class labels than segmentations.')
    end%if
    
    if numel(grid_params.par1.range)~=2
        warning('grid_params.par1.range should have two elements. Using minimum and maximum entry.')
        grid_params.par1.range = [min(grid_params.par1.range), max(grid_params.par1.range)];
    end%if
    if grid_params.mode == 2
        if numel(grid_params.par2.range)~=2
            warning('grid_params.par2.range should have two elements. Using minimum and maximum entry.')
            grid_params.par2.range = [min(grid_params.par2.range), max(grid_params.par2.range)];
        end%if
    end%if
    
    if (grid_params.q < 2) || (round(grid_params.q)~=grid_params.q)
        error('grid_params.q has to be an integer > 1')
    end%if
    
    if ~isfield(grid_params,'spacing')
        switch lower(grid_params.spacing)
            case {'lin','linear','quad','quadratic','log','logarithmic'}
            otherwise
                warning('Grid spacing invalid. Set to linear.')
                grid_params.spacing = 'lin';
        end%switch
        warning('Grid spacing not specified. Set to linear.')
        grid_params.spacing = 'lin';
    end%if
    
    for jj=1:numel(tissue_comps)
        if any(tissue_comps{jj}.par1<=0) || any(tissue_comps{jj}.f<=0)
            error('Tissue component parameters and volume fractions cannot be <= 0.')
        end%if
        if grid_params.mode == 2
            if any(tissue_comps{jj}.par2<=0)
                error('Tissue component parameters cannot be <= 0.')
            end%if
        end%if
        if sum(tissue_comps{jj}.f) ~= 1
            warning(['Tissue component volume fractions have to sum up to 1. Rescaling tissue_comps{',num2str(jj),'}.f'])
            tissue_comps{jj}.f = tissue_comps{jj}.f / sum(tissue_comps{jj}.f);
        end%if
        if any(tissue_comps{jj}.par1 > grid_params.par1.range(2)) || any(tissue_comps{jj}.par1 < grid_params.par1.range(1))
            error('Tissue component parameters have to be within the range specified for the grid.')
        end%if
        if grid_params.mode == 2
            if any(tissue_comps{jj}.par2 > grid_params.par2.range(2)) || any(tissue_comps{jj}.par2 < grid_params.par2.range(1))
                error('Tissue component parameters have to be within the range specified for the grid.')
            end%if
        end%if
        if numel(tissue_comps{jj}.subcomp_label)~=numel(tissue_comps{jj}.par1) || ...
                numel(tissue_comps{jj}.par1)~=numel(tissue_comps{jj}.f)
            error('Tissue components need the same number of parameters, volume fractions and labels (within a tissue class).')
        end%if
        if grid_params.mode == 2
            if numel(tissue_comps{jj}.subcomp_label)~=numel(tissue_comps{jj}.par2)
                error('Tissue components need the same number of parameters, volume fractions and labels (within a tissue class).')
            end%if
        end%if
    end%for
    
    if adjustments.compartmental_gauss_std <= 0 || adjustments.compartmental_gauss_std > 100
        error('adjustments.compartmental_gauss_std must be in ]0,100].')
    end%if
    if adjustments.compartmental_gauss_std_increase < 1
        warning('adjustments.compartmental_gauss_std_increase should be >= 1. Set to 1.')
        adjustments.compartmental_gauss_std_increase = 1;
    end%if
    if adjustments.vxl_vxl_std_par < 0
        warning('adjustments.vxl_vxl_std_par must be > 0. Set to 0.')
        adjustments.vxl_vxl_std_par = 0;
    end%if
    if adjustments.vxl_vxl_std_f < 0
        warning('adjustments.vxl_vxl_std_f must be > 0. Set to 0.')
        adjustments.vxl_vxl_std_f = 0;
    end%if

end%function