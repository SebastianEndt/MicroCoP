function [] = plot_microcop(tissue_comps, cvf, seg, phantom, grid_params)
    % plot only one slice
    slice_to_plot = round(size(seg.seg,3)/2);
    
    % sanity check: groundtruth compartmental volume fraction maps
    % subplots: rows: tissue classes, columns: sub-compartments
    %           last column: sub of sub-compartments
    n_subcomps = 0; for j_tissue_class=1:numel(tissue_comps), n_subcomps = max(n_subcomps,numel(tissue_comps{j_tissue_class}.f)); end
    figure();
    for j_tissue_class=1:numel(tissue_comps)
        for j_subcomp=1:numel(tissue_comps{j_tissue_class}.f)
            subplot(numel(tissue_comps),n_subcomps+1,sub2ind([n_subcomps+1,numel(tissue_comps)],j_subcomp,j_tissue_class))
            imagesc(cvf{j_tissue_class}(:,:,slice_to_plot,j_subcomp))
            title([seg.seg_label{j_tissue_class},' ',tissue_comps{j_tissue_class}.subcomp_label{j_subcomp}])
            colorbar;
            axis equal;
            xlim([1,size(phantom,2)]);
            ylim([1,size(phantom,1)]);
            axis off;        
        end% for j_subcomp
        subplot(numel(tissue_comps),n_subcomps+1,sub2ind([n_subcomps+1,numel(tissue_comps)],n_subcomps+1,j_tissue_class))
        imagesc(sum(cvf{j_tissue_class}(:,:,slice_to_plot,:),4))
        title([seg.seg_label{j_tissue_class},' total'])
        colorbar;
        axis equal;
        xlim([1,size(phantom,2)]);
        ylim([1,size(phantom,1)]);
        axis off;
    end% for j_tissue_class
    sgtitle('Compartmental volume fraction maps')

    % sanity check: ROI-average spectra
    figure();
    for j_tissue_class=1:numel(tissue_comps)
        subplot(floor(sqrt(numel(tissue_comps))),ceil(numel(tissue_comps)/floor(sqrt(numel(tissue_comps)))),j_tissue_class)
        if grid_params.mode == 2
            contourf(grid_params.par1.grid,grid_params.par2.grid,masked_average(phantom, seg.seg(:,:,:,j_tissue_class)>0.9*max(col(seg.seg(:,:,:,j_tissue_class)))),10);
            axis square;
            ylabel([grid_params.par2.name,' (',grid_params.par2.unit,')']);
        else
            plot(grid_params.par1.values, masked_average(phantom, seg.seg(:,:,:,j_tissue_class)>0.9*max(col(seg.seg(:,:,:,j_tissue_class)))));
        end%if
        title(seg.seg_label{j_tissue_class});
        xlabel([grid_params.par1.name,' (',grid_params.par1.unit,')']);
    end% for j_tissue_class
    subplot(floor(sqrt(numel(tissue_comps))),ceil(numel(tissue_comps)/floor(sqrt(numel(tissue_comps)))),floor(sqrt(numel(tissue_comps)))*(ceil(numel(tissue_comps)/floor(sqrt(numel(tissue_comps))))))
    if grid_params.mode == 2
        contourf(grid_params.par1.grid,grid_params.par2.grid,masked_average(phantom, sum(seg.seg,4)>0.9),10);
        axis square;
        ylabel([grid_params.par2.name,' (',grid_params.par2.unit,')']);
    else
        plot(grid_params.par1.values, masked_average(phantom, seg.seg(:,:,:,j_tissue_class)>0.9*max(col(seg.seg(:,:,:,j_tissue_class)))));
    end%if
    title('whole phantom');
    xlabel([grid_params.par1.name,' (',grid_params.par1.unit,')']);
    sgtitle('ROI-average spectra');

    % sanity check: parameter maps
    if grid_params.mode == 2
        par1map = sum(phantom,4);
        for jj=1:size(par1map,5)
            par1map(:,:,:,:,jj) = par1map(:,:,:,:,jj) * grid_params.par1.values(jj);
        end% for jj
        par1map = sum(par1map,5);
        par2map = sum(phantom,5);
        for jj=1:size(par2map,4)
            par2map(:,:,:,jj,:) = par2map(:,:,:,jj,:) * grid_params.par2.values(jj);
        end% for jj
        par2map = sum(par2map,4);
    else
        par1map = NaN(size(phantom));
        for jj=1:size(phantom,4)
            par1map(:,:,:,jj) = phantom(:,:,:,jj) * grid_params.par1.values(jj);
        end% for jj
        par1map = sum(par1map,4);
    end%if

    figure();
    subplot(1,grid_params.mode,1)
    imagesc(par1map(:,:,slice_to_plot),[0,max(grid_params.par1.values)]);
    cb = colorbar;
    ylabel(cb,[grid_params.par1.name,' (',grid_params.par1.unit,')']);
    title(['mean ',grid_params.par1.name,' map']);
    axis equal;
    xlim([1,size(phantom,2)]);
    ylim([1,size(phantom,1)]);
    axis off;
    
    if grid_params.mode == 2
        subplot(1,2,2)
        imagesc(par2map(:,:,slice_to_plot),[0,max(grid_params.par2.values)]);
        cb = colorbar;
        ylabel(cb,[grid_params.par2.name,' (',grid_params.par2.unit,')']);
        title(['mean ',grid_params.par2.name,' map']);
        axis equal;
        xlim([1,size(phantom,2)]);
        ylim([1,size(phantom,1)]);
        axis off;
    end%if
end%function