function gaussstack = stack_of_gaussians(x, means3D, stds3D)

    if numel(stds3D) == 1
        stds3D = stds3D*ones(size(means3D));
    end%if

    x4D = repmat(permute(x(:),[2,3,4,1]), [size(means3D), 1]);
    means4D = repmat(means3D, [1,1,1,size(x4D,4)]);
    stds4D = repmat(stds3D, [1,1,1,size(x4D,4)]);
    
    gaussstack = exp( (-(x4D-means4D).^2) ./ (2*stds4D.^2));

end%function