function F_av = masked_average(F, mask)
    F_av = squeeze(mean(F.*repmat(mask,[1,1,1,size(F,4),size(F,5)]),1:3));
end%function