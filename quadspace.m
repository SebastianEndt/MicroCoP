function vec = quadspace(start_value, end_value, steps)
    vec = linspace(sqrt(start_value),sqrt(end_value),steps).^2;
end%function