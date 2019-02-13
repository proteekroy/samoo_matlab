function x = denormalize(xnorm, xmin, xmax)

    x = xmin + xnorm*(xmax-xmin);
    
end