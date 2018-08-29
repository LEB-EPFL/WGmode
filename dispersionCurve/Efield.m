function [out] = Efield(x, corDim, ki, sigma, phi, xi) 

    % Electric field distribution: TE mode profile in the plane perpendicular
    % to the propagation direction
    if x > corDim/2
        out = cos(ki*corDim/2 - phi)*exp(-sigma*(x - corDim/2));
    elseif x >= - corDim/2 && x <= corDim/2
        out = cos(ki*x - phi);
    else
        out = cos(ki*corDim/2 + phi)*exp(xi*(x + corDim/2));
    end

end

