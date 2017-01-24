function phi = compute_basis(center, point, h, degree, deriv_order)
    ex = (point-center)/h;  % evaluation point
    
    phi = zeros(degree+1,1);
    
    switch deriv_order
        case 0
            for ii = 0:degree
                phi(ii+1) = ex^ii;
            end
            
        case 1
            for ii = 1:degree
                phi(ii+1) = (ii/h)*(ex^(ii-1));
            end
            
        case 2
            for ii = 2:degree
                phi(ii+1) = (ii*(ii-1)/(h*h))*(ex^(ii-2));
            end
            
        case 3
            for ii = 3:degree
                phi(ii+1) = (ii*(ii-1)*(ii-2)/(h*h*h))*(ex^(ii-3));
            end
            
        otherwise
            disp('Basis: derivatives > 3rd are not implemented');
    end
end