function eigval_hho()
    pd = initialize_hho();
    
    for ii = 1:pd.N
        data = compute_cell_data(pd, ii);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the struct with the problem parameters
function pd = initialize_hho()
    pd = struct;
    pd.N        = 16;
    pd.K        = 2;
    pd.gamma0   = 1;
    pd.gamma1   = 1;
    pd.max_stab = 2;
    
    %%%%%%%%%%%%%%%%%%%%%
    pd.h        = 1/pd.N;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quadrature rule to integrate polynomials
function [nodes, weights] = gauss_quadrature(order, h, elem)
    if rem(order, 2) == 0
        order = order+1;
    end
    
    if rem(order,2) == 0
        num_nodes = floor((order+1)/2);
    else
        num_nodes = floor((order+2)/2);
    end
    num_nodes = max(num_nodes, 1);

    if num_nodes == 1
        nodes = h*(elem-0.5);
        weights = h;
        return;
    end
    M = zeros(num_nodes, num_nodes);
    for ii = 1:num_nodes-1
        v = sqrt(1/(4-1/((ii)^2)));
        M(ii+1,ii) = v;
        M(ii,ii+1) = v;
    end

    [d,v] = eig(M);
    tnodes = diag(v)';
    tweights = (d(1,:).^2); 
    
    weights = tweights*h;
    all_ones = ones(size(tnodes));
    nodes = h*( (0.5*tnodes) + ((elem-0.5)*all_ones));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute basis and derivatives
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Util functions
function xT = cell_center(pd, ii)
    xT = (ii-0.5)*pd.h;
end

function [xF1, xF2] = face_centers(pd, ii)
    xF1 = (ii-1)*pd.h;
    xF2 = ii * pd.h;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Precomputed data for cell ii
function cell_data = compute_cell_data(pd, ii)
    [nodes,weights]     = gauss_quadrature(2*(pd.K+1), pd.h, ii);
    
    xT = cell_center(pd, ii);
    MM = zeros(pd.K+2, pd.K+2);
    SM = zeros(pd.K+2, pd.K+2);
    for wi = 1:length(weights)
        phi = compute_basis(xT, nodes(wi), pd.h, pd.K+1, 0);
        dphi = compute_basis(xT, nodes(wi), pd.h, pd.K+1, 1);
        MM = MM + weights(wi)*phi*phi';
        SM = SM + weights(wi)*dphi*dphi';
    end
    
    cell_data           = struct;
    cell_data.ii        = ii;
    cell_data.MM        = MM;
    cell_data.SM        = SM;
    cell_data.qnodes    = nodes;
    cell_data.qweights  = weights;
    cell_data.R0        = reconstruction0(pd, cell_data);
    %cell_data.R1        = reconstruction1(pd, ii, SM);
    cell_data.S0        = stabilization0(pd, cell_data);
    %cell_data.S1        = stabilization1();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = reconstruction0(pd, cd)
    ii          = cd.ii;
    stiff_mat   = cd.SM;
    
    xT = cell_center(pd, ii);
    [xF1, xF2] = face_centers(pd, ii);

    MG = stiff_mat(2:end, 2:end);
    BG(:,1:pd.K+1) = stiff_mat(2:end,1:pd.K+1);
    
    phiF1   = compute_basis(xT, xF1, pd.h, pd.K+1, 0);
    dphiF1  = compute_basis(xT, xF1, pd.h, pd.K+1, 1);
    
    phiF2   = compute_basis(xT, xF2, pd.h, pd.K+1, 0);
    dphiF2  = compute_basis(xT, xF2, pd.h, pd.K+1, 1);
    
    BG(1:end, 1:pd.K+1) = BG(1:end, 1:pd.K+1) + dphiF1(2:end)*phiF1(1:pd.K+1)';
    BG(1:end, 1:pd.K+1) = BG(1:end, 1:pd.K+1) - dphiF2(2:end)*phiF2(1:pd.K+1)';
    BG(1:end, pd.K+2) = - dphiF1(2:end);
    BG(1:end, pd.K+3) = + dphiF2(2:end);
    
    R = MG\BG;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = stabilization0(pd, cd)
    ii          = cd.ii;
    mass_mat    = cd.MM;
    R           = cd.R0;

    xT = cell_center(pd, ii);
    [xF1, xF2] = face_centers(pd, ii);
    
    RR = R(1:end, 1:pd.K+1+2);
    
    M1 = mass_mat(1:pd.K+1,1:pd.K+1);
    M2 = mass_mat(1:pd.K+1,2:pd.K+2);
    P1 = - M1\(M2*RR);
    P1(1:pd.K+1, 1:pd.K+1) = P1(1:pd.K+1, 1:pd.K+1) + eye(pd.K+1);
    
    phiF1 = compute_basis(xT, xF1, pd.h, pd.K+1, 0);
    MFF = 1;
    MFT = phiF1;
    P2 = MFF \ (MFT(2:end)'*RR);
    P2(pd.K+2) = P2(pd.K+2)-1;
    P3 = MFF \ (MFT(1:pd.K+1)'*P1);
    B = P2 + P3;
    
    S = B' * MFF * B / pd.h;
    
    phiF2 = compute_basis(xT, xF2, pd.h, pd.K+1, 0);
    MFF = 1;
    MFT = phiF2;
    P2 = MFF \ (MFT(2:end)'*RR);
    P2(pd.K+3) = P2(pd.K+3)-1;
    P3 = MFF \ (MFT(1:pd.K+1)'*P1);
    B = P2 + P3;
    
    S = S + B' * MFF * B / pd.h;
end



