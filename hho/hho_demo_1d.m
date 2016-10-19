%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       /\
%      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
%     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
%    /\    /\
%   /__\  /__\    1D Hybrid High-order demo code
%  /_\/_\/_\/_\
%

function hho()

    what = 2;

    if (what == 1)
        plot_hho_convergence();
    end
    
    N = 16;
    K = 2;
    pd = initialize_hho(N, K);
    
    if (what == 2)
        test_gradient_reconstruction(pd)
    end
        
    if (what == 3)
        [~, x, y] = hho(pd, @f_load, @f_solution);
        plot(x,y);
    end
    
    if (what == 4)
        hho_eigenvalues(pd);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the struct with the problem parameters
function pd = initialize_hho(N, K)
    pd = struct;
    pd.N = N;
    pd.K = K;
    pd.h = 1/pd.N;
    pd.tp = 8;
end

% Load function
function l = f_load(x)
    l = sin(pi*x) * pi^2;
    return;
end

% Analytical solution
function l = f_solution(x)
    l = sin(pi*x);
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_hho_convergence()
    
    for deg = 0:4
        maxii = 11;
        nodes = 4;
        for ii = 1:maxii
            disp(sprintf('**** N = %d, K = %d ****', nodes, deg));
            pd = initialize_hho(nodes, deg);
            [e, ~, ~, he] = hho(pd, @f_load, @f_solution);
            err(deg+1, ii) = e;
            herr(deg+1, ii) = he;
            h(deg+1, ii) = 1/nodes;
            disp(sprintf('h = %g, error = %g', 1/nodes, e));
            nodes = nodes*2;
            
            for fig = 1:deg+1
                loglog(h(fig,:),err(fig,:));
                hold on;
                loglog(h(fig,:),herr(fig,:), '--');
                set(gca,'XDir','Reverse');
                grid on;
                drawnow;
            end
            hold off;
        end
    end
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

% Basis functions
function [phi, dphi] = basis(centr, point, h, dgree)
    for ii = 0:dgree
        phi(ii+1) = ((point-centr)/h)^ii;
        
        if ii == 0
            dphi(ii+1) = 0;
        else
            dphi(ii+1) = (ii/h)*((point-centr)/h)^(ii-1);
        end
    end
    phi = phi';
    dphi = dphi';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xT = cell_center(pd, ii)
    xT = (ii-0.5)*pd.h;
end

function [xF1, xF2] = face_centers(pd, ii)
    xF1 = (ii-1)*pd.h;
    xF2 = ii * pd.h;
end

function P = projector(pd, ii, fun)
    [nodes,weights] = gauss_quadrature(2*pd.K, pd.h, ii);
    
    xT = cell_center(pd, ii);
    [xF1, xF2] = face_centers(pd, ii);
    
    mass_mat = zeros(pd.K+1, pd.K+1);
    rhs = zeros(pd.K+1, 1);
    for wi = 1:length(weights)
        [phi, ~] = basis(xT, nodes(wi), pd.h, pd.K);
        mass_mat = mass_mat + weights(wi)*phi*phi';
        rhs = rhs + weights(wi)*phi*fun(nodes(wi));
    end

    pT = mass_mat \ rhs;
    pF(1) = fun(xF1);
    pF(2) = fun(xF2);
    
    P = [pT;pF'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function error = compute_error(pd, ii, fun, dofs)
    [nodes,weights] = gauss_quadrature(2*pd.K, pd.h, ii);
    
    xT = cell_center(pd, ii);
    [xF1, xF2] = face_centers(pd, ii);
    
    mass_mat = zeros(pd.K+1, pd.K+1);
    rhs = zeros(pd.K+1, 1);
    for wi = 1:length(weights)
        [phi, ~] = basis(xT, nodes(wi), pd.h, pd.K);
        mass_mat = mass_mat + weights(wi)*phi*phi';
        rhs = rhs + weights(wi)*phi*fun(nodes(wi));
    end

    pT = mass_mat \ rhs;
    
    error = dot(pT-dofs, mass_mat*(pT-dofs));
end

function cell_data = compute_cell_data(pd, ii)
    [nodes,weights] = gauss_quadrature(2*(pd.K+1), pd.h, ii);
    [MM, SM]            = cell_matrices(pd, ii, nodes, weights);
    cell_data           = struct;
    cell_data.MM        = MM;
    cell_data.SM        = SM;
    cell_data.qnodes    = nodes;
    cell_data.qweights  = weights;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MM, SM] = cell_matrices(pd, ii, nodes, weights)
    xT = cell_center(pd, ii);

    MM = zeros(pd.K+2, pd.K+2);
    SM = zeros(pd.K+2, pd.K+2);
    for wi = 1:length(weights)
        [phi, dphi] = basis(xT, nodes(wi), pd.h, pd.K+1);
        MM = MM + weights(wi)*phi*phi';
        SM = SM + weights(wi)*dphi*dphi';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = reconstruction_operator(pd, ii, stiff_mat)
    xT = cell_center(pd, ii);
    [xF1, xF2] = face_centers(pd, ii);

    MG = stiff_mat(2:end, 2:end);
    BG(:,1:pd.K+1) = stiff_mat(2:end,1:pd.K+1);
    [phiF1, dphiF1] = basis(xT, xF1, pd.h, pd.K+1);
    [phiF2, dphiF2] = basis(xT, xF2, pd.h, pd.K+1);
    
    BG(1:end, 1:pd.K+1) = BG(1:end, 1:pd.K+1) + dphiF1(2:end)*phiF1(1:pd.K+1)';
    BG(1:end, 1:pd.K+1) = BG(1:end, 1:pd.K+1) - dphiF2(2:end)*phiF2(1:pd.K+1)';
    BG(1:end, pd.K+2) = - dphiF1(2:end);
    BG(1:end, pd.K+3) = + dphiF2(2:end);
    
    R = MG\BG;
end

function tps = make_test_points(pd, ii)
    for pi = 1:pd.tp
        tps(pi) = (ii-1)*pd.h + (pi-1)*pd.h/(pd.tp-1);
    end
end

function test_gradient_reconstruction(pd)
    pos = 1;
    for ii = 1:pd.N
        xT = cell_center(pd, ii);
        data = compute_cell_data(pd, ii);
        l = projector(pd, ii, @f_solution);
        R = reconstruction_operator(pd, ii, data.SM);
        g = R*l;
        
        tps = make_test_points(pd, ii);
        for tp_i = 1:length(tps)
            [phi, dphi] = basis(xT, tps(tp_i), pd.h, pd.K+1);
            x_val(pos) = tps(tp_i);
            potential(pos) = g'*phi(2:end) + l(1);
            gradient(pos) = g'*dphi(2:end);
            pos = pos+1;
        end
    
    end
    plot(x_val, potential);
    hold on;
    plot(x_val, gradient);
end

% Compute stabilization
function S = stabilization(pd, ii, R, mass_mat)
    
    xT = cell_center(pd, ii);
    [xF1, xF2] = face_centers(pd, ii);
    
    M1 = mass_mat(1:pd.K+1,1:pd.K+1);
    M2 = mass_mat(1:pd.K+1,2:pd.K+2);
    P1 = - M1\(M2*R);
    P1(1:pd.K+1, 1:pd.K+1) = P1(1:pd.K+1, 1:pd.K+1) + eye(pd.K+1);
    
    [phiF1, ~] = basis(xT, xF1, pd.h, pd.K+1);
    MFF = 1;
    MFT = phiF1;
    P2 = MFF \ (MFT(2:end)'*R);
    P2(pd.K+2) = P2(pd.K+2)-1;
    P3 = MFF \ (MFT(1:pd.K+1)'*P1);
    B = P2 + P3;
    S = B' * MFF * B / pd.h;
    
    [phiF2, ~] = basis(xT, xF2, pd.h, pd.K+1);
    MFF = 1;
    MFT = phiF2;
    P2 = MFF \ (MFT(2:end)'*R);
    P2(pd.K+3) = P2(pd.K+3)-1;
    P3 = MFF \ (MFT(1:pd.K+1)'*P1);
    B = P2 + P3;
    S = S + B' * MFF * B / pd.h;
end

% Compute cell right hand side
function rhs = cell_rhs(pd, ii, fun)
    [nodes,weights] = gauss_quadrature(2*pd.K+2, pd.h, ii);
    
    xT = cell_center(pd, ii);
    
    rhs = zeros(pd.K+1, 1);
    for wi = 1:length(weights)
        [phi, ~] = basis(xT, nodes(wi), pd.h, pd.K);
        rhs = rhs + weights(wi)*phi*fun(nodes(wi));
    end
end

% Do static condensation
function [AC, bC] = static_condensation(pd, A, bT)
    K_TT = A(1:pd.K+1, 1:pd.K+1);
    K_TF = A(1:pd.K+1, end-1:end);
    K_FT = A(end-1:end, 1:pd.K+1);
    K_FF = A(end-1:end, end-1:end);

    AL = K_TT \ K_TF;
    bL = K_TT \ bT;
    
    AC = K_FF - K_FT * AL;
    bC = - K_FT * bL;
end

% Run HHO method
function [L2error, x_val, potential, H1error] = hho(pd, loadfun, solfun)
    cell_size = pd.K+1;
    face_offset = cell_size*pd.N;
    num_cells = pd.N;
    num_faces = pd.N+1;
    system_size = num_cells*cell_size + num_faces;
    disp(sprintf('DoFs: %d', system_size));
    
    globA = sparse(num_faces+2, num_faces+2);
    globrhs = zeros(num_faces+2);
    
    % Assemble global matrix
    for ii = 1:pd.N
        data = compute_cell_data(pd, ii);
        R = reconstruction_operator(pd, ii, data.SM);
        S = stabilization(pd, ii, R, data.MM);
        A = R' * data.SM(2:end, 2:end) * R;
        LC = A + S;
        rhs = cell_rhs(pd, ii, loadfun);
        
        [AC, bC] = static_condensation(pd, LC, rhs);
        e = ii-1;
        for dof_i = 1:2
            for dof_j = 1:2
                globA(e+dof_i, e+dof_j) = globA(e+dof_i, e+dof_j) + AC(dof_i, dof_j);
            end
            globrhs(e+dof_i) = globrhs(e+dof_i) + bC(dof_i);
        end
    end
    
    % Apply boundary conditions
    globA(num_faces+1, 1) = 1;
    globA(num_faces+2, num_faces) = 1;
    globA(1, num_faces+1) = 1;
    globA(num_faces, num_faces+2) = 1;
    
    disp(sprintf('Condition number: %g', condest(globA)));
    
    % Solve linear system
    result = globA\globrhs;
    
    L2error = 0;
    H1error = 0;
    
    % Postprocess
    pos = 1;
    for ii = 1:pd.N
        xT = cell_center(pd, ii);
        data = compute_cell_data(pd, ii);
        R = reconstruction_operator(pd, ii, data.SM);
        S = stabilization(pd, ii, R, data.MM);
        A = R' * data.SM(2:end, 2:end) * R;
        LC = A + S;
        rhs = cell_rhs(pd, ii, loadfun);
        
        solF = result(ii:ii+1)';
        
        K_TT = LC(1:pd.K+1, 1:pd.K+1);
        K_TF = LC(1:pd.K+1, end-1:end);
        
        solT = K_TT \ (rhs - K_TF*solF);
        
        tps = make_test_points(pd, ii);
        for tp_i = 1:length(tps)
            [phi, ~] = basis(xT, tps(tp_i), pd.h, pd.K);
            x_val(pos) = tps(tp_i);
            potential(pos) = solT'*phi;
            pos = pos+1;
        end
        u = projector(pd, ii, solfun);
        uh = [solT; solF];
        diffH1 = u - uh;
        diffL2 = u(1:cell_size) - solT;
        H1error = H1error + dot(diffH1', LC*diffH1);
        %error = error + compute_error(pd, ii, solfun, solT);
        L2error = L2error + dot(diffL2', data.MM(1:cell_size, 1:cell_size)*diffL2);
    end
    
    H1error = sqrt(H1error);
    L2error = sqrt(L2error);
end

function hho_eigenvalues(pd)
    cell_size = pd.K+1;
    face_offset = cell_size*pd.N;
    num_cells = pd.N;
    num_faces = pd.N+1;
    system_size = num_cells*cell_size + num_faces;
    disp(sprintf('System size: %d', system_size));
    
    globA = sparse(system_size, system_size);
    globR = sparse(system_size, system_size);
    
    for ii = 1:pd.N
        data = compute_cell_data(pd, ii);
        R = reconstruction_operator(pd, ii, data.SM);
        S = stabilization(pd, ii, R, data.MM);
        A = R' * data.SM(2:end, 2:end) * R;
        LC = A + S;
        
        sz = size(R);
        const = zeros(1, sz(2));
        const(1) = 1;
        B = [const;R]'*data.MM*[const;R];
        
        
        for dof_i = 1:cell_size
            l2g(dof_i) = (ii-1)*cell_size+dof_i;
        end
        l2g(cell_size+1) = face_offset+ii;
        l2g(cell_size+2) = face_offset+ii+1;
        
        for dof_i = 1:cell_size+2
            for dof_j = 1:cell_size+2
                globA(l2g(dof_i),l2g(dof_j)) = globA(l2g(dof_i),l2g(dof_j)) + LC(dof_i, dof_j);
                %globR(l2g(dof_i),l2g(dof_j)) = globR(l2g(dof_i),l2g(dof_j)) + B(dof_i, dof_j);
            end
        end
        
        for dof_i = 1:cell_size
            for dof_j = 1:cell_size
                globR(l2g(dof_i),l2g(dof_j)) = globR(l2g(dof_i),l2g(dof_j)) + data.MM(dof_i, dof_j);
            end
        end
        
    end
    
    % Apply boundary conditions
    globA(face_offset+1,:) = [];
    globA(:,face_offset+1) = [];
    globA(face_offset+num_faces-1,:) = [];
    globA(:,face_offset+num_faces-1) = [];
    globR(face_offset+1,:) = [];
    globR(:,face_offset+1) = [];
    globR(face_offset+num_faces-1,:) = [];
    globR(:,face_offset+num_faces-1) = [];
    
    problem_formulation = 2;
    
    if (problem_formulation == 1)
        [V,D] = eig(full(globA), full(globR));
    else
        
        NA = globA(1:face_offset, 1:face_offset);
        NB = globA(1:face_offset, face_offset+1:end);
        ND = globA(face_offset+1:end, face_offset+1:end);
        NM = globR(1:face_offset, 1:face_offset);
        NS = NA - NB*(ND\NB');
        
    
        % Solve eigenvalue problem
        [V1,D] = eig(full(NS), full(NM));
        sz = size(V1);
        for vv = 1:sz(2)
            NV = -ND\(NB'*V1(:,vv));
            V(:,vv) = [V1(:,vv);NV];
        end
    end
    DD = diag(D);
    disp(sprintf('Num of eigenvalues: %g', length(DD)));
    
    %K = [DD,V];
    [~,I] = sort(DD);
    
    lambda_h = DD(I);
    n = linspace(1,length(DD),length(DD))';
    lambda = (pi*n).^2;
    eig_error = (lambda_h-lambda)./lambda;
    
    for ev_i = 1:length(DD)
        disp(sprintf('Eigenvalue: %g', lambda_h(ev_i)));
        eigv = V(:,I(ev_i));
        [L2e, H1e, Se, De] = reconstruct_eigenvector(pd, eigv, ev_i);
        L2error(ev_i) = L2e;
        H1error(ev_i) = H1e;
        Serror(ev_i) = Se;
        Derror(ev_i) = De;
    end
    
    plot(eig_error);
    hold on;
    plot(L2error);
    He = H1error'./lambda;
    plot(He);
    Se = Serror'./lambda;
    plot(Se);
    De = Derror'./lambda;
    plot(Derror);
    plot(eig_error + L2error' + He - Se - De)
    
    
   
    legend('Eig', 'L2', 'H1', 'S', 'D', 'Sum');
end

function [L2error, H1error, Serror, Derror] = reconstruct_eigenvector(pd, ev, ev_i)
    cell_size = pd.K+1;
    face_offset = cell_size*pd.N;
    Lpos = 1;
    Hpos = 1;
    H1error = 0;
    Serror = 0;
    Derror = 0;
    
    factor = 1;
%     if (pd.K == 0)
%         factor = sqrt(2);
%     else
%         factor = 1;
%     end
    
    for ii = 1:pd.N
        xT = cell_center(pd, ii);
        [xF1, xF2] = face_centers(pd, ii);
        
        data = compute_cell_data(pd, ii);
        R = reconstruction_operator(pd, ii, data.SM);
        S = stabilization(pd, ii, R, data.MM);
        A = R' * data.SM(2:end, 2:end) * R;
        
        start = (ii-1) * (pd.K+1) + 1;
        stop = ii * (pd.K+1);
        
        [nodes,weights] = gauss_quadrature(2*pd.K, pd.h, ii);
        for kk = 1:length(nodes)
            [phi, ~] = basis(xT, nodes(kk), pd.h, pd.K);
            x_val(Lpos) = nodes(kk);
            yh = dot( ev(start:stop), phi );
            y = factor*sin(ev_i*pi*nodes(kk));
            w(Lpos) = weights(kk);
            yh_val(Lpos) = yh;
            y_val(Lpos) = y;
            Lpos = Lpos+1;
        end
        
        if (ii == 1)
            uh = [ev(start:stop); 0; ev(face_offset+ii)];
        else if (ii == pd.N)
            uh = [ev(start:stop); ev(face_offset+ii-1); 0];
            else
                uh = [ev(start:stop); ev(face_offset+ii-1); ev(face_offset+ii)];
            end
        end
        
        Serror = Serror + uh'*S*uh;
        
        g = R*uh;
        [nodes,weights] = gauss_quadrature(2*pd.K+2, pd.h, ii);
        for kk = 1:length(nodes)
            [~, dphi] = basis(xT, nodes(kk), pd.h, pd.K+1);
            dx_val(Hpos) = nodes(kk);
            
            dyh = dot( g, dphi(2:end) );                    % approx gradient
            dy = factor*ev_i*pi*cos(ev_i*pi*nodes(kk));    % expected gradient
            dw(Hpos) = weights(kk);                         % weights
            
            dyh_val(Hpos) = dyh;
            dy_val(Hpos) = dy;
            Hpos = Hpos+1;
        end
        
        phiF1 = basis(xT, xF1, pd.h, pd.K); 
        
        Derror = Derror - uh(cell_size+1);
        Derror = Derror + uh(cell_size+2);
    end
    
    if (yh_val(2) < yh_val(1))
        yh_val = -yh_val;
    end
    
    if (dyh_val(6) > dyh_val(1))
        dyh_val = -dyh_val;
    end
    
    L2error = dot(w, (yh_val - y_val).^2);
    H1error = dot(dw, (dyh_val - dy_val).^2);
    
%       plot(x_val, yh_val);
%       hold on;
%       plot(x_val, y_val);
%       plot(dx_val, dyh_val);
%       plot(dx_val, dy_val);
%       hold off;
%       pause;
end


