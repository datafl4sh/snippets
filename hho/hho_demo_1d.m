%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       /\
%      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
%     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
%    /\    /\
%   /__\  /__\    1D Hybrid High-order demo code
%  /_\/_\/_\/_\
%

function hho()

    what = 4;

    if (what == 1)
        plot_hho_convergence();
    end
    
    N = 16;
    K = 1;
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

function P = projector_k_plus_one(pd, ii, fun)
    [nodes,weights] = gauss_quadrature(2*pd.K+2, pd.h, ii);
    
    xT = cell_center(pd, ii);
    
    mass_mat = zeros(pd.K+2, pd.K+2);
    rhs = zeros(pd.K+2, 1);
    for wi = 1:length(weights)
        [phi, ~] = basis(xT, nodes(wi), pd.h, pd.K+1);
        mass_mat = mass_mat + weights(wi)*phi*phi';
        rhs = rhs + weights(wi)*phi*fun(nodes(wi));
    end

    P = mass_mat \ rhs;
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
    cell_data.R         = reconstruction_operator(pd, ii, SM);
    cell_data.S         = stabilization(pd, ii, cell_data.R, MM);
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
        R = data.R;
        S = data.S;
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

function [x,y] = hho_eval(pd, ii, poly, N)
    xT = cell_center(pd, ii);
    [xF1, xF2] = face_centers(pd, ii);
    x = linspace(xF1, xF2, N);
    
    for kk = 1:length(x)
        [phi, ~] = basis(xT, x(kk), pd.h, pd.K+1);
        y(kk) = dot(poly, phi);
    end
end

function [x,y] = hho_eval_interior(pd, ii, poly, N)
    xT = cell_center(pd, ii);
    [xF1, xF2] = face_centers(pd, ii);
    x = linspace(xF1, xF2, N);
    
    for kk = 1:length(x)
        [phi, ~] = basis(xT, x(kk), pd.h, pd.K);
        y(kk) = dot(poly, phi);
    end
end



function [uh, Th] = element_dofs(pd, ii, solution)
    face_offset = (pd.K+1)*pd.N;
    start = (ii-1) * (pd.K+1) + 1;
    stop = ii * (pd.K+1);
    uh = solution(start:stop);
    
    if (ii == 1)
        Th = [0; solution(face_offset+ii)];
    else if (ii == pd.N)
        Th = [solution(face_offset+ii-1); 0];
        else
            Th = [solution(face_offset+ii-1); solution(face_offset+ii)];
        end
    end
end

function hho_plot_eigfun(pd, solution, ev_i)
    xvals = [];
    yvals = [];
    yivals = [];
    for ii = 1:pd.N
        [uh, Th] = element_dofs(pd, ii, solution);
        data = compute_cell_data(pd, ii);
        R = reconstruction_operator(pd, ii, data.SM);
        g = R*[uh;Th];
        g = [uh(1);g];
        
        eef = @(x) sqrt(2)*sin(ev_i*pi*x);
        pg = projector_k_plus_one(pd, ii, eef);
        
        if (dot(pg, data.MM*g) < 0)
            g = -g;
            uh = -uh;
        end
        
        [x,y] = hho_eval(pd, ii, g, 20);
        [x,yi] = hho_eval_interior(pd, ii, uh, 20);
        xvals = [xvals,x];
        yvals = [yvals,y];
        yivals = [yivals,yi];
    end
    
    subplot(2,1,1)
    plot(xvals, sqrt(2)*sin(ev_i*pi*xvals));
    hold on;
    plot(xvals, yvals);
    plot(xvals, yivals);
    hold off;
    subplot(2,1,2)
    plot(xvals, yvals-yivals);
    
    
    %w = hamming(length(yvals));
    %w = w/mean(w);
    %fftdata = abs(fft(yvals' .* w));
    %plot( fftdata(1:length(fftdata)/2+1) )
end

function [L2, H1, D] = compute_L2_H1_D_error(pd, ii, uhTh, ev_i, alpha)
    xT = cell_center(pd, ii);
    data = compute_cell_data(pd, ii);
    R = reconstruction_operator(pd, ii, data.SM);
    g = R*[uhTh];
    g = [uhTh(1);g];

    eef = @(x) sqrt(2)*sin(ev_i*pi*x);
    pg = projector_k_plus_one(pd, ii, eef);

    if (dot(pg, data.MM*g) < 0)
        g = -g;
    end

    L2 = 0;
    H1 = 0;
    [nodes,weights] = gauss_quadrature(2*pd.K+2, pd.h, ii);
    for kk = 1:length(nodes)
        [phi, dphi] = basis(xT, nodes(kk), pd.h, pd.K+1);
        y = sqrt(2)*sin(ev_i*pi*nodes(kk));
        yh = dot(g, phi);
        dy = sqrt(2)*ev_i*pi*cos(ev_i*pi*nodes(kk));
        dyh = dot(g(2:end), dphi(2:end));
        L2 = L2 + alpha * weights(kk) * (y-yh/alpha) * (y-yh/alpha);
        H1 = H1 + weights(kk) * (dy-dyh) * (dy-dyh);
    end
    
    xT = cell_center(pd, ii);
    [xF1, xF2] = face_centers(pd, ii);
    [phi, ~] = basis(xT, xF1, pd.h, pd.K+1);
    D = - sqrt(2)*ev_i*pi*cos(ev_i*pi*xF1)*dot(phi,g);
    [phi, ~] = basis(xT, xF2, pd.h, pd.K+1);
    D = D + sqrt(2)*ev_i*pi*cos(ev_i*pi*xF2)*dot(phi,g);
end

function h1 = h1_error(pd, ii, ev_i, cell_data, uhTh)
    h1 = 0;    
    [nodes,weights] = gauss_quadrature(2*pd.K+2, pd.h, ii);
    
    %compute reconstruction
    g = cell_data.R*uhTh;         
    gg = [uhTh(1);g];
    
    %determine correct sign
    eef = @(x) sqrt(2)*sin(ev_i*pi*x);
    pg = projector_k_plus_one(pd, ii, eef);

    if (dot(pg, cell_data.MM*gg) < 0)
        g = -g;
    end
    
    xT = cell_center(pd, ii);
    
    for kk = 1:length(nodes)
        [~, dphi] = basis(xT, nodes(kk), pd.h, pd.K+1);
        y = sqrt(2)*ev_i*pi*cos(ev_i*pi*nodes(kk));
        yh = dot(g, dphi(2:end));
        h1 = h1 + weights(kk) * (y-yh) * (y-yh);
    end
end

function l2 = l2_on_interior(pd, ii, ev_i, cell_data, uh)
    l2 = 0;    
    [nodes,weights] = gauss_quadrature(2*pd.K, pd.h, ii);
    
    %determine correct sign
    eef = @(x) sqrt(2)*sin(ev_i*pi*x);
    puh = projector(pd, ii, eef);

    if (dot(puh(1:pd.K+1), cell_data.MM(1:pd.K+1, 1:pd.K+1)*uh) < 0)
        uh = -uh;
    end
    
    xT = cell_center(pd, ii);
    
    for kk = 1:length(nodes)
        [phi, ~] = basis(xT, nodes(kk), pd.h, pd.K);
        y = sqrt(2)*sin(ev_i*pi*nodes(kk));
        yh = dot(uh, phi);
        l2 = l2 + weights(kk) * (y-yh) * (y-yh);
    end
end


function l2 = l2_on_reconstruction(pd, ii, ev_i, cell_data, uhTh)
    l2 = 0;    
    [nodes,weights] = gauss_quadrature(2*pd.K+2, pd.h, ii);
    
    %compute reconstruction
    g = cell_data.R*uhTh;         
    g = [uhTh(1);g];
    
    %determine correct sign
    eef = @(x) sqrt(2)*sin(ev_i*pi*x);
    pg = projector_k_plus_one(pd, ii, eef);

    if (dot(pg, cell_data.MM*g) < 0)
        g = -g;
    end
    
    xT = cell_center(pd, ii);
    
    for kk = 1:length(nodes)
        [phi, ~] = basis(xT, nodes(kk), pd.h, pd.K+1);
        y = sqrt(2)*sin(ev_i*pi*nodes(kk));
        yh = dot(g, phi);
        l2 = l2 + weights(kk) * (y-yh) * (y-yh);
    end
end

function serr = s_error(pd, ii, cell_data, uhTh)
    serr = dot(uhTh, cell_data.S*uhTh);
end

function derr = jump_error(pd, ii, ev_i, cell_data, uhTh)
    xT = cell_center(pd, ii);
    g = cell_data.R*[uhTh];
    g = [uhTh(1);g];

    eef = @(x) sqrt(2)*sin(ev_i*pi*x);
    pg = projector_k_plus_one(pd, ii, eef);

    if (dot(pg, cell_data.MM*g) < 0)
        g = -g;
    end

    [xF1, xF2] = face_centers(pd, ii);
    [phi, ~] = basis(xT, xF1, pd.h, pd.K+1);
    derr = - sqrt(2)*ev_i*pi*cos(ev_i*pi*xF1)*dot(phi,g);
    [phi, ~] = basis(xT, xF2, pd.h, pd.K+1);
    derr = derr + sqrt(2)*ev_i*pi*cos(ev_i*pi*xF2)*dot(phi,g);
end

function Rn = R_norm(pd, ii, cell_data, uhTh)
    g = cell_data.R*uhTh;
    g = [uhTh(1);g];
    Rn = dot(g, cell_data.MM*g);   
end

function [L2error, H1error, Serror, Rn, Derror] = compute_errors(pd, ev, ev_i, alpha)

    L2error = 0;
    H1error = 0;
    Serror = 0;
    Rn = 0;
    Derror = 0;
    for ii = 1:pd.N
        [uh, Th] = element_dofs(pd, ii, ev);
        [L2, H1, D] = compute_L2_H1_D_error(pd, ii, [uh;Th], ev_i, alpha);
        L2error = L2error + L2;
        H1error = H1error + H1;
        Derror = Derror + D;
        
        data = compute_cell_data(pd, ii);
        R = reconstruction_operator(pd, ii, data.SM);
        S = stabilization(pd, ii, R, data.MM);
        Serror = Serror + dot([uh;Th]', S*[uh;Th]);
        
        g = R*[uh;Th];
        g = [uh(1);g];
        Rn = Rn + dot(g, data.MM*g)/alpha;
        
    end
end

function err = compute_errs(pd, ev, ev_i)
    err = struct;
    
    L2R     = 0;
    L2I     = 0;
    H1      = 0;
    S       = 0;
    D       = 0;
    Rn      = 0;
    
    for ii = 1:pd.N
        data        = compute_cell_data(pd, ii);
        [uh, Th]    = element_dofs(pd, ii, ev);
        uhTh        = [uh;Th];
        
        L2R     = L2R + l2_on_reconstruction(pd, ii, ev_i, data, uhTh);
        L2I     = L2I + l2_on_interior(pd, ii, ev_i, data, uh);
        H1      = H1 + h1_error(pd, ii, ev_i, data, uhTh);
        S       = S + s_error(pd, ii, data, uhTh);
        D       = D + jump_error(pd, ii, ev_i, data, uhTh);
        Rn      = Rn + R_norm(pd, ii, data, uhTh);
    end
    
    err.L2R = L2R;
    err.L2I = L2I;
    err.H1  = H1;
    err.S   = S;
    err.D   = D;
    err.Rn  = Rn;
end

function new_eigv = normalize_eigfun(pd, eigv)
    new_eigv = zeros(size(eigv));
    for ii = 1:pd.N
%         [uh, Th, cb, ce, f1o, f2o] = element_dofs(pd, ii, eigv);
%         data = compute_cell_data(pd, ii);
%         R = reconstruction_operator(pd, ii, data.MM);
%         g = R*[uh;Th];
%         g = [uh(1);g];
%         nf = dot(g, data.SM*g);
%         uh = uh/sqrt(nf);
%         Th = Th/sqrt(nf);
%         new_eigv(cb:ce) = uh;
%         new_eigv(f1o) = Th(1);
%         new_eigv(f2o) = Th(2);

         [uh, Th, cb, ce, f1o, f2o] = element_dofs(pd, ii, eigv);
         data = compute_cell_data(pd, ii);
         %nf = dot(uh, data.MM(1:pd.K+1, 1:pd.K+1)*uh);
         %uh = uh;%/sqrt(nf);
         %Th = Th;%/sqrt(nf);
         %new_eigv(cb:ce) = uh;
         %new_eigv(f1o) = Th(1);
         %new_eigv(f2o) = Th(2);


    end
end

function hho_eigenvalues(pd)
    cell_size = pd.K+1;
    face_offset = cell_size*pd.N;
    num_cells = pd.N;
    num_faces = pd.N+1;
    system_size = num_cells*cell_size + num_faces;
    disp(sprintf('System size: %d', system_size));
    
    globK = sparse(system_size, system_size);
    globM = sparse(system_size, system_size);
    %globM_f = sparse(system_size, system_size);
    
    for ii = 1:pd.N
        data = compute_cell_data(pd, ii);
        R = reconstruction_operator(pd, ii, data.SM);
        S = stabilization(pd, ii, R, data.MM);
        K = R' * data.SM(2:end, 2:end) * R;
        LC = K + S;
        
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
                globK(l2g(dof_i),l2g(dof_j)) = globK(l2g(dof_i),l2g(dof_j)) + LC(dof_i, dof_j);
                %globM_f(l2g(dof_i),l2g(dof_j)) = globM_f(l2g(dof_i),l2g(dof_j)) + B(dof_i, dof_j);
            end
        end
        
        %for dof_i = 1:cell_size
        %    for dof_j = 1:cell_size
        %        globM_f(l2g(dof_i),l2g(dof_j)) = globM_f(l2g(dof_i),l2g(dof_j)) + data.MM(dof_i, dof_j);
        %    end
        %end
        %globM_f(l2g(cell_size+1),l2g(cell_size+1)) =  0;
        %globM_f(l2g(cell_size+2),l2g(cell_size+2)) =  0;
        
        for dof_i = 1:cell_size
            for dof_j = 1:cell_size
                globM(l2g(dof_i),l2g(dof_j)) = globM(l2g(dof_i),l2g(dof_j)) + data.MM(dof_i, dof_j);
            end
        end
        
    end
    
    % Apply boundary conditions
    globK(face_offset+1,:) = [];
    globK(:,face_offset+1) = [];
    globK(face_offset+num_faces-1,:) = [];
    globK(:,face_offset+num_faces-1) = [];
    globM(face_offset+1,:) = [];
    globM(:,face_offset+1) = [];
    globM(face_offset+num_faces-1,:) = [];
    globM(:,face_offset+num_faces-1) = [];
    %globM_f(face_offset+1,:) = [];
    %globM_f(:,face_offset+1) = [];
    %globM_f(face_offset+num_faces-1,:) = [];
    %globM_f(:,face_offset+num_faces-1) = [];
    
    
        
    NA = globK(1:face_offset, 1:face_offset);
    NB = globK(1:face_offset, face_offset+1:end);
    ND = globK(face_offset+1:end, face_offset+1:end);
    NM = globM(1:face_offset, 1:face_offset);
    NS = NA - NB*(ND\NB');


    % Solve eigenvalue problem
    [V1,D] = eig(full(NS), full(NM));
    sz = size(V1);
    for vv = 1:sz(2)
        NV = -ND\(NB'*V1(:,vv));
        ev = V1(:,vv);
        V(:,vv) = [ev;NV];
    end

    DD = diag(D);
    disp(sprintf('Num of eigenvalues: %g', length(DD)));
    
    %K = [DD,V];
    [~,I] = sort(DD);
    
    lambda_h = DD(I);
    n = linspace(1,length(DD),length(DD))';
    lambda = (pi*n).^2;
    
    lambda_err = (lambda_h - lambda)./lambda;
    
    for ev_i = 1:length(DD)
        disp(sprintf('Eigenvalue: %g', lambda_h(ev_i)));
        eigv = V(:,I(ev_i));                %get eigvector
        nf = dot(eigv,globM*eigv);          %compute eigenfunction norm
        eigv = eigv/sqrt(nf);               %normalize eigenfunction
        %hho_plot_eigfun(pd, eigv, ev_i);
        %pause;
        errs = compute_errs(pd, eigv, ev_i);
        
        L2R(ev_i)   = errs.L2R;
        L2I(ev_i)   = errs.L2I;
        H1(ev_i)    = errs.H1;
        Se(ev_i)    = errs.S;
        De(ev_i)    = errs.D;
        Rn(ev_i)    = errs.Rn;
    end
    
    H1 = H1' ./ lambda;
    Se = Se' ./ lambda;
    De = De' ./ lambda;
    
    subplot(2,1,1);
    plot(lambda_err);
    hold on;
    plot(L2R);
    %plot(L2I);
    plot(Rn);
    plot(H1);
    plot(Se);
    plot(De);
    
    Sum = lambda_err + L2R' + 1 - Rn' - Se - 2*De;
    plot(Sum);
    grid on;
    legend('Eig', 'L2R', 'Rn', 'H1', 'S', 'Jump', 'Sum');
    subplot(2,1,2);
    semilogy(abs(H1-Sum));
    grid on;
    legend('H1 - Sum');
    
    
    %He = H1error'./lambda;                  % H1 norm error
    %lambdae = (lambda_h-lambda)./lambda;    % Eigenvalue error
    %lambdae = (lambda.*(1-2*alpha') + lambda_h)./lambda;
    %Le = L2error';                          % L2 norm error
    %Se = Serror'./(lambda .* (alpha.^2)');                   % Stabilization error
    %De = Derror'./(lambda .* alpha'); 
    % Jump error

    %subplot(2,1,1);
    %plot(He);
    %hold on;
    %plot(Le);
    %plot(lambdae);
    %plot(Se);
    %plot(Rns);
    %plot(De);
    %Sum = lambdae + Le + (1 - Rns').*alpha' - Se - De;
    %plot(Sum);
    %grid on;
    %legend('H1', 'L2', 'Eig', 'S', 'Rn', 'Delta', 'Sum');
    %subplot(2,1,2);
    %semilogy(abs(He-Sum));
    %grid on;
    %legend('H1 - Sum');
end

