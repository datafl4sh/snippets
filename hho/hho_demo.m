function  hho_demo()
    
     for deg = 0:4
         maxii = 7;
         nodes = 4;
         for ii = 1:maxii
             disp(sprintf('**** N = %d, K = %d ****', nodes, deg));
             pd = initialize_hho(nodes, deg);
             e = hho(pd, @f_load, @f_solution);
             err(deg+1, ii) = e;
             h(deg+1, ii) = 1/nodes;
             disp(sprintf('h = %g, error = %g', 1/nodes, e));
             nodes = nodes*2;
             
             for fig = 1:deg+1
                 loglog(h(fig,:),err(fig,:));
                 set(gca,'XDir','Reverse');
                 hold on;
                 grid on;
                 drawnow;
             end
             hold off;
         end
     end

##nodes = 4;
##deg = 0;
##pd = initialize_hho(nodes, deg);
##show_gr(pd);
end
    

function pd = initialize_hho(N, K)
    pd = struct;
    pd.N = N;
    pd.K = K;
    pd.h = 1/pd.N;
    pd.tp = 8;
end

function l = f_load(x)
    l = sin(pi*x) * pi^2;
    return;
end

function l = f_solution(x)
    l = sin(pi*x);
    return;
end

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

function [phi, dphi] = basis(centr, point, h, dgree)
    h_scale = 2;
    for ii = 0:dgree
        phi(ii+1) = ((point-centr)/(h/h_scale))^ii;
        
        if ii == 0
            dphi(ii+1) = 0;
        else
            dphi(ii+1) = (h_scale*ii/h)*((point-centr)/(h/h_scale))^(ii-1);
        end
    end
    phi = phi';
    dphi = dphi';
end

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

function [GT, A] = gradient_reconstruction(pd, ii)

    [nodes,weights] = gauss_quadrature(2*(pd.K+1), pd.h, ii);
    
    xT = cell_center(pd, ii);
    [xF1, xF2] = face_centers(pd, ii);

    stiff_mat = zeros(pd.K+2, pd.K+2);
    MG = zeros(pd.K+1, pd.K+1);
    BG = zeros(pd.K+1, pd.K+3);
    
    for wi = 1:length(weights)
        [~, dphi] = basis(xT, nodes(wi), pd.h, pd.K+1);
        stiff_mat = stiff_mat + weights(wi)*dphi*dphi';
    end
    ##disp(stiff_mat);
    MG = stiff_mat(2:end, 2:end);
    BG(:,1:pd.K+1) = stiff_mat(2:end,1:pd.K+1);
    [phiF1, dphiF1] = basis(xT, xF1, pd.h, pd.K+1);
    [phiF2, dphiF2] = basis(xT, xF2, pd.h, pd.K+1);
    
    BG(1:end, 1:pd.K+1) = BG(1:end, 1:pd.K+1) + dphiF1(2:end)*phiF1(1:pd.K+1)';
    BG(1:end, 1:pd.K+1) = BG(1:end, 1:pd.K+1) - dphiF2(2:end)*phiF2(1:pd.K+1)';
    BG(1:end, pd.K+2) = - dphiF1(2:end);
    BG(1:end, pd.K+3) = + dphiF2(2:end);
    
    ##disp(BG);
    
    GT = MG\BG;
    A = BG'*GT;
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
        l = projector(pd, ii, @f_solution);
        [GT, ~, ~] = gradient_reconstruction(pd, ii);
        g = GT*l;
        
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

function show_gr(pd)
    hold on;
    pos = 1;
    for ii = 1:pd.N
        xT = cell_center(pd, ii);
        l = projector(pd, ii, @f_solution);
        [GT, ~] = gradient_reconstruction(pd, ii);
        g = GT*l;
        
        tps = make_test_points(pd, ii);
        for tp_i = 1:length(tps)
            [phi, ~] = basis(xT, tps(tp_i), pd.h, pd.K+1);
            x_val(pos) = tps(tp_i);
            potential(pos) = g'*phi(2:end) + l(1);
            proj(pos) = dot( l(1:end-2), phi(1:end-1));
            pos = pos+1;
        end
        [xF1, xF2] = face_centers(pd, ii);
        plot(xF1, l(end-1), '*k');
        plot(xF2, l(end), '*k');
    end
    plot(x_val, potential);
    plot(x_val, proj);
end

function S = stabilization(pd, ii, GT)
    [nodes,weights] = gauss_quadrature(2*pd.K+2, pd.h, ii);
    
    xT = cell_center(pd, ii);
    [xF1, xF2] = face_centers(pd, ii);
    
    mass_mat = zeros(pd.K+2, pd.K+2);
    for wi = 1:length(weights)
        [phi, ~] = basis(xT, nodes(wi), pd.h, pd.K+1);
        mass_mat = mass_mat + weights(wi)*phi*phi';
    end
    
    M1 = mass_mat(1:pd.K+1,1:pd.K+1);
    M2 = mass_mat(1:pd.K+1,2:pd.K+2);
    P1 = - M1\(M2*GT);
    P1(1:pd.K+1, 1:pd.K+1) = P1(1:pd.K+1, 1:pd.K+1) + eye(pd.K+1);
    
    [phiF1, ~] = basis(xT, xF1, pd.h, pd.K+1);
    MFF = 1;
    MFT = phiF1;
    P2 = MFF \ (MFT(2:end)'*GT);
    P2(pd.K+2) = P2(pd.K+2)-1;
    P3 = MFF \ (MFT(1:pd.K+1)'*P1);
    B = P2 + P3;
    S = B' * MFF * B / pd.h;
    
    [phiF2, ~] = basis(xT, xF2, pd.h, pd.K+1);
    MFF = 1;
    MFT = phiF2;
    P2 = MFF \ (MFT(2:end)'*GT);
    P2(pd.K+3) = P2(pd.K+3)-1;
    P3 = MFF \ (MFT(1:pd.K+1)'*P1);
    B = P2 + P3;
    S = S + B' * MFF * B / pd.h;
end

function rhs = cell_rhs(pd, ii, fun)
    [nodes,weights] = gauss_quadrature(2*pd.K+2, pd.h, ii);
    
    xT = cell_center(pd, ii);
    
    rhs = zeros(pd.K+1, 1);
    for wi = 1:length(weights)
        [phi, ~] = basis(xT, nodes(wi), pd.h, pd.K);
        rhs = rhs + weights(wi)*phi*fun(nodes(wi));
    end
end

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


function globerr = hho(pd, loadfun, solfun)
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
        [GT, A] = gradient_reconstruction(pd, ii);
        S = stabilization(pd, ii, GT);
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
    
    % Solve linear system
    result = globA\globrhs;
    
    % Postprocess
    error = 0;
    pos = 1;
    for ii = 1:pd.N
        xT = cell_center(pd, ii);
        [GT, A] = gradient_reconstruction(pd, ii);
        S = stabilization(pd, ii, GT);
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
        
        error = error + compute_error(pd, ii, solfun, solT);
    end
    
    globerr = sqrt(error);
end
