#include <iostream>
#include <vector>
#include <stdexcept>
#include <cstdint>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cassert>

#include <silo.h>

template<typename T>
class bounding_box
{
    T m_min_x, m_min_y, m_max_x, m_max_y;

public:
    bounding_box()
        : m_min_x(0), m_min_y(0), m_max_x(1), m_max_y(1)
    {}

    bounding_box(T p_min_x, T p_min_y, T p_max_x, T p_max_y)
        : m_min_x(p_min_x), m_min_y(p_min_y), m_max_x(p_max_x), m_max_y(p_max_y)
    {}

    T min_x() const { return m_min_x; }
    T max_x() const { return m_max_x; }
    T min_y() const { return m_min_y; }
    T max_y() const { return m_max_y; }

    T x_size() const { return m_max_x - m_min_x; }
    T y_size() const { return m_max_y - m_min_y; }
};

template<typename CoordT, typename ValueT = CoordT>
class fd_grid
{
    std::vector<ValueT>     m_u;
    bounding_box<CoordT>    m_bb;
    int32_t                 m_x_size, m_y_size;

public:
    using coordinate_type = CoordT;
    using value_type =      ValueT;
    
    fd_grid()
        : m_x_size(0), m_y_size(0), m_u(9)
    {}

    fd_grid(int32_t x_size, int32_t y_size)
        : m_x_size(x_size), m_y_size(y_size)
    {
        auto x_points = m_x_size+3;
        auto y_points = m_y_size+3;

        m_u.resize(x_points*y_points);
        for (auto& v : m_u)
            v = 0.0;
    }

    fd_grid(int32_t x_size, int32_t y_size, const bounding_box<CoordT>& p_bb)
        : m_x_size(x_size), m_y_size(y_size)
    {
        auto x_points = m_x_size+3;
        auto y_points = m_y_size+3;

        m_u.resize(x_points*y_points);
        for (auto& v : m_u)
            v = 0.0;

        m_bb = p_bb;
    }

    fd_grid(fd_grid&& other)
    {
        m_x_size = other.m_x_size;
        m_y_size = other.m_y_size;
        m_u = std::move(other.m_u);
        m_bb = other.m_bb;
    }

    fd_grid& operator=(fd_grid&& other)
    {
        m_x_size = other.m_x_size;
        m_y_size = other.m_y_size;
        m_u = std::move(other.m_u);
        m_bb = other.m_bb;
        return *this;
    }

    ValueT operator()(int32_t i, int32_t j) const
    {
        if (i < -1 or j < -1 or i > m_x_size+1 or j > m_y_size+1)
            throw std::out_of_range("access out of range");

        i = i+1;
        j = j+1;
        auto y_points = m_y_size+3;
    
        auto ofs = i*y_points + j;
        assert(ofs < m_u.size());
        return m_u[ofs];
    }

    ValueT& operator()(int32_t i, int32_t j)
    {
        if (i < -1 or j < -1 or i > m_x_size+1 or j > m_y_size+1)
            throw std::out_of_range("access out of range");

        i = i+1;
        j = j+1;
        auto y_points = m_y_size+3;

        auto ofs = i*y_points + j;
        assert(ofs < m_u.size());
        return m_u[ofs];
    }

    int32_t x_size() const { return m_x_size; }
    int32_t y_size() const { return m_y_size; }

    CoordT dx() const { return m_bb.x_size()/m_x_size; }
    CoordT dy() const { return m_bb.y_size()/m_y_size; }

    CoordT grid_to_domain_x(int32_t nx) const {
        if (nx < 0 or nx > m_x_size+1)
            throw std::invalid_argument("node not in grid");

        return m_bb.min_x() + nx*dx();
    }

    CoordT grid_to_domain_y(int32_t ny) const {
        if (ny < 0 or ny > m_y_size+1)
            throw std::invalid_argument("node not in grid");

        return m_bb.min_y() + ny*dy();
    }

    bool is_compatible_with(const fd_grid& other) const
    {
        return (m_x_size == other.m_x_size) and (m_y_size == other.m_y_size);
    }
};

enum class bc_type {
    DIRICHLET,
    NEUMANN,
//    ROBIN
};

enum class boundary {
    TOP,
    BOTTOM,
    LEFT,
    RIGHT
};

template<typename FDGridT, typename Function>
void apply_bc_top(FDGridT& grid, boundary b, bc_type bc, const Function& f)
{
    auto x_size = grid.x_size();
    auto y_size = grid.y_size();
    for (size_t i = 0; i < x_size; i++)
    {
        auto x = grid.grid_to_domain_x(i);
        auto y = grid.grid_to_domain_y(y_size);
        auto val = f(x,y);
     
        switch (bc)
        {
            case bc_type::DIRICHLET:
                grid(i, y_size) = val;
                grid(i, y_size+1) = val;
                break;

            case bc_type::NEUMANN:
                grid(i, y_size+1) = grid(i, y_size);
                break;
        }
    }
}

template<typename FDGridT, typename Function>
void apply_bc_bottom(FDGridT& grid, boundary b, bc_type bc, const Function& f)
{
    auto x_size = grid.x_size();
    for (size_t i = 0; i < x_size; i++)
    {
        auto x = grid.grid_to_domain_x(i);
        auto y = grid.grid_to_domain_y(0);
        auto val = f(x,y);
        
        switch (bc)
        {
            case bc_type::DIRICHLET:
                grid(i, 0) = val;
                grid(i, -1) = val;
                break;

            case bc_type::NEUMANN:
                grid(i, -1) = grid(i, 0);
                break;
        }
    }
}

template<typename FDGridT, typename Function>
void apply_bc_left(FDGridT& grid, boundary b, bc_type bc, const Function& f)
{
    auto y_size = grid.y_size();
    for (size_t j = 0; j < y_size; j++)
    {
        auto x = grid.grid_to_domain_x(0);
        auto y = grid.grid_to_domain_y(j);
        auto val = f(x,y);
        
        switch (bc)
        {
            case bc_type::DIRICHLET:
                grid(0, j) = val;
                grid(-1, j) = val;
                break;

            case bc_type::NEUMANN:
                grid(-1, j) = grid(0, j);
                break;
        }
    }
}

template<typename FDGridT, typename Function>
void apply_bc_right(FDGridT& grid, boundary b, bc_type bc, const Function& f)
{
    auto x_size = grid.x_size();
    auto y_size = grid.y_size();
    for (size_t j = 0; j < y_size; j++)
    {
        auto x = grid.grid_to_domain_x(x_size);
        auto y = grid.grid_to_domain_y(j);
        auto val = f(x,y);
        
        switch (bc)
        {
            case bc_type::DIRICHLET:
                grid(x_size, j) = val;
                grid(x_size+1, j) = val;
                break;

            case bc_type::NEUMANN:
                grid(x_size+1, j) = grid(x_size, j);
                break;
        }
    }
}

template<typename FDGridT, typename Function>
void apply_bc(FDGridT& grid, boundary b, bc_type bc, const Function& f)
{
    switch(b)
    {
        case boundary::TOP:
            apply_bc_top(grid, b, bc, f);
            break;

        case boundary::BOTTOM:
            apply_bc_bottom(grid, b, bc, f);
            break;

        case boundary::LEFT:
            apply_bc_left(grid, b, bc, f);
            break;

        case boundary::RIGHT:
            apply_bc_right(grid, b, bc, f);
            break;
    }
}

template<typename FDGridT>
void apply_dirichlet(FDGridT& grid, boundary b, const typename FDGridT::value_type& val)
{
    using coordinate_type = typename FDGridT::coordinate_type;

    auto fun = [&](const coordinate_type& x, const coordinate_type& y) -> auto {
        return val;
    };

    apply_bc(grid, b, bc_type::DIRICHLET, fun);
}

template<typename FDGridT>
void apply_neumann(FDGridT& grid, boundary b, const typename FDGridT::value_type& val)
{
    using coordinate_type = typename FDGridT::coordinate_type;

    auto fun = [&](const coordinate_type& x, const coordinate_type& y) -> auto {
        return val;
    };

    apply_bc(grid, b, bc_type::NEUMANN, fun);
}

template<typename T>
struct vect
{
    T x, y;

    vect() : x(0), y(0) {}
    vect(T v) : x(v), y(v) {}
    vect(T px, T py) : x(px), y(py) {}
};

template<typename T>
using scalar_fd_grid = fd_grid<T,T>;

template<typename T>
using vector_fd_grid = fd_grid<T, vect<T>>;

template<typename T>
int
visit_dump(const fd_grid<T>& g, int32_t t)
{
    std::stringstream ss;
    ss << "timestep_" << t << ".silo";
    
    DBfile *db = nullptr;
    db = DBCreate(ss.str().c_str(), DB_CLOBBER, DB_LOCAL, "Wave equation", DB_HDF5);
    if (!db)
    {
        std::cout << "Cannot write simulation output" << std::endl;
        return -1;
    }
    
    std::vector<double> x(g.x_size()+1);
    std::vector<double> y(g.y_size()+1);
    
    for (size_t i = 0; i < x.size(); i++)
        x.at(i) = double(i)/g.x_size();
        
    for (size_t i = 0; i < y.size(); i++)
        y.at(i) = double(i)/g.y_size();
        
    int dims[] = { int(x.size()), int(y.size())};
    int ndims = 2;
    double *coords[] = {x.data(), y.data()};
    
    DBPutQuadmesh(db, "mesh", NULL, coords, dims, ndims, DB_DOUBLE, DB_COLLINEAR, NULL);
    
    std::vector<double> data(x.size() * y.size());
    
    for (size_t i = 0; i < x.size(); i++)
        for (size_t j = 0; j < y.size(); j++)
            data.at(j*x.size()+i) = g(i,j);
        
    DBPutQuadvar1(db, "wave", "mesh", data.data(), dims, ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
    
    DBClose(db);
    return 0;
}

template<typename T>
void
gnuplot_dump(const fd_grid<T>& g, int32_t t)
{
    std::stringstream ss;
    ss << "timestep_" << t << ".dat";
    std::ofstream ofs(ss.str());

    ofs << "X Y Z value" << std::endl;
    for (int32_t i = 0; i < g.x_size(); i++)
        for (int32_t j = 0; j < g.y_size(); j++)
            ofs << i << " " << j << " " << "0.0 " << g(i,j) << std::endl;

    ofs.close();
}

int heat(void)
{
    using T = double;

    const int32_t Xdim = 200, Ydim = 200;

    const T alpha = 0.1;
    const T dt = 0.00005;
    const T dx = 1./Xdim;
    const T dy = 1./Ydim;

    const T cx = alpha*dt/(dx*dx);
    const T cy = alpha*dt/(dy*dy);

    std::cout << "CFL: " << cx << " " << cy << std::endl;
    if (cx > 0.5 or cy > 0.5)
    {
        std::cout << "CFL too high" << std::endl;
        return -1;
    }

    fd_grid<T> u(Xdim, Ydim), un(Xdim, Ydim);

    for (int32_t i = 0; i <= Xdim; i++)
    {   
        for (int32_t j = 0; j <= Ydim; j++)
        {
            T x = dx*i - 0.5;
            T y = dy*j - 0.5;
            u(i,j) = 10*std::exp(-20*(x*x + y*y));
        }
    }

    for (int32_t t = 0; t < 400; t++)
    {
        
        for (int32_t i = 0; i <= Xdim; i++)
        {
            
            u(i,-1) = 5;
            u(i, 0) = 5;
            
            //u(i,-1) = u(i,0);
            //u(i,Ydim+1) = u(i,Ydim);

            u(i,Ydim) = -5;
            u(i,Ydim+1) = -5;
        }
        

        for (int32_t j = 0; j <= Ydim; j++)
        {
            u(-1,j) = u(0,j);
            u(Xdim+1,j) = u(Xdim,j);
        }

        if (t%100 == 0)
        {
            std::cout << "Timestep " << t << std::endl;
            visit_dump(u,t);
        }

        for (int32_t i = 0; i <= Xdim; i++)
        {
            for(int32_t j = 0; j <= Ydim; j++)
            {
                un(i,j) = cx*(u(i+1,j) - 2*u(i,j) + u(i-1,j)) +
                          cy*(u(i,j+1) - 2*u(i,j) + u(i,j-1)) +
                          u(i,j);
            }
        }
        
        std::swap(un, u);

    }

    return 0;

}




int wave(void)
{
    using T = double;

    const int32_t Xdim = 100, Ydim = 100;

    auto c_fun = [&](int32_t x, int32_t y) -> auto {
        if (x > Xdim/2)
            return 10;
        return 3;
    };

    const T c = 50; //wave speed
    const T b = 0.1; //damping factor
    const T dt = 0.0005;
    const T dx = 1./Xdim;
    const T dy = 1./Ydim;

    const T cx = dt*dt/(c*dx*dx);
    const T cy = dt*dt/(c*dy*dy);

    std::cout << "CFL: " << cx << " " << cy << std::endl;
    if (cx > 0.5 or cy > 0.5)
    {
        std::cout << "CFL too high" << std::endl;
        return -1;
    }

    fd_grid<T> u(Xdim, Ydim), up(Xdim, Ydim), upp(Xdim, Ydim);

    for (int32_t i = 0; i <= Xdim; i++)
    {   
        for (int32_t j = 0; j <= Ydim; j++)
        {
            T x = dx*i - 0.3;
            T y = dy*j - 0.1;
            upp(i,j) = -std::exp(-2400*(x*x + y*y));
            up(i,j) = 2*dt*upp(i,j);
        }
    }

    auto D1 = 1./(1. + 0.5*b*dt);
    auto D2 = 0.5*b*dt - 1;

    for (int32_t t = 0; t < 100000; t++)
    {
        
        for (int32_t i = 0; i <= Xdim; i++)
        {
            u(i,-1) = 2*dy*u(i,0); /* neumann, without mult is abc */
            u(i,Ydim+1) = 2*dy*u(i,Ydim);
            //u(i, Ydim) = u(i,-1);
            //u(i, 0) = u(i, Ydim);
        }
        

        for (int32_t j = 0; j <= Ydim; j++)
        {
            u(-1,j) = 2*dx*u(0,j);
            u(Xdim+1,j) = 2*dx*u(Xdim,j);
        }

        if (t%100 == 0)
        {
            std::cout << "Timestep " << t << std::endl;
            visit_dump(u,t);
        }

        for (int32_t i = 0; i <= Xdim; i++)
        {
            for(int32_t j = 0; j <= Ydim; j++)
            {
                T f = 1./4.;
                if (i > Xdim/2)
                    f = 1./30.;

                u(i,j) = D1*( D2*upp(i,j) +  2*up(i,j) + f*cx*(up(i+1,j) - 2*up(i,j) + up(i-1,j)) + f*cy*(up(i,j+1) - 2*up(i,j) + up(i,j-1)) );

            }
        }
        
        std::swap(up, upp);
        std::swap(u, up);

    }

    return 0;

}



template<typename T>
int
visit_dump_ns(const scalar_fd_grid<T>& p, const vector_fd_grid<T>& u, int32_t t)
{
    std::stringstream ss;
    ss << "timestep_" << t << ".silo";
    
    DBfile *db = nullptr;
    db = DBCreate(ss.str().c_str(), DB_CLOBBER, DB_LOCAL, "Wave equation", DB_PDB);
    if (!db)
    {
        std::cout << "Cannot write simulation output" << std::endl;
        return -1;
    }
    
    std::vector<double> x(p.x_size()+1);
    std::vector<double> y(p.y_size()+1);
    
    for (size_t i = 0; i < x.size(); i++)
        x.at(i) = double(i)/p.x_size();
        
    for (size_t i = 0; i < y.size(); i++)
        y.at(i) = double(i)/p.y_size();
        
    int dims[] = { int(x.size()), int(y.size())};
    int ndims = 2;
    double *coords[] = {x.data(), y.data()};
    
    DBPutQuadmesh(db, "mesh", NULL, coords, dims, ndims, DB_DOUBLE, DB_COLLINEAR, NULL);
    
    std::vector<double> pv(x.size() * y.size());
    std::vector<double> xv(x.size() * y.size());
    std::vector<double> yv(x.size() * y.size());
    
    for (size_t i = 0; i < x.size(); i++)
    {
        for (size_t j = 0; j < y.size(); j++)
        {
            pv.at(j*x.size()+i) = p(i,j);
            xv.at(j*x.size()+i) = u(i,j).x;
            yv.at(j*x.size()+i) = u(i,j).y;
        }
    }
        
    DBPutQuadvar1(db, "p", "mesh", pv.data(), dims, ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
    DBPutQuadvar1(db, "ux", "mesh", xv.data(), dims, ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
    DBPutQuadvar1(db, "uy", "mesh", yv.data(), dims, ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);

    const char *names[] = {"u", "vorticity"};
    const char *defs[] = {"{ux, uy}", "curl(u)"};
    int types[] = {DB_VARTYPE_VECTOR, DB_VARTYPE_SCALAR};
    DBPutDefvars(db, "defvars", 2, names, types, defs, NULL);
    
    DBClose(db);
    return 0;
}

template<typename T>
T l2_error(const fd_grid<T>& u, const fd_grid<T>& v)
{
    if ( !u.is_compatible_with(v) )
        throw std::invalid_argument("Incompatible grid dimensions");

    const int32_t xdim = u.x_size();    
    const int32_t ydim = u.y_size();    

    const T dx = u.dx();
    const T dy = u.dy();

    const T A = (dx*dy);

    T l2err = 0.0;
    for(int32_t i = 0; i < xdim; i++)
    {
        for(int32_t j = 0; j < ydim; j++)
        {   // 13 mults, 14 adds
            auto d00 = u(i,j) - v(i,j);
            auto d01 = u(i,j+1) - v(i,j+1);
            auto d10 = u(i+1,j) - v(i+1,j);
            auto d11 = u(i+1,j+1) - v(i+1,j+1);
            auto x1 = dx*i; auto y1 = dy*j;
            auto x2 = dx*(i+1); auto y2 = dy*(j+1);
            auto x = x1+0.5*dx; auto y = y1+0.5*dy;

            auto f = (x2-x)*(d00*(y2-y) + d01*(y-y1)) + (x-x1)*(d10*(y2-y) + d11*(y-y1));

            l2err += f*f;
        }
    }

    return std::sqrt(l2err/A);
}

int navier_stokes(void)
{
    using T = double;

    const int32_t Xdim = 300, Ydim = 300;
    const bounding_box<T> bb(0., 0., 2., 2.);

    const T dt = 0.001;
    
    const T rho = 1;
    const T nu = 0.01;
    
    fd_grid<T> ux(Xdim, Ydim, bb), ux_p(Xdim, Ydim, bb);
    fd_grid<T> uy(Xdim, Ydim, bb), uy_p(Xdim, Ydim, bb);
    fd_grid<T> p(Xdim, Ydim, bb), pp(Xdim, Ydim, bb);//, b(Xdim, Ydim, bb);

    const T dx = ux.dx();
    const T dy = ux.dy();

    std::cout << nu*dt/(dx*dy) << std::endl;
    std::cout << 1/nu << std::endl;
    
    for (int32_t t = 0; t < 10; t++)
    {
        std::cout << "Timestep " << t << ", t = " << t*dt << std::endl;
        /*
        for (int32_t i = 0; i <= Xdim; i++)
        {
            for(int32_t j = 0; j <= Ydim; j++)
            {
                auto p1x = (ux(i+1,j) - ux(i-1,j))/(2*dx);
                auto p1y = (uy(i,j+1) - uy(i,j-1))/(2*dy);
                b(i,j) = (p1x + p1y)/dt - p1x*p1x
                       - 2 * ((ux(i,j+1) - ux(i,j-1))/(2*dy)) * ((uy(i+1,j) - uy(i-1,j))/(2*dx))
                       - p1y*p1y;
            }
        }
        */
        
        T err;
        int32_t li;
        for (li = 0; li < 10000; li++)
        {
            std::swap(p, pp);
            for (int32_t i = 0; i <= Xdim; i++)
            {
                for(int32_t j = 0; j <= Ydim; j++)
                {   
                    auto p1x = (ux(i+1,j) - ux(i-1,j))/(2*dx);
                    auto p1y = (uy(i,j+1) - uy(i,j-1))/(2*dy);
                    auto b = (p1x + p1y)/dt - p1x*p1x
                         - 2 * ((ux(i,j+1) - ux(i,j-1))/(2*dy)) * ((uy(i+1,j) - uy(i-1,j))/(2*dx))
                         - p1y*p1y;

                    auto n = dx*dx*dy*dy;
                    auto d = 2*(dx*dx + dy*dy);
                    auto pa = dy*dy*(pp(i+1,j) + pp(i-1,j))/d;     
                    auto pb = dx*dx*(pp(i,j+1) + pp(i,j-1))/d;
                    auto pc = rho * (n/d) * b;
                    
                    p(i,j) = pa + pb - pc;
                }
            }

            apply_neumann(p, boundary::LEFT, 0);
            apply_neumann(p, boundary::RIGHT, 0);
            apply_neumann(p, boundary::BOTTOM, 0);
            apply_dirichlet(p, boundary::TOP, 0);

            err = l2_error(p,pp);
            if (err < 1e-6)
                break;
        }

        std::cout << "Pressure solver iterations at step " << t << ": " << li << " err: " << err << std::endl;
        
        std::swap(ux, ux_p);
        std::swap(uy, uy_p);
        
        for (int32_t i = 0; i <= Xdim; i++)
        {
            for(int32_t j = 0; j <= Ydim; j++)
            {
                ux(i,j) = ux_p(i,j)
                        - ux_p(i,j)*(dt/dx)*(ux_p(i,j) - ux_p(i-1,j))
                        - uy_p(i,j)*(dt/dy)*(ux_p(i,j) - ux_p(i,j-1))
                        - (dt/(2*rho*dx))*(p(i+1,j) - p(i-1,j))
                        + nu * (dt/(dx*dx)) * (ux_p(i+1,j) - 2*ux_p(i,j) + ux_p(i-1,j))
                        + nu * (dt/(dy*dy)) * (ux_p(i,j+1) - 2*ux_p(i,j) + ux_p(i,j-1));
                        
                uy(i,j) = uy_p(i,j)
                        - ux_p(i,j)*(dt/dx)*(uy_p(i,j) - uy_p(i-1,j))
                        - uy_p(i,j)*(dt/dy)*(uy_p(i,j) - uy_p(i,j-1))
                        - (dt/(2*rho*dy))*(p(i,j+1) - p(i,j-1))
                        + nu * (dt/(dx*dx)) * (uy_p(i+1,j) - 2*uy_p(i,j) + uy_p(i-1,j))
                        + nu * (dt/(dy*dy)) * (uy_p(i,j+1) - 2*uy_p(i,j) + uy_p(i,j-1));
                        
            }
        }

        apply_dirichlet(ux, boundary::TOP, 1);
        apply_dirichlet(ux, boundary::LEFT, 0);
        apply_dirichlet(ux, boundary::RIGHT, 0);
        apply_dirichlet(ux, boundary::BOTTOM, 0);
        
        apply_dirichlet(uy, boundary::TOP, 0);
        apply_dirichlet(uy, boundary::LEFT, 0);
        apply_dirichlet(uy, boundary::RIGHT, 0);
        apply_dirichlet(uy, boundary::BOTTOM, 0);

        //if (t%10 == 0)
        //    visit_dump_ns(p,ux,uy,t);
    }

    return 0;

}

int navier_stokes_vgrid(void)
{
    using T = double;

    const int32_t Xdim = 400, Ydim = 400;
    const int32_t Press_iters = 10000;
    const int32_t NS_iters = 1000000;
    const bounding_box<T> bb(0., 0., 1., 1.);

    const T dt = 0.00001;
    
    const T rho = 1;
    const T nu = 0.001;
    
    vector_fd_grid<T> u(Xdim, Ydim, bb), up(Xdim, Ydim, bb);
    scalar_fd_grid<T> p(Xdim, Ydim, bb), pp(Xdim, Ydim, bb);

    const T dx = u.dx();
    const T dy = u.dy();
    const T dxdx = dx*dx;
    const T dydy = dy*dy;
    const T inv_twodx = 1./(2*dx);
    const T inv_twody = 1./(2*dy);
    const T inv_dt = 1./dt;
    const T dt_over_dx = dt/dx;
    const T dt_over_dy = dt/dy;
    const T dt_over_dxdx = dt/(dx*dx);
    const T dt_over_dydy = dt/(dy*dy);
    const T dt_over_2rhodx = dt/(2*rho*dx);
    const T dt_over_2rhody = dt/(2*rho*dy);

    const T n = dxdx*dydy;
    const T d = 2*(dxdx + dydy);
    const T inv_d = 1./d;
    const T n_over_d = n/d;

    std::cout << nu*dt/(dx*dy) << std::endl;
    std::cout << 1/nu << std::endl;
    
    size_t flops = 0;

    for (int32_t t = 0; t < NS_iters; t++)
    {
        std::cout << "Timestep " << t << ", t = " << t*dt << std::endl;
        
        T err;
        int32_t li;
        for (li = 0; li < Press_iters; li++)
        {
            std::swap(p, pp);
            for (int32_t i = 0; i <= Xdim; i++)
            {
                for(int32_t j = 0; j <= Ydim; j++)
                {
                    // 14 mults, 12 adds
                    auto p1x = (u(i+1,j).x - u(i-1,j).x) * inv_twodx;
                    auto p1y = (u(i,j+1).y - u(i,j-1).y) * inv_twody;
                    auto b = (p1x + p1y)*inv_dt - p1x*p1x
                         - 2 * ((u(i,j+1).x - u(i,j-1).x) * inv_twody) * ((u(i+1,j).y - u(i-1,j).y) * inv_twodx)
                         - p1y*p1y;

                    auto pa = dydy*(pp(i+1,j) + pp(i-1,j))*inv_d;     
                    auto pb = dxdx*(pp(i,j+1) + pp(i,j-1))*inv_d;
                    auto pc = rho * n_over_d * b;
                    
                    p(i,j) = pa + pb - pc;
                }
            }

            apply_neumann(p, boundary::LEFT, 0);
            apply_neumann(p, boundary::RIGHT, 0);
            apply_neumann(p, boundary::BOTTOM, 0);
            apply_dirichlet(p, boundary::TOP, 0);

            err = l2_error(p,pp);
            if (err < 1e-7)
                break;
        }

        flops += (26+27)*li*Xdim*Ydim;

        std::cout << "Pressure solver iterations at step " << t << ": " << li << " err: " << err << std::endl;
        
        std::swap(u, up);
        
        for (int32_t i = 0; i <= Xdim; i++)
        {
            for(int32_t j = 0; j <= Ydim; j++)
            {   // 22 mults, 24 adds
                u(i,j).x = up(i,j).x
                         - up(i,j).x*dt_over_dx*(up(i,j).x - up(i-1,j).x)
                         - up(i,j).y*dt_over_dy*(up(i,j).y - up(i,j-1).y)
                         - dt_over_2rhodx*(p(i+1,j) - p(i-1,j))
                         + nu * dt_over_dxdx * (up(i+1,j).x - 2*up(i,j).x + up(i-1,j).x)
                         + nu * dt_over_dydy * (up(i,j+1).x - 2*up(i,j).x + up(i,j-1).x);
                        
                u(i,j).y = up(i,j).y
                         - up(i,j).x*dt_over_dx*(up(i,j).y - up(i-1,j).y)
                         - up(i,j).y*dt_over_dy*(up(i,j).y - up(i,j-1).y)
                         - dt_over_2rhody*(p(i,j+1) - p(i,j-1))
                         + nu * dt_over_dxdx * (up(i+1,j).y - 2*up(i,j).y + up(i-1,j).y)
                         + nu * dt_over_dydy * (up(i,j+1).y - 2*up(i,j).y + up(i,j-1).y);
                        
            }
        }

        flops += 46*Xdim*Ydim;

        apply_dirichlet(u, boundary::TOP, vect<T>(1,0));
        apply_dirichlet(u, boundary::LEFT, vect<T>(0,0));
        apply_dirichlet(u, boundary::RIGHT, vect<T>(0,0));
        apply_dirichlet(u, boundary::BOTTOM, vect<T>(0,0));

        if (t%1000 == 0)
            visit_dump_ns(p,u,t);
    }

    std::cout << "flops: " << flops << std::endl;

    return 0;

}


int main(void)
{
#ifdef VGRID
    navier_stokes_vgrid();
#else
    navier_stokes();
#endif
}







