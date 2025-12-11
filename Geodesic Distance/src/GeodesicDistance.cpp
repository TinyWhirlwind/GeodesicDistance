#include "GeodesicDistance.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
class GeodesicDistance::PImpl
{
public:
    PImpl(std::shared_ptr<OMesh> mesh):_mesh(mesh)
    {
        _time = 0.0f;
        _init = true;
    }
    ~PImpl() {}
public:
    void GetLaplaceMat();
    void BuildU0(int source);
    void SolveHeatFlow();
    void ComputeGradU();
    void ComputeNormalizedX();
    void ComputeDivergence();
    void SolvePoisson();
    void NormalizePhi(int source);
    void Show();

public:
    bool _init = false;
    std::shared_ptr<OMesh> _mesh;
    float _time = 0.0;
    Eigen::VectorXd _u0;
    Eigen::VectorXd _u;
    std::vector<Eigen::Vector3d> _gradU;
    std::vector<Eigen::Vector3d> _X;
    Eigen::Vector3d _div;
    Eigen::VectorXd _phi;
    Eigen::SparseMatrix<double> _A;
    Eigen::SparseMatrix<double> _Lc;

};
GeodesicDistance::GeodesicDistance(std::shared_ptr<OMesh> mesh, OPoint sourcePoint)
{
    impl_.reset(new PImpl(mesh));
}

GeodesicDistance::~GeodesicDistance()
{
}

void GeodesicDistance::PImpl::GetLaplaceMat()
{
    size_t vn = _mesh->n_vertices();
    _A.resize(vn, vn);
    _Lc.resize(vn, vn);

    std::vector<Eigen::Triplet<double>> tripletsA;
    std::vector<Eigen::Triplet<double>> tripletsLc;
    for (int i = 0; i < vn; ++i)
    {
        float Ai = 0.0;
        auto v_it = _mesh->vertex_handle(i);
        for (auto vf_it = _mesh->vf_iter(v_it); !vf_it->is_valid(); ++vf_it)
        {
            float area = _mesh->calc_face_area(*vf_it);
            Ai += area;
        }
        Ai /= 3.0;
        tripletsA.emplace_back(i,i, Ai);
    }

    float averge_edge_length = 0.0;
    for (auto e_it = _mesh->edges_begin(); e_it != _mesh->edges_end(); ++e_it)
    {
        double Wij = 0.0;
        for (int i = 0; i < 2; ++i)
        {
            auto he_it = _mesh->halfedge_handle(*e_it, i);
            if (!he_it.is_valid())continue;
            auto next_he_it = _mesh->next_halfedge_handle(he_it);
            auto p0 = _mesh->point(_mesh->to_vertex_handle(he_it));
            auto p1 = _mesh->point(_mesh->from_vertex_handle(he_it));
            auto p2 = _mesh->point(_mesh->from_vertex_handle(next_he_it));
            auto v0 = p1 - p0;
            auto v1 = p2 - p0;
            Wij += v0.dot(v1) / v0.cross(v1).norm();
        }
        Wij /= 2;
        int pid0 = e_it->v0().idx();
        int pid1 = e_it->v1().idx();
        tripletsLc.emplace_back( pid0,pid1,-Wij );
        tripletsLc.emplace_back(pid1,pid0,-Wij);
        /*tripletsLc.push_back({ pid0,pid0,Wij });
        tripletsLc.push_back({ pid1,pid1,Wij });*/
        
        float per_edge_length = _mesh->calc_edge_length(*e_it);
        averge_edge_length += per_edge_length;
    }
    averge_edge_length /= _mesh->n_edges();
    _time = averge_edge_length * averge_edge_length;

    _A.setFromTriplets(tripletsA.begin(), tripletsA.end());
    _Lc.setFromTriplets(tripletsLc.begin(), tripletsLc.end());
}

void GeodesicDistance::PImpl::BuildU0(int source)
{
    size_t vn = _mesh->n_vertices();
    _u0 = Eigen::VectorXd::Zero(vn);
    _u0[source] = 1.0;
}

void GeodesicDistance::PImpl::SolveHeatFlow()
{
    Eigen::SparseMatrix<double> M = _A - _time * _Lc;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(M);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Heat Solve failed\n";
        return;
    }

    _u = solver.solve(_A * _u0);  // 顶点标量
}

auto transfer = [](OPoint p)->Eigen::Vector3d
    {
        return Eigen::Vector3d(p[0], p[1], p[2]);
    };

void GeodesicDistance::PImpl::ComputeGradU()
{
    size_t fn = _mesh->n_faces();
    _gradU.resize(fn);
    for (int i = 0; i < fn; ++i)
    {
        auto f_it = _mesh->face_handle(i);
        double u[3];
        OPoint p[3];
        for (auto fv_it = _mesh->fv_iter(f_it); !fv_it->is_valid(); ++fv_it)
        {
            int vid = fv_it->idx();
            u[i] = _u(vid);
            p[i] = _mesh->point(*fv_it);
        }

        Eigen::Vector3d e0 = transfer(p[1] - p[0]);
        Eigen::Vector3d e1 = transfer(p[2] - p[1]);
        Eigen::Vector3d e2 = transfer(p[0] - p[2]);

        Eigen::Vector3d n = e0.cross(-e2);
        double area = n.norm() / 2.0;
        Eigen::Vector3d grad_u = (n.cross(e0)*u[0] + n.cross(e1) * u[1] + n.cross(e2) * u[2]) / (2.0 * area);
        _gradU[f_it.idx()] = grad_u;
    }
}

void GeodesicDistance::PImpl::ComputeNormalizedX()
{
    size_t fn = _mesh->n_faces();
    _X.resize(fn);

    for (size_t i = 0; i < fn; ++i)
    {
        Eigen::Vector3d g = _gradU[i];

        if (g.norm() > 1e-10)
            _X[i] = -g.normalized();    // 负梯度方向 = geodesic 方向场
        else
            _X[i] = Eigen::Vector3d::Zero();
    }
}

//计算散度
void GeodesicDistance::PImpl::ComputeDivergence()
{
    size_t vn = _mesh->n_vertices();
    _div = Eigen::VectorXd::Zero(vn);

    for (auto v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); ++v_it)
    {
        double div_val = 0.0;
        for (auto vf_it = _mesh->vf_iter(*v_it); vf_it.is_valid(); ++vf_it)
        {
            int fid = vf_it->idx();
            Eigen::Vector3d X = _X[fid];
            Eigen::Vector3d normal = transfer(_mesh->calc_face_normal(*vf_it));
            div_val += X.dot(normal) * 0.5;
        }

        _div[v_it->idx()] = div_val;
    }
}

void GeodesicDistance::PImpl::SolvePoisson()
{
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(_Lc);

    if (solver.info() != Eigen::Success)
    {
        std::cerr << "Poisson solve failed!\n";
        return;
    }

    _phi = solver.solve(_div);
}

void GeodesicDistance::PImpl::NormalizePhi(int source)
{
    double offset = _phi[source];
    _phi.array() -= offset;
}




