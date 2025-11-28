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
    void SetTime();
    void GetLaplaceMat();
public:
    bool _init = false;
    std::shared_ptr<OMesh> _mesh;
    float _time = 0.0;
};
GeodesicDistance::GeodesicDistance(std::shared_ptr<OMesh> mesh, OPoint sourcePoint)
{
    impl_.reset(new PImpl(mesh));
}

GeodesicDistance::~GeodesicDistance()
{
}

void GeodesicDistance::SetTime()
{
    impl_->SetTime();
}

void GeodesicDistance::PImpl::SetTime()
{
    if (!_init)return;
    int enb = 0;
    double avergeEdgeLength = 0.0;
    /*for (auto fit = _mesh.face.begin(); fit != _mesh.face.end(); ++fit)
    {
        if (fit->IsD())
        {
            continue;
        }
        float edgeLength = 0.0;
        for (int i = 0; i < 3; ++i)
        {
            auto& v0 = fit->V(i);
            auto& v1 = fit->V((i + 1) % 3);
            if (v0->IsV() && v1->IsV())
            {
                continue;
            }
            edgeLength = (v0->P() - v1->P()).Norm();
            avergeEdgeLength += edgeLength;
            enb++;
        }
    }
    avergeEdgeLength /= enb;
    _time = avergeEdgeLength * avergeEdgeLength;*/
}

void GeodesicDistance::PImpl::GetLaplaceMat()
{
    int vn = _mesh->n_faces();
    Eigen::SparseMatrix<double> A(vn, vn); // 质量矩阵
    Eigen::SparseMatrix<double> Lc(vn, vn); // 权重矩阵

    std::vector<Eigen::Triplet<double>> tripletsA;
    std::vector<Eigen::Triplet<double>> tripletsLc;
    for (int i = 0; i < vn; ++i)
    {
        float Ai = 0.0;
        auto v_it = _mesh->vertex_handle(i);
        for (auto vf_it = _mesh->vf_iter(v_it); vf_it->is_valid(); ++vf_it)
        {
            float area = _mesh->calc_face_area(vf_it);
            Ai += area;
        }
        Ai /= 3.0;
        tripletsA.push_back({ i,i, Ai });
    }

    float averge_edge_length = 0.0;
    for (auto e_it = _mesh->edges_begin(); e_it != _mesh->edges_end(); ++e_it)
    {
        double Wij = 0.0;
        for (int i = 0; i < 2; ++i)
        {
            auto he_it = _mesh->halfedge_handle(e_it, i);
            if (!he_it.is_valid())continue;
            auto next_he_it = _mesh->next_halfedge_handle(he_it);
            OPoint p0 = _mesh->point(_mesh->to_vertex_handle(he_it));
            OPoint p1 = _mesh->point(_mesh->from_vertex_handle(he_it));
            OPoint p2 = _mesh->point(_mesh->from_vertex_handle(next_he_it));
            Wij += ((p1 - p0).dot(p2 - p0))/((p1 - p0).cross(p2 - p0)).norm();
        }
        Wij /= 2;
        int pid0 = e_it->v0().idx();
        int pid1 = e_it->v1().idx();
        tripletsLc.push_back({ pid0,pid1,-Wij });
        tripletsLc.push_back({ pid1,pid0,-Wij });
        /*tripletsLc.push_back({ pid0,pid0,Wij });
        tripletsLc.push_back({ pid1,pid1,Wij });*/
        
        float per_edge_length = _mesh->calc_edge_length(e_it);
        averge_edge_length += per_edge_length;
    }
    averge_edge_length /= vn;
    _time = averge_edge_length * averge_edge_length;
}


