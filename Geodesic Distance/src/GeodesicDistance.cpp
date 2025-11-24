#include "GeodesicDistance.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
class GeodesicDistance::PImpl
{
public:
    PImpl(CMeshO& mesh)
    {
        _mesh = mesh;
        _time = 0.0f;
        _mesh.vert.EnableMark();
        _mesh.face.EnableMark();
        vcg::tri::UnMarkAll(_mesh);
        vcg::tri::InitFaceIMark(_mesh);
        vcg::tri::InitVertexIMark(_mesh);
        _init = true;
    }
    ~PImpl() {}
public:
    void SetTime();
    void GetLaplaceMat();
public:
    bool _init = false;
    CMeshO _mesh;
    float _time = 0.0;
};
GeodesicDistance::GeodesicDistance(CMeshO& mesh, CVertexO* sourcePoint)
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
    for (auto fit = _mesh.face.begin(); fit != _mesh.face.end(); ++fit)
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
    _time = avergeEdgeLength * avergeEdgeLength;
    _mesh.ClearAttributes();
}

void GeodesicDistance::PImpl::GetLaplaceMat()
{
    int vn = _mesh.vn;
    Eigen::SparseMatrix<double> K(vn, vn); // 刚度矩阵
    Eigen::SparseMatrix<double> M(vn, vn); // 质量矩阵

    std::vector<Eigen::Triplet<double>> tripletsK;
    std::vector<Eigen::Triplet<double>> tripletsM;
    for (auto fit = _mesh.face.begin(); fit != _mesh.face.end(); ++fit)
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

        }
    }
}


