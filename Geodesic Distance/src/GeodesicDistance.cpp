#include "GeodesicDistance.h"

class GeodesicDistance::PImpl
{
    public:
        PImpl(CMeshO& mesh)
        {
            _mesh = mesh;
            _avergeEdgeLength = 0.0f;
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
    void computeAverageEdgeLength();

public:
    bool _init = false;
    float _avergeEdgeLength = 0.0f;
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

void GeodesicDistance::PImpl::computeAverageEdgeLength()
{
    if (!_init)return;
    int enb = 0;
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
            _avergeEdgeLength += edgeLength;
            enb++;
        }
    }
    _avergeEdgeLength /= enb;
}

void GeodesicDistance::SetTime()
{
    impl_->_time = impl_->_avergeEdgeLength;
}