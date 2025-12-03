#ifndef GEODESIC_DISTANCE_H
#define GEODESIC_DISTANCE_H
#include <mymesh.h>
typedef OMesh::Point OPoint;
typedef OMesh::VertexHandle OVertex;
typedef OMesh::FaceHandle OFace;
class GeodesicDistance
{
public:
    GeodesicDistance(std::shared_ptr<OMesh> mesh, OPoint sourcePoint);
    ~GeodesicDistance();

    void SetTime();
private:
    class PImpl;
    std::shared_ptr<PImpl> impl_;
};
#endif