#ifndef GEODESIC_DISTANCE_H
#define GEODESIC_DISTANCE_H
#include <mymesh.h>
class GeodesicDistance
{
public:
    GeodesicDistance(CMeshO& mesh, CVertexO* sourcePoint);
    ~GeodesicDistance();

    void SetTime();
private:
    class PImpl;
    std::shared_ptr<PImpl> impl_;
};
#endif