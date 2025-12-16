#include "mymesh.h"
#include <cmath>
template<typename T>
struct Vec3
{
    T x, y, z;
    Vec3() :x(0), y(0), z(0) {}
    Vec3(T _x, T _y, T _z) :x(_x), y(_y), z(_z) {}

    Vec3 operator+ (const Vec3& v) const {
        return Vec3(x + v.x, y + v.y, z + v.z);
    }

    Vec3 operator-(const Vec3& v) const {
        return Vec3(x - v.x, y - v.y, z - v.z);
    }

    Vec3 operator*(T s) const {
        return Vec3(x * s, y * s, z * s);
    }

    Vec3 operator/(T s) const {
        return Vec3(x / s, y / s, z / s);
    }

    Vec3 operator-() const {
        return Vec3(-x, -y, -z);
    }

    T dot(const Vec3& v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    Vec3 cross(const Vec3& v) const {
        return Vec3(
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x
        );
    }
    T length() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    Vec3 normalized() const {
        T len = length();
        T eps = std::numeric_limits<T>::epsilon();
        if (len < eps) 
            return Vec3(0, 0, 0);
        return Vec3(x / len, y / len, z / len);
    }

    void normalize() {
        T len = length();
        T eps = std::numeric_limits<T>::epsilon();
        if (len < eps)
        {
            x = 0; y = 0; z = 0;
            return;
        }
        x /= len;
        y /= len;
        z /= len);
    }
};

template<typename T>
struct ABBox {

    Vec3<T> min;
    Vec3<T> max;
    ABBox() :min(Vec3(std::numeric_limits<T>::max(), std::numeric_limits<T>::max(), std::numeric_limits<T>::max())),
        max(Vec3(std::numeric_limits<T>::min(), std::numeric_limits<T>::min(), std::numeric_limits<T>::min())) {
    };
    ABBox(const Vec3& min, const Vec3& max) :min(min), max(max) {}

    void expand(const Vec3& p)
    {
        min.x = std::min(min.x, p.x);
        min.y = std::min(min.y, p.y);
        min.z = std::min(min.z, p.z);
        max.x = std::max(max.x, p.x);
        max.y = std::max(max.y, p.y);
        max.z = std::max(max.z, p.z);
    }

    void merge(const ABBox& other)
    {
        min.x = std::min(min.x, other.min.x);
        min.y = std::min(min.y, other.min.y);
        min.z = std::min(min.z, other.min.z);

        max.x = std::max(max.x, other.max.x);
        max.y = std::max(max.y, other.max.y);
        max.z = std::max(max.z, other.max.z);
    }

    //sah cost
    T surfaceArea() const {
        Vec3 d = max - min;
        return 2 * (d.x * d.y + d.y * d.z + d.z * d.x);
    }

    Vec3 center() const
    {
        return (min + max) / 2; 
    }
};

template<typename T>
struct Ray
{
    Vec3 origin;
    Vec3 direction;
    Vec3 invDirection;

    Ray(const Vec3& ori, const Vec3& dir) :origin(ori), direction(dir) 
    {
        invDirection.x = 1.0f / direction.x;
        invDirection.y = 1.0f / direction.y;
        invDirection.z = 1.0f / direction.z;
    };

    Vec3 pointAt(T t)
    {
        return origin + direction * t;
    }
};

template<typename T>
struct Triangle 
{
    Vec3 v0, v1, v2;
    int id;//vcg/openmesh中的faceId
    Vec3 centroid;
    ABBox bbox;

    Triangle(const Vec3& a, const Vec3& b, const Vec3& c, int id = 0) :v0(a), v1(b), v2(c), id(id) 
    {
        centroid = (v0 + v1 + v2) * (1.0f / 3.0f);
        bbox = ABBox();
        bbox.expand(v0);
        bbox.expand(v1);
        bbox.expand(v2);
    };
    
    const ABBox& Box()->const
    {
        return bbox;
    }

    const Vec3& centroid()->const
    {
        return centroid;
    }

    bool intersect(const Ray& ray,float& t,float& u,float& v) const 
    {
        const float EPSILON = 1e-6f;

        Vec3 edge1 = v1 - v0;
        Vec3 edge2 = v2 - v0;
        Vec3 h = ray.direction.cross(edge2);
        float a = edge1.dot(h);
        if (std::fabs(a) < EPSILON) return false;

        float f = 1.0f / a;
        Vec3 s = ray.origin - v0;
        u = f * s.dot(h);
        if (u < 0.0f || u > 1.0f) return false;

        Vec3 q = s.cross(edge1);
        v = f * ray.direction.dot(q);
        if (v < 0.0f || u + v > 1.0f) return false;

        t = f * edge2.dot(q);
        return t > EPSILON;
    }
};

template<typename T>
struct BVHNode
{
    ABBox<T> bbox;
    std::unique_ptr<BVHNode> left;
    std::unique_ptr<BVHNode> right;
    int start = 0;
    int count = 0;

    BVHNode() :start(0), count(0);
    bool isLeaf() const {
        return (left == nullptr) && (right == nullptr);
    }
};

enum class BVHSplitMode {
    MIDDLE,
    EQUAL_COUNT,
    SAH
};

struct BVHBuildConfig {
    BVHSplitMode mode = BVHSplitMode::MIDDLE;
    int sahThreshold = 20000;
    int maxLeafSize = 4;
    int numBuckets = 12;
};

template<typename T>
class BVH 
{
public:
    explicit BVH(std::vector<Triangle>&& tris, BVHBuildConfig cfg = {}):triangles(std::move(tris)),config(cfg)
    {
        primIndices.resize(triangles.size());
        for (int i = 0; i < (int)triangles.size(); ++i)
        {
            primIndices[i] = i;
        }
        if (!triangles.empty())
        {
            buildMode = ((int)triangles.size() >= config.sahThreshold)
                ? BVHSplitMode::SAH_Bucket
                : config.defaultMode;

            root = buildRecursive(0, (int)triangles.size(), 0);
        }
    }
private:
    std::unique_ptr<BVHNode> root;
    std::vector<Trangle> triangles;
    std::vector<int> primIndices;
    BVHBuildConfig config;
    BVHSplitMode buildMode = BVHSplitMode::Median;

private:
    static int longestAxis(const ABBox<T>& box)
    {
        Vec3<T> d = box.max - box.min;
        if (d.x > d.y && d.x > d.z)
            return 0;
        else if (d.y > d.z)
            return 1;
        else
            return 2;
    }

    std::unique_ptr<BVHNode> buidlRecursive(int start, int end, int depth)
    {
        auto node = std::make_unique<BVHNode>();
        for (int i = start; i < end; ++i)
        {
            node->bbox.merge(triangles[primIndices[i]].Box());
        }
        int count = end - start;
        if (count <= config.maxLeafSize)
        {
            node->start = start;
            node->count = count;
            return node;
        }


    }
};