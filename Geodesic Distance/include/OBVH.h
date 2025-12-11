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

    void add(const Vec3& p)
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
    Vec3<T> origin;
    Vec3<T> direction;

    Ray(const Vec3& ori, const Vec3& dir) :origin(ori), direction(dir) {}
    
};