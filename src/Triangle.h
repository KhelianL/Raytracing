#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Vec3.h"
#include "Ray.h"
#include "Plane.h"

using namespace std;

struct RayTriangleIntersection
{
    bool intersectionExists;
    float t;
    float w0, w1, w2;
    unsigned int tIndex;
    Vec3 intersection;
    Vec3 normal;
};

class Triangle
{
private:
    Vec3 m_c[3], m_normal;
    float area;
    Plane m_plane;

public:
    Triangle() {}
    Triangle(Vec3 const &c0, Vec3 const &c1, Vec3 const &c2)
    {
        m_c[0] = c0;
        m_c[1] = c1;
        m_c[2] = c2;
        updateAreaAndNormal();
        m_plane = Plane(m_c[0], m_normal);
    }
    void updateAreaAndNormal()
    {
        Vec3 nNotNormalized = Vec3::cross(m_c[1] - m_c[0], m_c[2] - m_c[0]);
        float norm = nNotNormalized.length();
        m_normal = nNotNormalized / norm;
        m_normal.normalize();
        area = norm / 2.f;
    }
    void setC0(Vec3 const &c0) { m_c[0] = c0; } // remember to update the area and normal afterwards!
    void setC1(Vec3 const &c1) { m_c[1] = c1; } // remember to update the area and normal afterwards!
    void setC2(Vec3 const &c2) { m_c[2] = c2; } // remember to update the area and normal afterwards!
    Vec3 const &normal() const { return m_normal; }

    float distanceToSupportPlane(Vec3 const &p) const
    {
        return sqrt(squareDistanceToSupportPlane(p));
    }

    Vec3 projectOnSupportPlane(Vec3 const &p) const
    {
        return m_plane.project(p);
    }

    float squareDistanceToSupportPlane(Vec3 const &p) const
    {
        return m_plane.squareDistance(p);
    }

    bool isParallelTo(Line const &L) const
    {
        return m_plane.isParallelTo(L);
    }

    Vec3 getIntersectionPointWithSupportPlane(Line const &L) const
    {
        return m_plane.getIntersectionPoint(L);
    }

    void computeBarycentricCoordinates(Vec3 const &p, float &u0, float &u1, float &u2) const
    {
        u0 = Vec3::dot(m_normal, Vec3::cross(m_c[1] - m_c[0], p - m_c[0]));
        u1 = Vec3::dot(m_normal, Vec3::cross(m_c[2] - m_c[1], p - m_c[1]));
        u2 = Vec3::dot(m_normal, Vec3::cross(m_c[0] - m_c[2], p - m_c[2]));
        // ((u0 > 0) && (u1 > 0) && (u2 > 0)) ? "P in triangle" : "P not in Triangle")
    }

    RayTriangleIntersection getIntersection(Ray const &ray) const
    {
        RayTriangleIntersection result;
        result.intersectionExists = false;

        // Test isParallelTo()
        if (!isParallelTo(ray))
        {
            // Test triangle devant le rayon (t >= 0)
            float t = Vec3::dot(m_c[0] - ray.origin(), m_normal) / Vec3::dot(ray.direction(), m_normal);
            if (t >= 0)
            {
                // Test intersection triangle avec computeBarycentricCoordinates()
                float u0, u1, u2;
                Vec3 intersect = getIntersectionPointWithSupportPlane(ray);
                computeBarycentricCoordinates(intersect, u0, u1, u2);

                // Si (u0 > 0) && (u1 > 0) && (u2 > 0) alors intersection
                if ((u0 >= 0) && (u1 >= 0) && (u2 >= 0))
                {
                    result.intersectionExists = true;
                    result.intersection = intersect;
                    result.w0 = u0;
                    result.w1 = u1;
                    result.w2 = u2;
                    result.t = t;
                    result.normal = m_normal;
                }
            }
        }
        return result;
    }
};
#endif
