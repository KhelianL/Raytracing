#ifndef Sphere_H
#define Sphere_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySphereIntersection
{
    bool intersectionExists;
    float t;
    float theta, phi;
    Vec3 intersection;
    Vec3 secondintersection;
    Vec3 normal;
};

Vec3 SphericalCoordinatesToEuclidean(Vec3 ThetaPhiR)
{
    return ThetaPhiR[2] * Vec3(cos(ThetaPhiR[0]) * cos(ThetaPhiR[1]), sin(ThetaPhiR[0]) * cos(ThetaPhiR[1]), sin(ThetaPhiR[1]));
}
Vec3 SphericalCoordinatesToEuclidean(float theta, float phi)
{
    return Vec3(cos(theta) * cos(phi), sin(theta) * cos(phi), sin(phi));
}

static Vec3 EuclideanCoordinatesToSpherical(Vec3 xyz)
{
    float R = xyz.length();
    float phi = asin(xyz[2] / R);
    float theta = atan2(xyz[1], xyz[0]);
    return Vec3(theta, phi, R);
}

class Sphere : public Mesh
{
public:
    Vec3 m_center;
    float m_radius;

    Sphere() : Mesh() {}
    Sphere(Vec3 c, float r) : Mesh(), m_center(c), m_radius(r) {}

    void build_arrays()
    {
        unsigned int nTheta = 20, nPhi = 20;
        positions_array.resize(3 * nTheta * nPhi);
        normalsArray.resize(3 * nTheta * nPhi);
        uvs_array.resize(2 * nTheta * nPhi);
        for (unsigned int thetaIt = 0; thetaIt < nTheta; ++thetaIt)
        {
            float u = (float)(thetaIt) / (float)(nTheta - 1);
            float theta = u * 2 * M_PI;
            for (unsigned int phiIt = 0; phiIt < nPhi; ++phiIt)
            {
                unsigned int vertexIndex = thetaIt + phiIt * nTheta;
                float v = (float)(phiIt) / (float)(nPhi - 1);
                float phi = -M_PI / 2.0 + v * M_PI;
                Vec3 xyz = SphericalCoordinatesToEuclidean(theta, phi);
                positions_array[3 * vertexIndex + 0] = m_center[0] + m_radius * xyz[0];
                positions_array[3 * vertexIndex + 1] = m_center[1] + m_radius * xyz[1];
                positions_array[3 * vertexIndex + 2] = m_center[2] + m_radius * xyz[2];
                normalsArray[3 * vertexIndex + 0] = xyz[0];
                normalsArray[3 * vertexIndex + 1] = xyz[1];
                normalsArray[3 * vertexIndex + 2] = xyz[2];
                uvs_array[2 * vertexIndex + 0] = u;
                uvs_array[2 * vertexIndex + 1] = v;
            }
        }
        triangles_array.clear();
        for (unsigned int thetaIt = 0; thetaIt < nTheta - 1; ++thetaIt)
        {
            for (unsigned int phiIt = 0; phiIt < nPhi - 1; ++phiIt)
            {
                unsigned int vertexuv = thetaIt + phiIt * nTheta;
                unsigned int vertexUv = thetaIt + 1 + phiIt * nTheta;
                unsigned int vertexuV = thetaIt + (phiIt + 1) * nTheta;
                unsigned int vertexUV = thetaIt + 1 + (phiIt + 1) * nTheta;
                triangles_array.push_back(vertexuv);
                triangles_array.push_back(vertexUv);
                triangles_array.push_back(vertexUV);
                triangles_array.push_back(vertexuv);
                triangles_array.push_back(vertexUV);
                triangles_array.push_back(vertexuV);
            }
        }
    }

    RaySphereIntersection intersect(const Ray &ray) const
    {
        RaySphereIntersection intersection;

        Vec3 oc = ray.origin() - m_center;
        float a = Vec3::dot(ray.direction(), ray.direction());
        float b = 2.0 * Vec3::dot(oc, ray.direction());
        float c = Vec3::dot(oc, oc) - m_radius * m_radius;
        float discriminant = b * b - 4 * a * c;

        // 2 Solutions
        if (discriminant > 0)
        {
            intersection.intersectionExists = true;
            float soluce1 = (-b - sqrt(discriminant)) / (2.0 * a);
            float soluce2 = (-b + sqrt(discriminant)) / (2.0 * a);

            intersection.t = (soluce1 > soluce2 ? soluce2 : soluce1);
            intersection.intersection = ray.origin() + (ray.direction() * intersection.t);
            intersection.secondintersection = ray.origin() + (ray.direction() * (soluce1 > soluce2 ? soluce1 : soluce2));
            intersection.normal = intersection.intersection - m_center;
            intersection.normal.normalize();

            // Coords Polaires
            Vec3 res = EuclideanCoordinatesToSpherical(intersection.intersection);
            intersection.theta = res[0];
            intersection.phi = res[1];
        }

        // 1 Solution
        if (discriminant == 0)
        {
            intersection.intersectionExists = true;
            intersection.t = (-b / (2.0 * a));
            intersection.intersection = ray.origin() + (ray.direction() * intersection.t);
            intersection.normal = intersection.intersection - m_center;
            intersection.normal.normalize();

            // Coords Polaires
            Vec3 res = EuclideanCoordinatesToSpherical(intersection.intersection);
            intersection.theta = res[0];
            intersection.phi = res[1];
        }

        // 0 Solution
        if (discriminant < 0)
        {
            intersection.intersectionExists = false;
        }
        return intersection;
    }
};
#endif
