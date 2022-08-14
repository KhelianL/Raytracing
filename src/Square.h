#ifndef SQUARE_H
#define SQUARE_H

#include <cmath>
#include <vector>
#include "Mesh.h"
#include "Vec3.h"

struct RaySquareIntersection
{
    bool intersectionExists;
    float t;
    float u, v;
    Vec3 intersection;
    Vec3 normal;
};

class Square : public Mesh
{
public:
    Vec3 m_normal;
    Vec3 m_bottom_left;
    Vec3 m_right_vector;
    Vec3 m_up_vector;

    Square() : Mesh() {}
    Square(Vec3 const &bottomLeft, Vec3 const &rightVector, Vec3 const &upVector, float width = 1., float height = 1.,
           float uMin = 0.f, float uMax = 1.f, float vMin = 0.f, float vMax = 1.f) : Mesh()
    {
        setQuad(bottomLeft, rightVector, upVector, width, height, uMin, uMax, vMin, vMax);
    }

    void setQuad(Vec3 const &bottomLeft, Vec3 const &rightVector, Vec3 const &upVector, float width = 1., float height = 1.,
                 float uMin = 0.f, float uMax = 1.f, float vMin = 0.f, float vMax = 1.f)
    {
        m_right_vector = rightVector;
        m_up_vector = upVector;
        m_normal = Vec3::cross(rightVector, upVector);
        m_bottom_left = bottomLeft;

        m_normal.normalize();
        m_right_vector.normalize();
        m_up_vector.normalize();

        m_right_vector = m_right_vector * width;
        m_up_vector = m_up_vector * height;

        vertices.clear();
        vertices.resize(4);
        vertices[0].position = bottomLeft;
        vertices[0].u = uMin;
        vertices[0].v = vMin;
        vertices[1].position = bottomLeft + m_right_vector;
        vertices[1].u = uMax;
        vertices[1].v = vMin;
        vertices[2].position = bottomLeft + m_right_vector + m_up_vector;
        vertices[2].u = uMax;
        vertices[2].v = vMax;
        vertices[3].position = bottomLeft + m_up_vector;
        vertices[3].u = uMin;
        vertices[3].v = vMax;
        vertices[0].normal = vertices[1].normal = vertices[2].normal = vertices[3].normal = m_normal;
        triangles.clear();
        triangles.resize(2);
        triangles[0][0] = 0;
        triangles[0][1] = 1;
        triangles[0][2] = 2;
        triangles[1][0] = 0;
        triangles[1][1] = 2;
        triangles[1][2] = 3;
    }

    RaySquareIntersection intersect(const Ray &ray) const
    {
        // Init
        RaySquareIntersection intersection;
        intersection.intersectionExists = false;
        intersection.t = FLT_MAX;

        Triangle t1 = Triangle(vertices[0].position, vertices[1].position, vertices[2].position);
        Triangle t2 = Triangle(vertices[0].position, vertices[2].position, vertices[3].position);
        RayTriangleIntersection dataT1 = t1.getIntersection(ray);
        RayTriangleIntersection dataT2 = t2.getIntersection(ray);

        if (dataT1.intersectionExists)
        {
            intersection.intersectionExists = true;
            intersection.t = dataT1.t;
            intersection.u = dataT1.w0;
            intersection.v = dataT1.w1;
            intersection.intersection = dataT1.intersection;
            intersection.normal = dataT1.normal;
        }
        if (dataT2.intersectionExists)
        {
            intersection.intersectionExists = true;
            intersection.t = dataT2.t;
            intersection.u = dataT2.w0;
            intersection.v = dataT2.w1;
            intersection.intersection = dataT2.intersection;
            intersection.normal = dataT2.normal;
        }

        return intersection;
    }
};
#endif // SQUARE_H
