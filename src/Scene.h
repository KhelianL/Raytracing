#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include <GL/glut.h>
#include <stdlib.h>
#include <cmath>

#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"
#include "raytracingParam.h"

using namespace std;

enum LightType
{
    LightType_Spherical,
    LightType_Quad
};

struct Light
{
    Vec3 material;
    bool isInCamSpace;
    LightType type;

    Vec3 pos;
    float radius;

    Square quad;

    float powerCorrection;

    Light() : powerCorrection(1.0) {}

    vector<Vec3> randPointSphere()
    {
        vector<Vec3> vectPointsLights(PARAM_NUMBER_SHADOWS);

        if (type == LightType::LightType_Spherical)
        {
            for (int i = 0; i < (int)PARAM_NUMBER_SHADOWS; i++)
            {
                double randPhi = ((((double)(rand() % 100)) / 100.0) * M_PI) - (M_PI / 2);
                double randTheta = (((double)(rand() % 100)) / 100.0) * (2.0 * M_PI);

                Vec3 tmpPoint;
                tmpPoint[0] = pos[0] + radius * sin(randPhi) * cos(randTheta);
                tmpPoint[1] = pos[1] + radius * sin(randPhi) * sin(randTheta);
                tmpPoint[2] = pos[2] + radius * sin(randPhi);
                vectPointsLights.push_back(tmpPoint);
            }
        }
        else if (type == LightType::LightType_Quad)
        {
            for (int i = 0; i < (int)PARAM_NUMBER_SHADOWS; i++)
            {
                double randUnitS = ((double)(rand() % 100)) / 100.0; // [0,1]
                double randUnitT = ((double)(rand() % 100)) / 100.0; // [0,1]

                Vec3 ptE = randUnitS * quad.vertices[0].position + (1 - randUnitS) * quad.vertices[1].position;
                Vec3 ptF = randUnitS * quad.vertices[2].position + (1 - randUnitS) * quad.vertices[3].position;

                Vec3 tmpPoint;
                tmpPoint = randUnitT * ptE + (1 - randUnitT) * ptF;
                vectPointsLights.push_back(tmpPoint);
            }
        }
        return vectPointsLights;
    }
};

enum TypeObject
{
    sphere,
    square,
    mesh
};

struct RaySceneIntersection
{
    bool intersectionExists;
    TypeObject typeOfIntersectedObject;
    unsigned int objectIndex;
    float t;
    RayTriangleIntersection rayMeshIntersection;
    RaySphereIntersection raySphereIntersection;
    RaySquareIntersection raySquareIntersection;
    RaySceneIntersection() : intersectionExists(false), t(FLT_MAX) {}
};

class Scene
{
    std::vector<Mesh> meshes;
    std::vector<Sphere> spheres;
    std::vector<Square> squares;
    std::vector<Light> lights;

public:
    Scene()
    {
    }

    void draw()
    {
        // iterer sur l'ensemble des objets, et faire leur rendu :
        for (unsigned int It = 0; It < meshes.size(); ++It)
        {
            Mesh const &mesh = meshes[It];
            mesh.draw();
        }
        for (unsigned int It = 0; It < spheres.size(); ++It)
        {
            Sphere const &sphere = spheres[It];
            sphere.draw();
        }
        for (unsigned int It = 0; It < squares.size(); ++It)
        {
            Square const &square = squares[It];
            square.draw();
        }
    }

    RaySceneIntersection computeIntersection(Ray const &ray)
    {
        RaySceneIntersection result;
        float tMin = FLT_MAX;

        // Mesh (Triangle)
        for (size_t i = 0, maxTab = meshes.size(); i < maxTab; i++)
        {
            RayTriangleIntersection tmp = meshes[i].intersect(ray);

            if (tmp.intersectionExists)
            {
                if (tmp.t < tMin && tmp.t >= 0)
                {
                    result.intersectionExists = true;
                    result.rayMeshIntersection = tmp;
                    result.typeOfIntersectedObject = TypeObject::mesh;
                    result.objectIndex = i;
                    result.t = tmp.t;
                    tMin = tmp.t;
                }
            }
        }

        // Spheres
        for (size_t i = 0, maxTab = spheres.size(); i < maxTab; i++)
        {
            RaySphereIntersection tmp = spheres[i].intersect(ray);

            if (tmp.intersectionExists)
            {
                if (tmp.t < tMin && tmp.t >= 0)
                {
                    result.intersectionExists = true;
                    result.raySphereIntersection = tmp;
                    result.typeOfIntersectedObject = TypeObject::sphere;
                    result.objectIndex = i;
                    result.t = tmp.t;
                    tMin = tmp.t;
                }
            }
        }

        // Squares
        for (size_t i = 0, maxTab = squares.size(); i < maxTab; i++)
        {
            RaySquareIntersection tmp = squares[i].intersect(ray);
            if (tmp.intersectionExists)
            {
                if (tmp.t < tMin && tmp.t >= 0)
                {
                    result.intersectionExists = true;
                    result.raySquareIntersection = tmp;
                    result.typeOfIntersectedObject = TypeObject::square;
                    result.objectIndex = i;
                    result.t = tmp.t;
                    tMin = tmp.t;
                }
            }
        }

        return result;
    }

    Vec3 rayTraceRecursive(Ray ray, int NRemainingBounces)
    {
        RaySceneIntersection raySceneIntersection = computeIntersection(ray);
        Vec3 color;
        float bias = 1e-2;

        if (NRemainingBounces > 0)
        {
            if (raySceneIntersection.intersectionExists)
            {
                Mesh tmpMesh;
                Vec3 N;
                Vec3 P;

                switch (raySceneIntersection.typeOfIntersectedObject)
                {
                // Mesh
                case TypeObject::mesh:
                    tmpMesh = meshes[raySceneIntersection.objectIndex];
                    N = raySceneIntersection.rayMeshIntersection.normal;
                    P = raySceneIntersection.rayMeshIntersection.intersection;
                    break;

                // Sphere
                case TypeObject::sphere:
                    tmpMesh = spheres[raySceneIntersection.objectIndex];
                    N = raySceneIntersection.raySphereIntersection.normal;
                    P = raySceneIntersection.raySphereIntersection.intersection;
                    break;

                // Square
                case TypeObject::square:
                    tmpMesh = squares[raySceneIntersection.objectIndex];
                    N = raySceneIntersection.raySquareIntersection.normal;
                    P = raySceneIntersection.raySquareIntersection.intersection;
                    break;
                }

                Vec3 refractDirection;
                Vec3 reflectDirection;
                Ray refractRay;
                Ray reflectRay;

                switch (tmpMesh.material.type)
                {
                case MaterialType::Material_Diffuse_Blinn_Phong:
                    color = phongShading(raySceneIntersection, ray, tmpMesh, N, P);
                    break;

                case MaterialType::Material_Glass:
                    // Refraction Ray
                    if (tmpMesh.material.transparency > 0)
                    {
                        refractDirection = refract(ray.direction(), N, tmpMesh.material.index_medium);
                        refractRay = Ray(P + (bias * refractDirection), refractDirection);
                        color = ((1 - tmpMesh.material.transparency) * phongShading(raySceneIntersection, ray, tmpMesh, N, P)) + (tmpMesh.material.transparency * rayTraceRecursive(refractRay, NRemainingBounces - 1));
                    }
                    else
                    {
                        color = phongShading(raySceneIntersection, ray, tmpMesh, N, P);
                    }
                    break;

                case MaterialType::Material_Mirror:
                    // Reflexion Ray
                    if (tmpMesh.material.metalness > 0)
                    {
                        reflectDirection = reflect(ray.direction(), N);
                        reflectRay = Ray(P + (bias * reflectDirection), reflectDirection);
                        color = ((1 - tmpMesh.material.metalness) * phongShading(raySceneIntersection, ray, tmpMesh, N, P)) + (tmpMesh.material.metalness * rayTraceRecursive(reflectRay, NRemainingBounces - 1));
                    }
                    else
                    {
                        color = phongShading(raySceneIntersection, ray, tmpMesh, N, P);
                    }
                    break;

                default:
                    color = Vec3(0.f, 0.f, 0.f);
                    break;
                }
            }
        }
        return color;
    }

    Vec3 reflect(Vec3 v, Vec3 n)
    {
        Vec3 res = v - 2 * Vec3::dot(v, n) * n;
        res.normalize();
        return res;
    }

    Vec3 refract(Vec3 v, Vec3 n, float ior)
    {
        Vec3 res = Vec3(0.f, 0.f, 0.f);
        float dot_nv = Vec3::dot(n, v);
        float k = (ior * ior) * (dot_nv * dot_nv);
        if (k >= 0.f)
        {
            res = (ior * v) + ((ior * (dot_nv + sqrt(k))) * n);
            res.normalize();
        }
        return res;
    }

    Vec3 phongShading(RaySceneIntersection raySceneIntersection, Ray ray, Mesh tmpMesh, Vec3 N, Vec3 P)
    {
        Vec3 color;
        float bias = 1e-3;

        // Vector intersection -> observateur
        Vec3 V = ray.origin() - P;
        V.normalize();

        // Pour chaque lumières (ADD color)
        int maxLight = lights.size();
        for (int i = 0; i < maxLight; i++)
        {
            // Vector intersection -> lumière
            Vec3 L = lights[i].pos - P;
            L.normalize();

            float dotLN = Vec3::dot(L, N);

            // Vector Réflexion Spéculaire
            Vec3 R = 2 * dotLN * N - L;
            R.normalize();

            // Init : diffuse / specular / ambiant
            Vec3 diffuse;
            Vec3 specular;
            Vec3 ambiant = Vec3::compProduct(lights[i].material, tmpMesh.material.ambient_material);

            // Shadow Ray for each point in Light
            vector<Vec3> randPtsInLight = lights[i].randPointSphere();
            Vec3 shadowColor;
            bool hitLight = false;

            // Pour toutes les lumières
            for (size_t i = 0, maxSizeLight = randPtsInLight.size(); i < maxSizeLight && !hitLight; i++)
            {
                // Vector intersection -> lumière
                Vec3 Ltmp = randPtsInLight[i] - P;
                Ltmp.normalize();

                Ray shadowRay = Ray(P + (bias * Ltmp), Ltmp);
                RaySceneIntersection shadowIntersect = computeIntersection(shadowRay);

                if (shadowIntersect.t > Ltmp.length())
                {
                    hitLight = true;
                }
            }

            if (dotLN >= 0)
            {
                if (hitLight)
                {
                    diffuse = Vec3::compProduct(lights[i].material, tmpMesh.material.diffuse_material) * dotLN;
                    specular = Vec3::compProduct(lights[i].material, tmpMesh.material.specular_material) * pow(Vec3::dot(R, V), tmpMesh.material.shininess);
                }
            }
            color += ambiant + diffuse + specular;
        }
        return color;
    }

    Vec3 rayTrace(Ray const &rayStart)
    {
        Vec3 color = rayTraceRecursive(rayStart, PARAM_NUMBER_BOUNCES);
        return color;
    }

    void setup_single_sphere()
    {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(-5, 5, 5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }
        {
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(0., 0., 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(1., 1., 1);
            s.material.specular_material = Vec3(0.2, 0.2, 0.2);
            s.material.shininess = 20;
        }
    }

    void setup_single_square()
    {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(-5, 5, 5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }

        {
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.build_arrays();
            s.material.diffuse_material = Vec3(0.8, 0.8, 0.8);
            s.material.specular_material = Vec3(0.8, 0.8, 0.8);
            s.material.shininess = 20;
        }
    }

    void setup_double_sphere()
    {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(-5, 5, 5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }
        {
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(0., 0., 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(1., 1., 1);
            s.material.specular_material = Vec3(0.2, 0.2, 0.2);
            s.material.shininess = 20;
        }
        {
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(2., 0., -2.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(1., 1., 1);
            s.material.specular_material = Vec3(0.2, 0.2, 0.2);
            s.material.shininess = 20;
        }
    }

    void setup_cornell_box()
    {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(0.0, 1.5, 0.0);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }
        { // Back Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
        }
        { // Left Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 0., 0.);
            s.material.specular_material = Vec3(1., 0., 0.);
            s.material.shininess = 16;
        }
        { // Right Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(0.0, 1.0, 0.0);
            s.material.specular_material = Vec3(0.0, 1.0, 0.0);
            s.material.shininess = 16;
        }
        { // Floor
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
            s.material.shininess = 16;
        }
        { // Ceiling
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
            s.material.shininess = 16;
        }
        { // GLASS Sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3(1., 0., 0.);
            s.material.specular_material = Vec3(1., 0., 0.);
            s.material.shininess = 16;
            s.material.transparency = 0.75;
            s.material.index_medium = 1.5;
        }
        { // MIRRORED Sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
            s.material.metalness = 1.0;
        }
    }

    void setup_giant_cornell_box()
    {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(0.0, 1.5, 0.0);
            light.radius = 0.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }

        { // Back Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., -1., -2.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 4., 4.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
            s.material.metalness = 1.0;
        }
        { // Left Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-5., -1., -2.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 7., 4.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 0., 0.);
            s.material.specular_material = Vec3(1., 0., 0.);
            s.material.shininess = 16;
        }
        { // Right Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., -1., -2.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 7., 4.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(0.0, 1.0, 0.0);
            s.material.specular_material = Vec3(0.0, 1.0, 0.0);
            s.material.shininess = 16;
        }
        { // Floor
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., -5., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 4., 7.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
            s.material.shininess = 16;
        }
        { // Ceiling
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., -2., -4.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 4., 7.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
            s.material.shininess = 16;
        }
        { // Front Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., -1., -8.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 4., 4.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
            s.material.metalness = 1.0;
        }
        { // GLASS Sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3(1., 0., 0.);
            s.material.specular_material = Vec3(1., 0., 0.);
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }
        { // MIRRORED Sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
            s.material.metalness = 1.0;
        }
    }

    void setup_double_light()
    {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        { // Light Left
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(1.0, 1.5, 0.0);
            light.radius = 0.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 0, 1);
            light.isInCamSpace = false;
        }

        { // Light Right
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(-1.0, 1.5, 0.0);
            light.radius = 0.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(0, 1, 1);
            light.isInCamSpace = false;
        }

        { // Back Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
        }

        { // Left Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(0., 1., 1.);
            s.material.specular_material = Vec3(0., 1., 1.);
            s.material.shininess = 16;
        }

        { // Right Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1.0, 0.0, 1.0);
            s.material.specular_material = Vec3(1.0, 0.0, 1.0);
            s.material.shininess = 16;
        }

        { // Floor
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
            s.material.shininess = 16;
        }

        { // Ceiling
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
            s.material.shininess = 16;
        }
        { // GLASS Sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3(1., 0., 0.);
            s.material.specular_material = Vec3(1., 0., 0.);
            s.material.shininess = 16;
            s.material.transparency = 0.5;
            s.material.index_medium = 1.0;
        }
        { // MIRRORED Sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
            s.material.metalness = 1.0;
        }
    }

    void setup_mesh()
    {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(0.0, 1.5, 0.0);
            light.radius = 0.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }

        { // Back Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
        }

        { // Left Wall

            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 0., 0.);
            s.material.specular_material = Vec3(1., 0., 0.);
            s.material.shininess = 16;
        }

        { // Right Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(0.0, 1.0, 0.0);
            s.material.specular_material = Vec3(0.0, 1.0, 0.0);
            s.material.shininess = 16;
        }

        { // Floor
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
            s.material.shininess = 16;
        }

        { // Ceiling
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
            s.material.shininess = 16;
        }

        { // Avion
            meshes.resize(meshes.size() + 1);
            Mesh &m = meshes[meshes.size() - 1];
            m.loadOFF("./Maillages/avion_64.off");
            m.centerAndScaleToUnit();

            m.rotate_x(-45);
            m.rotate_y(45);

            m.recomputeNormals();
            m.build_arrays();

            m.material.type = Material_Mirror;
            m.material.diffuse_material = Vec3(1., 0., 0.);
            m.material.specular_material = Vec3(1., 0., 0.);
            m.material.shininess = 16;
            m.material.transparency = 1.0;
            m.material.index_medium = 1.4;
        }
    }
};

#endif
