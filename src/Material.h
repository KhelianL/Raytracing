#ifndef MATERIAL_H
#define MATERIAL_H

#include "imageLoader.h"
#include "Vec3.h"
#include <cmath>

#include <GL/glut.h>

enum MaterialType
{
    Material_Diffuse_Blinn_Phong,
    Material_Glass,
    Material_Mirror
};

struct Material
{
    MaterialType type;
    Vec3 ambient_material;
    Vec3 diffuse_material;
    Vec3 specular_material;

    double shininess;
    float index_medium;     // IOR
    float transparency;     // [0 - 1]
    float metalness;        // [0 - 1]

    Material()
    {
        type = Material_Diffuse_Blinn_Phong;
        ambient_material = Vec3(0., 0., 0.);
        diffuse_material = Vec3(0., 0., 0.);
        specular_material = Vec3(0., 0., 0.);

        shininess = 1.0;

        index_medium = 1.0;
        transparency = 0.0;
        metalness = 0.0;
    }
};

#endif // MATERIAL_H
