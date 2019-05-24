/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Roughness.h
 * Author: ytsong
 *
 * Created on April 11, 2019, 12:31 PM
 */

#ifndef ROUGHNESS_H
#define ROUGHNESS_H

#include "ia/mesh.h"
#include "ia/vertex.h"
#include "ia/triangle.h"

class Roughness {
public:
    Roughness(Spatial_Mesh &mesh){      roughness.assign(mesh.get_vertices_num(),0);time_for_compute=0;};

    void compute_values(Spatial_Mesh &mesh);
    void print_stats(Spatial_Mesh &mesh);
    void store_result(Spatial_Mesh &mesh);
private:
     coord_type compute(itype vid, Spatial_Mesh &mesh);
     dvect roughness;

         double time_for_compute;
};

#endif /* ROUGHNESS_H */

