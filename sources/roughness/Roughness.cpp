/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Roughness.cpp
 * Author: ytsong
 * 
 * Created on April 11, 2019, 12:31 PM
 */

#include "Roughness.h"
#include "utilities/timer.h"

void Roughness::compute_values(Spatial_Mesh& mesh)
{

    
        for(itype i=0; i<mesh.get_vertices_num(); i++)
    {
        roughness[i]=compute(i,mesh); 
    }
}


coord_type Roughness::compute(itype vid, Spatial_Mesh &mesh){


    coord_type zDistSum = 0.0;
    coord_type zSum=0.0;
    coord_type rough=0;

    ivect vv = mesh.VV(vid);
               
    //t==-1 cannot happen if the mesh has no isolated vertices
    if (vv.size() == 0) return 0.0;

    for(auto v1_id:vv)
    {

        Vertex& v1= mesh.get_vertex(v1_id);
        zSum+=v1.get_c(2);
    
    }
    coord_type zAVG=zSum/vv.size();
    for(auto v1_id:vv)
    {
        Vertex& v1= mesh.get_vertex(v1_id);
        zDistSum+=(v1.get_c(2)-zAVG)*(v1.get_c(2)-zAVG);
    
    }
    rough=sqrt(zDistSum/vv.size());
 

    return rough;

}


void Roughness::print_stats(Spatial_Mesh& mesh){
    
    coord_type min=INFINITY, max=-INFINITY, avg=0;
    for(itype v=0; v<mesh.get_vertices_num(); v++)
    {
        if(roughness[v] < min)
            min = roughness[v];
        if(roughness[v] > max)
            max = roughness[v];
        avg += roughness[v];
    }
    cerr<<"[STATS] roughness min: "<<min<<" avg: "<<avg/(coord_type)mesh.get_vertices_num()<<" max: "<<max<<endl;
    
}

void Roughness::store_result(Spatial_Mesh& mesh){

    for(itype v=0;v<mesh.get_vertices_num();v++)
    {
        mesh.get_vertex(v).add_field(roughness[v]);
       
    }

}