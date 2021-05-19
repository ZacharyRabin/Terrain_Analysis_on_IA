/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Gradient.h
 * Author: ytsong
 *
 * Created on April 11, 2019, 2:59 PM
 */

#ifndef GRADIENT_H
#define GRADIENT_H

#include "ia/mesh.h"
#include "ia/vertex.h"
#include "ia/triangle.h"

#ifdef _WIN32   // Windows system specific
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#else          // Unix based system specific
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>
#include <eigen3/Eigen/Eigenvalues>
#endif



#include <complex>
#include <limits>
#include "utilities/timer.h"

using namespace Eigen;



class Gradient {
public:
    Gradient(){

    //    multifield.assign(mesh.get_vertices_num(),0);
    
    };
    Gradient(string mode){
       if(mode=="All"){
        fields.push_back(0);
        fields.push_back(1);
        fields.push_back(2);
        fields.push_back(3);
        fields.push_back(4);
        block_time=0;
       }
       else if (mode=="rRGB")
       {
        fields.push_back(1);
        fields.push_back(2);
        fields.push_back(3);
        fields.push_back(4);
       
       }
       else if (mode=="rz"){
           fields.push_back(0);
           fields.push_back(1);
       }
    };
    
    //Here we don't need to compute the gradient of all the triangles in the mesh. So we don't need PCE_gradient.
    
    void vertex_gradient(Spatial_Mesh& mesh, int field_index);
    
    void multi_field(Spatial_Mesh& mesh);
    void print_time(){
        cout<<"The time of this block is:"<<block_time<<endl;
        
    }
    
    
     inline void compute_field_stats(Spatial_Mesh &mesh,int factor)
    {
        for(auto f:fields){
        coord_type min=INFINITY, max=-INFINITY;
        
        for(itype v=0; v<mesh.get_vertices_num(); v++)
        {
            coord_type c = mesh.get_vertex(v).get_field(f);

            if(c < min)
                min = c;
            if(c > max)
                max = c;
          //  avg += c;
        }
           coord_type range=(max-min)/factor;
        this->ranges.push_back(range);
        }
    }
    
    void print_stats(Spatial_Mesh& mesh){
    
    coord_type min=INFINITY, max=-INFINITY, avg=0;
    itype fid= mesh.get_vertex(10).get_fields_num()-1;
    
    for(itype v=0; v<mesh.get_vertices_num(); v++)
    {
        if(mesh.get_vertex(v).get_field(fid) < min)
            min = mesh.get_vertex(v).get_field(fid);
        if(mesh.get_vertex(v).get_field(fid) > max)
            max = mesh.get_vertex(v).get_field(fid);
        avg += mesh.get_vertex(v).get_field(fid);
    }
    cerr<<"[STATS] multifield min: "<<min<<" avg: "<<avg/(coord_type)mesh.get_vertices_num()<<" max: "<<max<<" number of field:"<<fid+1<<endl;
    
}

private:

    dvect PCE_compute(itype tid, Spatial_Mesh &mesh, int field_index);
    FG vertex_compute(itype vid,Spatial_Mesh &mesh, int field_index);
    void multifield_compute(itype vid, Spatial_Mesh &mesh);
    coord_type block_time;
    
    
    coord_type cross_2d(dvect i,dvect j);
    coord_type dot_2d(dvect i,dvect j);
    

    
    dvect ranges;
    ivect fields;
    
};
#endif /* GRADIENT_H */

