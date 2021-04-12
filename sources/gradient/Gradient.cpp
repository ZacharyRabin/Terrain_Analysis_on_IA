/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   Gradient.cpp
 * Author: ytsong
 *
 * Created on April 11, 2019, 2:59 PM
 */

#include "Gradient.h"


dvect Gradient::PCE_compute(itype tid, Spatial_Mesh &mesh, int field_index){

    coord_type area;

    Triangle &t = mesh.get_triangle(tid);



    Vertex &vi = mesh.get_vertex(t.TV(0));
    Vertex &vj = mesh.get_vertex(t.TV(1));
    Vertex &vk = mesh.get_vertex(t.TV(2));

    dvect ki={vi.get_c(0)-vk.get_c(0),vi.get_c(1)-vk.get_c(1)};
    dvect ij={vj.get_c(0)-vi.get_c(0),vj.get_c(1)-vi.get_c(1)};
        // compute the vectors of the edges after rotate by 90 degrees
    dvect ki_vert= { vk.get_c(1)-vi.get_c(1) , vi.get_c(0)-vk.get_c(0)};
    dvect ij_vert= {vi.get_c(1)-vj.get_c(1) , vj.get_c(0)-vi.get_c(0)};
    // compute the area of the triangle
    area=abs(0.5*(ki[0]*ij[1]-ki[1]*ij[0]));
    int field_pos=fields[field_index];

    dvect area_multiply_ki={(vj.get_field(field_pos)-vi.get_field(field_pos))*ki_vert[0]/(2*area),(vj.get_field(field_pos)-vi.get_field(field_pos))*ki_vert[1]/(2*area)};
    dvect area_multiply_ij={(vk.get_field(field_pos)-vi.get_field(field_pos))*ij_vert[0]/(2*area),(vk.get_field(field_pos)-vi.get_field(field_pos))*ij_vert[1]/(2*area)};

    coord_type range=ranges[field_index];
    dvect g={(area_multiply_ki[0]+area_multiply_ij[0])/range,(area_multiply_ki[1]+area_multiply_ij[1])/range};
    return g;

}


void Gradient::multi_field(Spatial_Mesh& mesh){
 for(itype i=0; i<mesh.get_vertices_num(); i++)
    {
        //coord_type multifield;
        /*multifield=*/multifield_compute(i,mesh);

    }

}


FG Gradient::vertex_compute(itype vid, Spatial_Mesh& mesh, int field_index){

        bool is_border=false;

        ivect vt = mesh.VT(vid,is_border);
    
        coord_type sum_area=0.0;
        dvect sum_gradient={0,0};
        dvect gradient_v={0,0};

        Vertex &vp=mesh.get_vertex(vid);

        for(auto tid:vt){
            dvect gradient_t;
            Triangle &t=mesh.get_triangle(tid);
            Vertex &vq=(t.TV(0)!=vid)?mesh.get_vertex(t.TV(0)):mesh.get_vertex(t.TV(1));
            Vertex &vr=((t.TV(1)!=vid)&&(t.TV(0)!=vid))?mesh.get_vertex(t.TV(1)):mesh.get_vertex(t.TV(2));

            gradient_t=this->PCE_compute(tid,mesh,field_index);

            dvect qp={vq.get_c(0)-vp.get_c(0),vq.get_c(1)-vp.get_c(1)};
            dvect qr={vq.get_c(0)-vr.get_c(0),vq.get_c(1)-vr.get_c(1)};
            dvect rp={vr.get_c(0)-vp.get_c(0),vr.get_c(1)-vp.get_c(1)};
            dvect rq={-qr[0],-qr[1]};

             coord_type Voronoi_area=0;
            if(dot_2d(qp,rp)<=0)
            { Voronoi_area=abs(0.25*(qp[0]*rp[1]-qp[1]*rp[0]));

            }
            else if((dot_2d(rp,rq)<0)||(dot_2d(qp,qr)<0))
            {Voronoi_area=abs(0.125*(qp[0]*rp[1]-qp[1]*rp[0]));
          //  cout<<"obtuse at Q: "<<dot_2d(qp,qr)<<"obtuse at R: "<<dot_2d(rp,rq)<<endl;
            }
            else
            {Voronoi_area=0.125*(abs((rp[0]*rp[0]+rp[1]*rp[1])*dot_2d(qp,qr)/(cross_2d(qp,qr)))+abs((qp[0]*qp[0]+qp[1]*qp[1])*dot_2d(rp,rq)/cross_2d(rp,rq)));
           // cout<<"no obtuse"<<endl;
            }
         //  cout<<"Voronoi Area: "<<Voronoi_area<<endl;
            sum_area+=Voronoi_area;

            sum_gradient[0]+=gradient_t[0]*Voronoi_area;
            sum_gradient[1]+=gradient_t[1]*Voronoi_area;

        }
        gradient_v[0]=sum_gradient[0]/sum_area;
        gradient_v[1]=sum_gradient[1]/sum_area;



        FG gradient;
        gradient.push_back(gradient_v[0]);
        gradient.push_back(gradient_v[1]);
        return gradient;
}


void Gradient::multifield_compute(itype vid,Spatial_Mesh &mesh){

      vect_FG FG_matrix;


      for(int i=0;i<fields.size();i++){
          FG_matrix.push_back(vertex_compute(vid,mesh,i));
      }



      itype field_num=FG_matrix.size();
        MatrixX2d gradientMatrix(field_num,2);

        for(int i=0;i<field_num;i++)
            for(int j=0;j<2;j++)
            {
                gradientMatrix(i,j)=FG_matrix[i][j];

            }

        MatrixXd input=gradientMatrix.transpose()*gradientMatrix;
        EigenSolver<MatrixXd> solver;
        solver.compute(input,false);

        double max = numeric_limits<double>::min();
	for(int i=0; i<solver.eigenvalues().rows(); i++)
	{
		complex<double> c = solver.eigenvalues().coeff(i,0);
        	// cout << abs(c) << endl;
		if(max < abs(c))
			max = abs(c);
	}
    //    cout<<"the result is:"<<sqrt(max)<<endl;
        coord_type multifield=0;
        multifield=sqrt(max);
        mesh.get_vertex(vid).add_field(multifield);
}

coord_type Gradient::dot_2d(dvect i,dvect j){
        return (i[0]*j[0]+i[1]*j[1]);
    }

coord_type Gradient::cross_2d(dvect i,dvect j){
        return (i[0]*j[1]-i[1]*j[0]);
    }
