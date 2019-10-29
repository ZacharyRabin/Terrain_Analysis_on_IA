#include "formangradientvector.h"

vector<int> FormanGradientVector::descending_2cells_extraction(bool with_geometry){

    vector<int> labeling = vector<int>(mesh->get_triangles_num(),-1);
    for(int i=0; i<mesh->get_triangles_num(); i++){

        if(is_face_critical(i)){

            queue<Edge* > coda;
            labeling[i]=i;

            vector<int> simplexi;
            int v=-1;
            for(int j=0; j<3; j++){
                    simplexi.push_back(mesh->get_triangle(i).TV(j));

            }

            for(int j=0; j<3; j++){
                vector<int> simplex = simplexi;
                simplex.erase(simplex.begin()+j);
                coda.push(new Edge(simplex[0], simplex[1]));
            }
            simplexi.clear();

            while(!coda.empty()){
                Edge* edge = coda.front();
                coda.pop();

                int tri = getEF(edge);
                if(tri != -1){

                    labeling[tri]=i;

                    for(int j=0; j<3; j++){
                        Edge* new_e = mesh->get_triangle(tri).TE_p(j);
                        if(*new_e != *edge){
                            coda.push(new_e);
                        }
                        else delete new_e;
                    }

                    delete edge;
                }
            }
        }
    }


    if(with_geometry)
        writeVTK_2cells("descending2cells.vtk", labeling);

    return labeling;
}


void FormanGradientVector::descending_1cells_extraction(bool with_geometry){

    queue<int> coda;

    set<int> vertici;
    set<pair<int, bool> > critici;
    set<pair<int,int> > edges;

    vector<bool> visited = vector<bool>(mesh->get_triangles_num(),false);
    for(int i=0; i<mesh->get_triangles_num(); i++){
        visited[i]=true;
        for(int j=0; j<3; j++){
            if(mesh->get_triangle(i).TT(j) == -1 || !visited[mesh->get_triangle(i).TT(j)]){

                Edge* edge = mesh->get_triangle(i).TE_p(j);
                if(is_edge_critical(edge->EV(0), edge->EV(1))){


                    if(with_geometry){
                        mesh->get_vertex(edge->EV(0)).get_c(2) > mesh->get_vertex(edge->EV(1)).get_c(2) ? edges.insert(pair<int,int>(edge->EV(0),edge->EV(1))) : edges.insert(pair<int,int>(edge->EV(1),edge->EV(0)));
                        vertici.insert(edge->EV(0));
                        vertici.insert(edge->EV(1));
                    }

                    coda.push(edge->EV(0));
                    coda.push(edge->EV(1));
                    delete edge;

                    while(!coda.empty()){

                        int vert = coda.front();
                        coda.pop();

                        edge = getVE(vert);
                        if(edge != NULL){
                            if(with_geometry){
                                mesh->get_vertex(edge->EV(0)).get_c(2) > mesh->get_vertex(edge->EV(1)).get_c(2) ? edges.insert(pair<int,int>(edge->EV(0),edge->EV(1))) : edges.insert(pair<int,int>(edge->EV(1),edge->EV(0)));
                                vertici.insert(edge->EV(0));
                                vertici.insert(edge->EV(1));
                            }
                            int v2 = edge->EV(0) == vert ? edge->EV(1) : edge->EV(0);
                            delete edge;
                            coda.push(v2);
                        }
                    }
                }
                else delete edge;
            }
        }
    }


    if(with_geometry)
        writeVTK_1cells("descending1cells.vtk", vertici, edges);
}

vector<int> FormanGradientVector::ascending_2cells_extraction(bool with_geometry){

    queue<int> coda;
    vector<int> label = vector<int>(mesh->get_vertices_num(), -1);

    int critical_edge =0;
    int minima=0;

    for(int i=0; i<mesh->get_vertices_num(); i++){
        if(is_vertex_critical(i)){
            minima++;
            label[i]=i;

            vector<Edge*> ve = mesh->VE_p(i);
            for(unsigned int j=0; j<ve.size(); j++){
                int v2 = ve[j]->EV(0) == i ? ve[j]->EV(1) : ve[j]->EV(0);
                Edge* edge = getVE(v2);
                if(edge != NULL && *edge == *(ve[j])){
                    coda.push(v2);
                }

                if(edge != NULL) delete edge;
                delete ve[j];
            }


            while(!coda.empty()){

                int v = coda.front();
                coda.pop();
                label[v]=i;

                ve = mesh->VE_p(v);
                for(unsigned int j=0; j<ve.size(); j++){
                    int v2 = ve[j]->EV(0) == v ? ve[j]->EV(1) : ve[j]->EV(0);
                    Edge* edge = getVE(v2);
                    if(edge != NULL && *edge == *(ve[j])){
                        coda.push(v2);
                    }
                    if(edge != NULL) delete edge;
                    delete ve[j];
                }
            }

        }
    }


    if(with_geometry)
        writeVTK_2cells_on_vert("ascending2cells.vtk", label);

    return label;
}

void FormanGradientVector::ascending_1cells_extraction(bool with_geometry){

    vector<bool> visited = vector<bool>(mesh->get_triangles_num(),false);
    map<int,int> visited_triangle;
    set<pair<Vertex,Vertex> > edges;

    for(int i=0; i<mesh->get_triangles_num(); i++){
        visited[i]=true;
        for(int j=0; j<3; j++){
            if(mesh->get_triangle(i).TT(j) == -1 || !visited[mesh->get_triangle(i).TT(j)]){

                queue<int> coda;
                Edge* edge = mesh->get_triangle(i).TE_p(j);
                if(is_edge_critical(edge->EV(0), edge->EV(1)))
                {

                    vector<int> et = mesh->ET(*edge);
                    for(int j=0; j<et.size(); j++){
                        coda.push(et[j]);

                        Vertex v1 = mesh->get_vertex(mesh->get_triangle(et[j]).TV(0));
                        v1 += mesh->get_vertex(mesh->get_triangle(et[j]).TV(1));
                        v1 += mesh->get_vertex(mesh->get_triangle(et[j]).TV(2));
                        v1 /= 3.0;

                        Vertex ve = mesh->get_vertex(edge->EV(0));
                        ve += mesh->get_vertex(edge->EV(1));
                        ve /= 2.0;
                        if(with_geometry){
                            edges.insert(pair<Vertex,Vertex>(v1,ve));
                        }
                    }

                    while(!coda.empty()){

                        int t = coda.front();
                        coda.pop();

                        visited_triangle[t]=i;

                        if(is_face_critical(t)){
                            //here critical face
                        }
                        else{
                            int e = getFE(t);

                            int t_adj = mesh->get_triangle(t).TT(e);
                            if(t_adj != -1){
                                coda.push(t_adj);

                                if(with_geometry){
                                    Vertex v1 = mesh->get_vertex(mesh->get_triangle(t_adj).TV(0));
                                    v1 += mesh->get_vertex(mesh->get_triangle(t_adj).TV(1));
                                    v1 += mesh->get_vertex(mesh->get_triangle(t_adj).TV(2));
                                    v1 /= 3.0;

                                    Edge* edg = mesh->get_triangle(t).TE_p(e);
                                    Vertex ve = mesh->get_vertex(edg->EV(0));
                                    ve += mesh->get_vertex(edg->EV(1));
                                    ve /= 2.0;

                                    Vertex v2 = mesh->get_vertex(mesh->get_triangle(t).TV(0));
                                    v2 += mesh->get_vertex(mesh->get_triangle(t).TV(1));
                                    v2 += mesh->get_vertex(mesh->get_triangle(t).TV(2));
                                    v2 /= 3.0;

                                    edges.insert(pair<Vertex,Vertex>(v1,ve));
                                    edges.insert(pair<Vertex,Vertex>(v2,ve));
                                    }
                            }
                        }
                    }
                }
                else delete edge;
            }
        }
    }

    if(with_geometry)
        writeVTK_1cells("ascending1cells.vtk",edges);
}



vector<int> FormanGradientVector::count_critical_simplexes(){
    vector<int> cp(3,0);

    for(int i=0; i<mesh->get_triangles_num(); i++){
        if(is_face_critical(i))
            cp[2]++;
    }

    for(int i=0; i<mesh->get_vertices_num(); i++){
        if(is_vertex_critical(i)){
            cp[0]++;
        }

        vector<int> vv = mesh->VV(i);
        for(auto v : vv){
            if(v < i){
                if(is_edge_critical(v,i))
                    cp[1]++;
            }
        }
    }

    return cp;
}
