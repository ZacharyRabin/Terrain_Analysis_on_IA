#include "formangradientvector.h"
#include <stdio.h>
#include "assert.h"

#include <climits>

FormanGradientVector::FormanGradientVector(Spatial_Mesh* mesh)
{
    this->mesh = mesh;

    forman_gradient = vector<unsigned int>(mesh->get_triangles_num(), NON_VALID_CONF);
    unsigned int tetCaseCountOuterAllowed=0;
    unsigned int tot_case=0;

    CompressedToExpandedLUT compressed_to_expanded_local = CompressedToExpandedLUT(512);
    ExpandedToCompressedLUT expanded_to_compressed_local = ExpandedToCompressedLUT();
    compressed_to_expanded = CompressedToExpandedLUT(512);
    expanded_to_compressed = ExpandedToCompressedLUT();

    for(unsigned int i = 0; i < 511; i++)
    {

        Arrows ga = static_cast<Arrows>( static_cast<unsigned int>( i ));
        if( isValidTetCase(ga) )
        {
            compressed_to_expanded_local[tetCaseCountOuterAllowed]       = ga;
            expanded_to_compressed_local[ga]                             = tetCaseCountOuterAllowed;
            tetCaseCountOuterAllowed++;
        }

        tot_case++;
    }

    compressed_to_expanded = compressed_to_expanded_local;
    expanded_to_compressed = expanded_to_compressed_local;

}

bool FormanGradientVector::isValidTetCase(Arrows arrow)
{
    // This disables the outer four bits (corresponding to the FT relations for the opposite T to each F) from being valid

    // test vertices
    if( BitTwiddle::popcount( arrow & CONTAINS_V0 ) > 1 )		return false;
    if( BitTwiddle::popcount( arrow & CONTAINS_V1 ) > 1 )		return false;
    if( BitTwiddle::popcount( arrow & CONTAINS_V2 ) > 1 )		return false;

    // test edges
    if( BitTwiddle::popcount( arrow & CONTAINS_E01 ) > 1 )		return false;
    if( BitTwiddle::popcount( arrow & CONTAINS_E02 ) > 1 )		return false;
    if( BitTwiddle::popcount( arrow & CONTAINS_E12 ) > 1 )		return false;

    // test faces
    if( BitTwiddle::popcount( arrow & CONTAINS_F012 ) > 1 )		return false;

    // valid!
    return true;
}


TriGradient FormanGradientVector::convert_compressed_to_expand(unsigned int ga){
    if(ga == NON_VALID_CONF) return TriGradient();
    return TriGradient(compressed_to_expanded[ga]);
}

unsigned int FormanGradientVector::convert_expand_to_compressed(Arrows ga){
    return expanded_to_compressed[ga];
}

void FormanGradientVector::setVE(int v, int v2){

    Edge edg=Edge(v,v2);
    vector<int> et = mesh->ET(edg);

    mesh->get_vertex(v).set_VTstar(et[0]);
    for(unsigned int i=0; i<et.size(); i++){
        TriGradient ga = convert_compressed_to_expand(forman_gradient[et[i]]);
        ga.setVE(mesh->get_triangle(et[i]).vertex_index(v), mesh->get_triangle(et[i]).vertex_index(v2));
        forman_gradient[et[i]] = convert_expand_to_compressed(ga.getArrow());
    }

}

void FormanGradientVector::setEF(int v, int f){

    TriGradient ga = convert_compressed_to_expand(forman_gradient[f]);
    ga.setEF(mesh->get_triangle(f).vertex_index(v));
    forman_gradient[f] = convert_expand_to_compressed(ga.getArrow());
}

Edge* FormanGradientVector::getVE(int v){

    int t = mesh->get_vertex(v).get_VTstar();
    int v2 = convert_compressed_to_expand(forman_gradient[t]).get_vertex_pair(mesh->get_triangle(t).vertex_index(v));

    if(v2 != -1){
        return new Edge(v, mesh->get_triangle(t).TV(v2));
    }

    return NULL;
}

int FormanGradientVector::getEF(Edge* e){

    vector<int> et = mesh->ET(*e);
    for(int i=0; i<et.size(); i++){
        int t = convert_compressed_to_expand(forman_gradient[et[i]]).get_edge_pair(mesh->get_triangle(et[i]).vertex_index(e->EV(0)),mesh->get_triangle(et[i]).vertex_index(e->EV(1)));
        if(t == 3) return et[i];
    }

    return -1;
}

int FormanGradientVector::getFE(int tri){

    return convert_compressed_to_expand(forman_gradient[tri]).get_face_pair();
}


void FormanGradientVector::freeVE(int v1, int v2){
    Edge edge=Edge(v1,v2);
    vector<int> et = mesh->ET(edge);
    for(unsigned int i=0; i<et.size(); i++){
        TriGradient ga = convert_compressed_to_expand(forman_gradient[et[i]]);
        ga.clearVE(mesh->get_triangle(et[i]).vertex_index(v1), mesh->get_triangle(et[i]).vertex_index(v2));
        forman_gradient[et[i]] = convert_expand_to_compressed(ga.getArrow());
    }
}

void FormanGradientVector::freeEF(int v, int f){
    TriGradient ga = convert_compressed_to_expand(forman_gradient[f]);
    ga.clearEF(mesh->get_triangle(f).vertex_index(v));
    forman_gradient[f] = convert_expand_to_compressed(ga.getArrow());
}

bool FormanGradientVector::is_vertex_critical(int v){
    return convert_compressed_to_expand(forman_gradient[mesh->get_vertex(v).get_VTstar()]).is_vertex_unpaired(mesh->get_triangle(mesh->get_vertex(v).get_VTstar()).vertex_index(v));
}

bool FormanGradientVector::is_edge_critical(int v1, int v2){

    int t1 = mesh->get_vertex(v1).get_VTstar();
    int t2 = mesh->get_vertex(v2).get_VTstar();

    if((mesh->get_triangle(t1).has_vertex(v2) && !(convert_compressed_to_expand(forman_gradient[t1]).is_edge_unpaired(mesh->get_triangle(t1).vertex_index(v1),mesh->get_triangle(t1).vertex_index(v2)))) ||
       (mesh->get_triangle(t2).has_vertex(v1) && !(convert_compressed_to_expand(forman_gradient[t2]).is_edge_unpaired(mesh->get_triangle(t2).vertex_index(v2),mesh->get_triangle(t2).vertex_index(v1))))) {
        return false;
    }
    Edge edge=Edge(v1,v2);
    vector<int> et = mesh->ET(edge);
    for(unsigned int i=0; i<et.size(); i++){
        if(!(convert_compressed_to_expand(forman_gradient[et[i]]).is_edge_unpaired(mesh->get_triangle(et[i]).vertex_index(v1),mesh->get_triangle(et[i]).vertex_index(v2))))
            return false;
    }

    return true;
}

bool FormanGradientVector::is_face_critical(int f){
    return convert_compressed_to_expand(forman_gradient[f]).is_face_unpaired();
}






