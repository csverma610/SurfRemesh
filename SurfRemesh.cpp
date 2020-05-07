#include "SurfRemesh.h"

#include <algorithm>
#include "trilib.hpp"
#include <cassert>
#include <sstream>
#include <iostream>
#include <fstream>
#include <limits>

using namespace std;
using namespace JMath;

#define ANSI_DECLARATORS
extern "C" {
#include <triangle.h>
}

SurfRemesh :: SurfRemesh()
{
    maxEdgeLength = std::numeric_limits<double>::max();
}
////////////////////////////////////////////////////////////////////////////////
void SurfRemesh :: readOFFMesh( const string &filename)
{
    ifstream ifile(filename.c_str(), ios::in);
    if( ifile.fail() ) {
        cout << "Error: Can't open input file " << endl;
        return;
    }

    string str;
    ifile >> str;
    if( str != "OFF") {
        cout << "Error: input file not in the off format" << endl;
        return;
    }

    int numPoints, numFaces, numEdges;
    ifile >> numPoints >> numFaces >> numEdges;

    if( numPoints < 1) {
        cout << "Warning: Input file has no points " << endl;
        return;
    }

    double x, y, z;
    mesh.nodes.resize(numPoints);

    for( size_t i = 0; i < numPoints; i++) {
        ifile >> x >> y >> z;
        auto v = std::make_shared<Node>();
        v->id  = i;
        v->xyz = {x,y,z};
        mesh.nodes[i] = v;
    }

    mesh.faces.resize(numFaces);

    int  nn, i0, i1, i2;
    for( size_t i = 0; i < numFaces; i++) {
        ifile >> nn; assert(nn == 3);
        ifile >> i0 >> i1 >> i2;
        auto n0 = mesh.nodes[i0];
        auto n1 = mesh.nodes[i1];
        auto n2 = mesh.nodes[i2];
        auto face = std::make_shared<Face>();
        face->nodes = {n0, n1, n2};
        mesh.faces[i] = face;
    }
}

////////////////////////////////////////////////////////////////////////////////

int SurfRemesh::Edge::refine( double required_length)
{
    newNodes.clear();

    double currlen = length( nodes[0]->xyz, nodes[1]->xyz);

    if( currlen < required_length) return 0;

    int nsegments = max(2.0, currlen/required_length);
    int npoints   = nsegments-1;

    newNodes.resize(npoints);

    double dt = 1.0/(double)nsegments;
    for( int i = 0; i < npoints; i++) {
        double t = (i+1)*dt;
        double x = (1-t)*nodes[0]->xyz[0] + t*nodes[1]->xyz[0];
        double y = (1-t)*nodes[0]->xyz[1] + t*nodes[1]->xyz[1];
        double z = (1-t)*nodes[0]->xyz[2] + t*nodes[1]->xyz[2];
        newNodes[i] = std::make_shared<Node>();
        newNodes[i]->xyz = {x,y,z};
    }
    return 1;
}

////////////////////////////////////////////////////////////////////////////////

SurfRemesh::EdgePtr SurfRemesh::Face::getEdgeAt( int i)
{
    auto v0 = nodes[i];
    auto v1 = nodes[(i+1)%3];
    auto vm = min(v0,v1);
    for( auto e : vm->edges) {
        if( e->isSame(v0,v1) ) return e;
    }

    assert(v0 != v1);
    auto newedge = std::make_shared<Edge>();
    newedge->nodes[0] = v0;
    newedge->nodes[1] = v1;
    vm->edges.push_back(newedge);
    return newedge;
}

////////////////////////////////////////////////////////////////////////////////

vector<double> SurfRemesh::Face::project()
{
    uvCorners[0] = {0.0, 0.0};

    double len0  = length(nodes[0]->xyz, nodes[1]->xyz);
    uvCorners[1] = {len0, 0.0};

    double angle  = angleAt(nodes[0]->xyz, nodes[1]->xyz, nodes[2]->xyz, ANGLE_IN_RADIANS);
    double len1   = length(nodes[0]->xyz, nodes[2]->xyz);
    uvCorners[2]  = {len1*cos(angle), len1*sin(angle)};

    vector<double> uvCoords;

    for( int i = 0; i < 3; i++) {
        uvCoords.push_back( uvCorners[i][0] );
        uvCoords.push_back( uvCorners[i][1] );
        auto edge = getEdgeAt(i);
        int  npoints =  edge->newNodes.size();
        if( npoints ) {
            int nsegments = npoints + 1;
            double dt = 1.0/(double)nsegments;
            for( int j = 0; j < npoints; j++) {
                double t = (j+1)*dt;
                double u = (1-t)*uvCorners[i][0] + t*uvCorners[(i+1)%3][0];
                double v = (1-t)*uvCorners[i][1] + t*uvCorners[(i+1)%3][1];
                uvCoords.push_back(u);
                uvCoords.push_back(v);
            }
        }
    }
    return uvCoords;
}
////////////////////////////////////////////////////////////////////////////////
void SurfRemesh::Face::unproject( const vector<double> &uvCoords)
{
    std::array<double,3> bCoords;
    std::array<Array3D,3> triCoords;
    triCoords[0] = nodes[0]->xyz;
    triCoords[1] = nodes[1]->xyz;
    triCoords[2] = nodes[2]->xyz;

    auto getXYZCoords = [] ( const std::array<Array3D,3> &tCoords, const Array3D &bCoord) {
        Array3D xyz;
        for( int i = 0; i < 3; i++)
            xyz[i] = bCoord[0]*tCoords[0][i] + bCoord[1]*tCoords[1][i] + bCoord[2]*tCoords[2][i];
        return  xyz;
    };

    int numNodes = uvCoords.size()/2;

    mesh = std::make_shared<Mesh>();
    mesh->nodes.reserve( numNodes);

    int istart =  3;
    for( int i = 0; i < 3; i++) {
        auto edge = getEdgeAt(i);
        int  np   = edge->newNodes.size();
        istart   += np;
        mesh->nodes.push_back( nodes[i] );
        int ori = edge->getOrientation(nodes[i], nodes[(i+1)%3] );
        assert( ori );
        if( ori > 0 ) {
            for( int j = 0; j < np; j++)
                mesh->nodes.push_back( edge->newNodes[j] ) ;
        } else {
            for( int j = 0; j < np; j++)
                mesh->nodes.push_back( edge->newNodes[np-j-1] ) ;
        }
    }

    Array2D qPoint;
    for( size_t i = istart; i < numNodes; i++) {
        qPoint[0]  = uvCoords[2*i];
        qPoint[1]  = uvCoords[2*i+1];
        bCoords    = barycoordinates(uvCorners[0], uvCorners[1], uvCorners[2], qPoint);
        auto vtx   = std::make_shared<Node>();
        vtx->xyz   = getXYZCoords( triCoords, bCoords);
        mesh->nodes.push_back(vtx);
    }
}
////////////////////////////////////////////////////////////////////////////////
void SurfRemesh::Face::delaunay(vector<int> &segments, vector<double> &uvCoords, vector<int> &triConnect)
{
    struct triangulateio in, out; // triangle sw data structure

    // The following stuff is for "Triangle" software ...
    in.numberofpoints = uvCoords.size()/2;
    in.pointlist = &uvCoords[0];

    double triarea[] = {area};

    in.numberofsegments = segments.size()/2;
    in.segmentlist  = &segments[0];

    in.numberofholes = 0;
    in.holelist = nullptr;

    in.numberofpointattributes = 0;
    in.pointattributelist = nullptr;
    in.pointmarkerlist = nullptr;
    in.segmentmarkerlist = nullptr;
    in.triangleattributelist = nullptr;
    in.numberofregions = 0;
    in.regionlist = nullptr;
    in.trianglearealist = triarea;

    // These needs to declared, I wasted full one day to find out this
    // restriction

    out.pointlist = nullptr;
    out.pointmarkerlist = nullptr;
    out.trianglelist = nullptr;
    out.neighborlist = nullptr;
    out.segmentlist = nullptr;
    out.segmentmarkerlist = nullptr;
    out.edgelist = nullptr;
    out.edgemarkerlist = nullptr;

    ostringstream oss;
    oss << "BYCPpzQq30";
    string options = oss.str();

    char *opt = const_cast<char*>( options.c_str() );
    triangulate(opt, &in, &out, (struct triangulateio *) nullptr);

    uvCoords.resize(2*out.numberofpoints);
    for (int i = 0; i < out.numberofpoints; i++) {
        uvCoords[2*i]   = out.pointlist[2*i + 0];
        uvCoords[2*i+1] = out.pointlist[2*i + 1];
    }

    triConnect.resize(3*out.numberoftriangles);
    for (int i = 0; i < out.numberoftriangles; i++) {
        triConnect[3*i+0] = out.trianglelist[3*i + 0];
        triConnect[3*i+1] = out.trianglelist[3*i + 1];
        triConnect[3*i+2] = out.trianglelist[3*i + 2];
    }

    if (out.pointlist) free(out.pointlist);
    if (out.trianglelist) free(out.trianglelist);
}

//////////////////////////////////////////////////////////////////////////////////

void SurfRemesh::Face::refine(double maxarea)
{
    area = maxarea;
    vector<double> uvCoords = project();

    int nPoints = uvCoords.size()/2;

    vector<int> segments(2*nPoints);

    int index = 0;
    for( int i = 0; i < nPoints; i++) {
        segments[2*i]   = index;
        segments[2*i+1] = (index+1)%nPoints;
        index++;
    }

    vector<int> trifaces;
    delaunay(segments, uvCoords, trifaces);
    unproject( uvCoords);

    int numTris = trifaces.size()/3;
    mesh->faces.resize(numTris);
    for( int i = 0; i < numTris; i++) {
        auto n0 = mesh->nodes[trifaces[3*i+0]];
        auto n1 = mesh->nodes[trifaces[3*i+1]];
        auto n2 = mesh->nodes[trifaces[3*i+2]];
        auto face = std::make_shared<Face>();
        face->nodes = {n0,n1,n2};
        mesh->faces[i] = face;
    }
}

//////////////////////////////////////////////////////////////////////////////////

void SurfRemesh::buildEdges()
{
    for(auto vtx : mesh.nodes) vtx->edges.clear();

    // Calling getEdgeAt will build new edges and store on the nodes ...
    for( auto face : mesh.faces) {
        for( int i = 0; i < 3; i++) face->getEdgeAt(i);
    }

    // Now collect unique edges from the nodes ...
    mesh.edges.clear();
    for( auto vtx : mesh.nodes) {
        for( auto e : vtx->edges) mesh.edges.push_back(e);
    }
}

//////////////////////////////////////////////////////////////////////////////////

void SurfRemesh::refine()
{
    if( mesh.nodes.empty() ) readOFFMesh( infilename);

    while(1) {

        // Generate unique edges from the faces..
        buildEdges();

        // First refine all the edges ...
        size_t ncount = 0;
        for( auto edge : mesh.edges) ncount += edge->refine(maxEdgeLength);

	cout << "#Edges refined " << ncount << endl;

        if( ncount == 0) break;

	double maxarea = sqrt(3)*0.25*maxEdgeLength*maxEdgeLength;
//#pragma omp parallel for
        for( size_t i = 0; i < mesh.faces.size(); i++) 
             mesh.faces[i]->refine(maxarea);

        // Count number of nodes in the output (orginal + edgenodes + facenodes)
        size_t countnodes = mesh.nodes.size(); // from the input mesh.
        for(auto edge : mesh.edges) countnodes += edge->newNodes.size();

        // Count number of faces in the output ...
        for(auto face : mesh.faces) countnodes += face->mesh->nodes.size();

        size_t countfaces = mesh.faces.size();
        for( auto face : mesh.faces) countfaces += face->mesh->faces.size();

        mesh.nodes.reserve( countnodes);

        // Collect new nodes on the edges ...
        int offset = mesh.nodes.size();
        for( auto edge : mesh.edges) {
            for( auto vtx : edge->newNodes) {
                vtx->id = offset++;
                mesh.nodes.push_back(vtx);
            }
        }
        mesh.edges.clear();

        // Collect new nodes on the faces ...
        vector<FacePtr> allfaces;
        for( auto face : mesh.faces) {
            assert( face->mesh );
            for( auto vtx : face->mesh->nodes) {
                vtx->id = offset++;
                mesh.nodes.push_back(vtx);
            }
            for( auto f : face->mesh->faces)
                allfaces.push_back(f);
            face->mesh = nullptr;
        }
        mesh.faces = allfaces;
    }
}

//////////////////////////////////////////////////////////////////////////////////

void SurfRemesh:: saveAs( const string &filename)
{
    ofstream ofile( filename.c_str(), ios::out);
    ofile << "OFF" << endl;

    ofile <<  mesh.nodes.size() << " " << mesh.faces.size() << " 0 " << endl;

    for( auto v : mesh.nodes)
        ofile << v->xyz[0] << " " << v->xyz[1] << " " << v->xyz[2] << endl;

    for( auto face : mesh.faces) {
        ofile << "3 " << face->nodes[0]->id << " "
              << face->nodes[1]->id << " "
              << face->nodes[2]->id << endl;
    }
}

//////////////////////////////////////////////////////////////////////////////////
