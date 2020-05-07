#include <array>
#include <vector>
#include <memory>

#include "veclib.hpp"

class SurfRemesh
{
    class Node;
    using NodePtr = std::shared_ptr<Node>;

    class Edge;
    using EdgePtr = std::shared_ptr<Edge>;

    class Face;
    using FacePtr = std::shared_ptr<Face>;

    class Mesh;
    using MeshPtr = std::shared_ptr<Mesh>;

    struct Node
    {
        size_t  id;
        Array3D xyz;
        std::vector<EdgePtr> edges;
    };

    struct Edge
    {
        int getOrientation( const NodePtr &v0, const NodePtr &v1) const
        {
            if( nodes[0] == v0  && nodes[1] == v1) return 1;
            if( nodes[0] == v1  && nodes[1] == v0) return -1;
            return 0;
        }

        bool   isSame( const NodePtr &v0, const NodePtr &v1) const
        {
            if( nodes[0] == v0  && nodes[1] == v1) return true;
            if( nodes[0] == v1  && nodes[1] == v0) return true;
            return false;
        }

        int refine( double len);

        std::array<NodePtr,2> nodes;
        std::vector<NodePtr>  newNodes;
    };

    struct Face
    {
        void refine(double len);
        void delaunay(std::vector<int> &segments, std::vector<double> &uvCoords, std::vector<int> &connect);

        EdgePtr getEdgeAt(int i);

        std::vector<double>   project();
        void unproject( const std::vector<double> &uv);

        MeshPtr mesh;
        std::array<Array2D,3> uvCorners;
        std::array<NodePtr,3> nodes;
	double area;
    };

    struct Mesh
    {
        void clear() {
            nodes.clear();
            edges.clear();
            faces.clear();
        }
        std::vector<NodePtr> nodes;
        std::vector<EdgePtr> edges;
        std::vector<FacePtr> faces;
    };

public:
    SurfRemesh();

    // Set the mesh (only off format) ...
    void setMesh( const std::string &f) {
        infilename = f;
    }

    // If we do not know the length of the model, then we can decide
    // the maxEdgeLength based on how many times the longest edge 
    // is divided. This can set the maxEdgeLength of the model.
    // Howwever, if the maxEdgeLength is spcified, then this parameter
    // is ignored.

    void setMaxSteps(int n) {
	    maxSteps = n;
    }

    // Set the desired maximum edge length ...
    void setMaxEdgeLength( double l) {
        maxEdgeLength = l;
    }

    // Set the desired minimum edge length ...
    void setMinEdgeLength( double l) {
        minEdgeLength = l;
    }

    // Set the crease angle..
    void setCreaseAngle( double a) {
        creaseAngle - a;
    }

    // Collected (minimum, mean, maximum) edge lengths of the model ...
    std::array<double,3>  getEdgeLengths() const;

    // Start executing the code ...
    void refine();

    // Return the refined mesh nodes and triangles ...
    std::vector<Array3D> getNodes() const;
    std::vector<Array3I> getTriangles() const;

    // Store the refined mesh ...
    void saveAs( const std::string &s);

private:
    Mesh mesh;
    std::string infilename;
    double maxEdgeLength   = 0.0;
    double minEdgeLength   = 0.0;
    double creaseAngle     = 0.0;
    int    maxSteps        = 0;

    void   readOFFMesh( const std::string &f);
    void   buildEdges();
};
