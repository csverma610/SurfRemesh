#include "SurfRemesh.h"
#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
    if( argc == 1) {
        cout << "Usage: " << argv[0] << " inmeshfile edgelengh outmeshfile " << endl;
        return 1;
    }

    SurfRemesh remesh;
    remesh.setMesh( argv[1] );
    double elen = atof(argv[2]);
    remesh.setMaxEdgeLength( elen );
    remesh.refine();
    remesh.saveAs( argv[3] );
    return 0;
}
