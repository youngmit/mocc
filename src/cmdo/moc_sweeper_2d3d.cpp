#include "moc_sweeper_2d3d.hpp"

#include <iostream>

using std::cout;
using std::cin;
using std::endl;

namespace mocc {
    MoCSweeper_2D3D::MoCSweeper_2D3D( const pugi::xml_node &input, 
            const CoreMesh &mesh ):
        MoCSweeper( input, mesh ),
        corrections_( nullptr )
    {
        
    };

    void MoCSweeper_2D3D::sweep1g_final( int group ) {
        cout << "2D3D final sweep" << endl;
        cin.ignore();
        return;
    }
}
