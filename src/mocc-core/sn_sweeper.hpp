#pragma once

#include "pugixml.hpp"

#include "mesh.hpp"
#include "sn_source.hpp"
#include "core_mesh.hpp"
#include "global_config.hpp"
#include "transport_sweeper.hpp"
#include "angular_quadrature.hpp"
#include "xs_mesh_homogenized.hpp"
#include "coarse_data.hpp"

namespace mocc {
    class SnSweeperBoundary {
    public:
        SnSweeperBoundary() { }
        SnSweeperBoundary( int n_grp, int n_ang , int nx, int ny, int nz):
            n_ang_( n_ang ),
            nx_( nx ),
            ny_( ny ),
            nz_( nz ),
            ang_stride_( nx*ny + nx*nz + ny*nz ),
            data_( (nx*ny + nx*nz + ny*nz) * n_ang * n_grp, 0.0 )
        {
            n_face_[(int)Normal::X_NORM] = ny_*nz_;
            n_face_[(int)Normal::Y_NORM] = nx_*nz_;
            n_face_[(int)Normal::Z_NORM] = nx_*ny_;
            int offset = 0;
            for( const auto norm: AllNormals ) {
                face_offset_[(int)norm] = offset;
                offset += n_face_[(int)norm];
            }
            return;
        }

        void get_face( int grp, int ang, Normal norm, float_t* out ) const {
            const float_t* data = data_.data() + ang_stride_*n_ang_*grp + 
                ang_stride_*ang + face_offset_[(int)norm];
            for( int i=0; i<n_face_[(int)norm]; i++ ) {
                out[i] = data[i];
            }
            return; 
        }

        const float_t* get_face_ptr( int grp, int ang, Normal norm ) const {
            return data_.data() + ang_stride_*ang*grp + 
                ang_stride_*ang + face_offset_[(int)norm];
        }

        void set_face( int grp, int ang, Normal norm, const float_t *in ) {
            float_t* data = data_.data() + ang_stride_*n_ang_*grp + 
                ang_stride_*ang + face_offset_[(int)norm];
            for( int i=0; i<n_face_[(int)norm]; i++ ) {
                data[i] = in[i];
            }
        }

        void zero_face( int grp, int ang, Normal norm ) {
            float_t* data = data_.data() + ang_stride_*n_ang_*grp + 
                ang_stride_*ang + face_offset_[(int)norm];
            for( int i=0; i<n_face_[(int)norm]; i++ ) {
                data[i] = 0.0;
            }
        }

        void initialize( float_t val ) {
            for( auto &v: data_ ) {
                v = val;
            }
        }

    private:
        int n_ang_;
        int nx_;
        int ny_;
        int nz_;
        int ang_stride_;
        int face_offset_[3];
        int n_face_[3];
        VecF data_;
    };

    class SnSweeper: public TransportSweeper {
    public:
        SnSweeper( const pugi::xml_node& input, const CoreMesh& mesh );

        ~SnSweeper() { }

        void sweep( int group );

        void initialize();
        void get_pin_flux( int ig, VecF& flux ) const;

        void output( H5File& file ) const;

        // Override the create_source() method to make an SnSource instead of
        // the regular
        UP_Source_t create_source() {
            Source *s = new SnSource( n_reg_, xs_mesh_.get(), this->cflux());
            UP_Source_t source( s );
            return source;
        }

        void homogenize( CoarseData &data ) const;
    private:
        // Update the boundary conditions 
        void update_boundary( int group );

        const CoreMesh& core_mesh_;
        unsigned int n_inner_;
        AngularQuadrature ang_quad_;

        // Boundary condition enumeration
        std::vector<Boundary> bc_type_;
        
        // Temporary storage for 1-group scalar flux
        ArrayX flux_1g_;

        // Mesh parameters
        int nx_;
        int ny_;
        int nz_;
        VecF hx_;
        VecF hy_;
        VecF hz_;

        // Temporary storage of the current-group transport cross section
        ArrayX xstr_;

        // Single-group isotropic source, should include in-scatter
        ArrayX q_;

        // Incomming boundary condition
        SnSweeperBoundary bc_in_;

        // Outgoing boundary condition. Only difined for one group
        SnSweeperBoundary bc_out_;

        void sweep_std( int group );
    };
}
