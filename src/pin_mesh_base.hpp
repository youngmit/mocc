#pragma once

#include <iostream>

#include "pugixml.hpp"

#include "global_config.hpp"

namespace mocc {
    class PinMesh{
    public:
        PinMesh( const pugi::xml_node &input );

        virtual ~PinMesh() {
        }
    	
        int id() const {
    		return id_;
    	}

        int n_reg() const {
            return n_reg_;
        }

        int n_xsreg() const {
            return n_xsreg_;
        }

        float_t pitch_x() const {
            return pitch_x_;
        }

        float_t pitch_y() const {
            return pitch_y_;
        }
    protected:
    	int id_;
        int n_reg_;
        int n_xsreg_;
    	float_t pitch_x_;
    	float_t pitch_y_;
    };
}
