#pragma once

#include <iostream>

#include "pugixml.hpp"

#include "global_config.hpp"

namespace mocc {
    class PinMesh{
    public:
        PinMesh( const pugi::xml_node &input );

        virtual ~PinMesh() {
            std::cout << "Destroying PinMesh" << std::endl;
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
    protected:
    	int id_;
        int n_reg_;
        int n_xsreg_;
    	float_t pitch_x_;
    	float_t pitch_y_;
    };
}
