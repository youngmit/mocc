#include "string_utils.hpp"

#include <sstream>
#include <ctype.h>

#include "core/global_config.hpp"
#include "core/error.hpp"

using mocc::Exception;

using std::cout;
using std::endl;

std::string print_range( const std::vector<bool> &input ) {
    std::stringstream output;
    bool left = false;
    int left_bound = 0;
    std::vector< std::pair<int, int> > ranges;
    for( int i=0; i<(int)input.size(); i++ ) {
        if( input[i] && !left ) {
            left = true;
            left_bound = i;
        }
        if( left && !input[i] ) {
            left = false;
            ranges.emplace_back(left_bound, i-1);
        }
        if( left && input[i] && (i == (int)input.size()-1) ) {
            ranges.emplace_back(left_bound, i);
        }
    }

    for( auto r: ranges ) {
        if( r.first != r.second ) {
            output << r.first << "-" << r.second;
        } else {
            output << r.first;
        }

        if( r != ranges.back() ) {
            output << ", ";
        }
    }
    return output.str();
}

template<typename T>
std::vector<T> explode_string(std::string data) {
    sanitize(data);
    auto is_invalid_char = [](char c) {
        return !(isspace(c) || isdigit(c));
    };

    // first, make sure there are no non-[numerals|whitespace]
    bool good = std::find_if(data.begin(), data.end(), is_invalid_char) 
        == data.end();
    if(!good) {
        throw EXCEPT("Malformed data");
    }


    std::vector<T> out;
    
    // Store the string as a stream and try to read all entries
    std::stringstream inBuf(data);
    T i;
    while ( !inBuf.eof() ) {
        inBuf >> i;
        out.push_back(i);
    }

    // Make sure everything went okay
    if ( inBuf.fail() ) {
        throw EXCEPT( "Trouble reading data" );
    }

    return out;
}

template std::vector<int> explode_string(std::string data);
template std::vector<mocc::real_t> explode_string(std::string data);
