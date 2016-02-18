#pragma once

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>


// From http://stackoverflow.com/a/25829233

// trim from left
inline std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(0, s.find_first_not_of(t));
    return s;
}

// trim from right
inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

// trim from left & right
inline std::string& trim(std::string& s, const char* t = " \t\n\r\f\v")
{
    return ltrim(rtrim(s, t), t);
}

// copying versions

inline std::string ltrim_copy(std::string s, const char* t = " \t\n\r\f\v")
{
    return ltrim(s, t);
}

inline std::string rtrim_copy(std::string s, const char* t = " \t\n\r\f\v")
{
    return rtrim(s, t);
}

inline std::string trim_copy(std::string s, const char* t = " \t\n\r\f\v")
{
    return trim(s, t);
}

/**
 * \brief Return a string, representing the ranges of a std::vector<bool> that
 * are true.
 */
inline std::string print_range( const std::vector<bool> &input ) {
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


// Sanitize a string: remove whitespace and cast to lowercase.
inline std::string& sanitize(std::string &s)
{
    std::transform( s.begin(), s.end(), s.begin(),
                       ::tolower );
    return trim(s);
}

