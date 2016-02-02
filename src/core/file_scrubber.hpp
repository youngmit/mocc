#pragma once
#include <string>
#include <fstream>

namespace mocc {
    /**
     * A simple class that extends ifstream to provide a safe/easy way to scrub
     * comments from an input file as it is read
     */
    class FileScrubber {
    public:
        FileScrubber(){}
        /**
         * Ininialize from a file name. Passes fName into the constructor of
         * the underlying ifstream object.
         *
         * \param fName the name of the file to read.
         * \param commentFlag a C-style string containing the sequence of
         * characters that signify a comment.
         *
         * This class only supports single-line comments, in which the \p
         * commentFlag and all characters following are ignored for the rest of
         * the line. Think '//' comments in C-style languages
         */
        FileScrubber(const char* fName, const char* commentFlag);

        /**
         * Extend getline to only return lines which are not empty after
         * removing comments, stripped of comments.
         */
        std::string getline();

        // Wrap the eof() function on the underlying ifstream object
        bool eof() {
            return stream_.eof();
        }

    private:
        std::ifstream stream_;
        std::string flag_;
    };
}
