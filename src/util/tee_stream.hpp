/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

/**
 * Adapted from <a href="http://wordaligned.org/articles/cpp-streambufs">
 * http://wordaligned.org/articles/cpp-streambufs</a>
 */
#pragma once

#include <iostream>
#include <streambuf>

template <typename char_type,
          typename traits = std::char_traits<char_type> >
class basic_teebuf:
    public std::basic_streambuf<char_type, traits>
{
public:
    typedef typename traits::int_type int_type;

    basic_teebuf(std::basic_streambuf<char_type, traits> * sb1,
                 std::basic_streambuf<char_type, traits> * sb2)
      : sb1_(sb1)
      , sb2_(sb2)
    {
        return;
    }

    void reset(std::basic_streambuf<char_type, traits> * sb1,
               std::basic_streambuf<char_type, traits> * sb2)
    {
        sb1_ = sb1;
        sb2_ = sb2;
        return;
    }
    
private:    
    virtual int sync()
    {
        int const r1 = sb1_->pubsync();
        int const r2 = sb2_->pubsync();
        return r1 == 0 && r2 == 0 ? 0 : -1;
    }
    
    virtual int_type overflow(int_type c)
    {
        int_type const eof = traits::eof();
        
        if (traits::eq_int_type(c, eof))
        {
            return traits::not_eof(c);
        }
        else
        {
            char_type const ch = traits::to_char_type(c);
            int_type const r1 = sb1_->sputc(ch);
            int_type const r2 = sb2_->sputc(ch);
            
            return
                traits::eq_int_type(r1, eof) ||
                traits::eq_int_type(r2, eof) ? eof : c;
        }
    }

    
private:
    std::basic_streambuf<char_type, traits> * sb1_;
    std::basic_streambuf<char_type, traits> * sb2_;
};

typedef basic_teebuf<char> teebuf;

/**
 * A simple output stream buffer which performs no actions. This is useful as an
 * argument to the constructor for \ref TeeStream and \ref TeeStream::reset()
 * for when no extra stream is needed (for instance, when there is no log file).
 *
 * Source:
 * http://stackoverflow.com/questions/760301/implementing-a-no-op-stdostream
 */
template <class cT, class traits = std::char_traits<cT> >
class basic_nullbuf: public std::basic_streambuf<cT, traits> {
    typename traits::int_type overflow(typename traits::int_type c)
    {
        return traits::not_eof(c); // indicate success
    }
};

template <class cT, class traits = std::char_traits<cT> >
class basic_onullstream: public std::basic_ostream<cT, traits> {
    public:
        basic_onullstream():
        std::basic_ios<cT, traits>(&m_sbuf),
        std::basic_ostream<cT, traits>(&m_sbuf)
        {
            this->init(&m_sbuf);
        }

    private:
        basic_nullbuf<cT, traits> m_sbuf;
};

typedef basic_onullstream<char> onullstream;
typedef basic_onullstream<wchar_t> wonullstream;


/**
 * An output stream that directs output to two stream buffers. This is super
 * useful for directing output to both standard output and a log file at the
 * same time.
 */
class TeeStream : public std::ostream
{
public:
    /**
     * \brief Construct a TeeStream using two output streams.
     */
    TeeStream(std::ostream & o1, std::ostream & o2) :
		std::ostream(&tbuf),
        tbuf(o1.rdbuf(), o2.rdbuf())
	{
        return;
	}

    void reset( std::ostream &o1, std::ostream &o2 ) {
        tbuf.reset(o1.rdbuf(), o2.rdbuf());
        return;
    }
private:
    teebuf tbuf;
};
