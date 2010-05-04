/**
 *  @file
 *  @brief Simple c++ interface to regular expressions.
 *  @author Doug Potter
 */

#ifndef REGEX_HPP
#define REGEX_HPP

#include "regex.h"
class Regex : protected re_pattern_buffer, protected re_registers, protected std::string {
    bool bCompiled;
public:
    Regex();
    Regex( const std::string &pattern );
    virtual ~Regex();

    int Search( const std::string &text );
    int Match(  const std::string &text );

    inline size_t Count() const { return re_nsub; }

    std::string operator[](size_t idx) const;
    std::string GetText() const;
};

Regex::Regex()
{
    bCompiled = false;
}

Regex::Regex( const std::string &pattern )
{
    bCompiled = true;
    re_set_syntax(RE_SYNTAX_POSIX_EXTENDED);

    buffer = 0;
    allocated = 0;
    fastmap = 0;
    translate = 0;
    no_sub = 0;

    re_compile_pattern(pattern.c_str(),(int)pattern.size(),this);
    re_set_registers( this, this, 0, 0, 0 );
}

Regex::~Regex()
{
    if ( bCompiled )
	regfree(this);
}

std::string Regex::operator[](size_t idx) const
{
    assert(idx < re_nsub);
    if ( start[idx+1] < 0 ) return "";
    return substr( start[idx+1],re_registers::end[idx+1]-start[idx+1] );
}

std::string Regex::GetText() const
{
    return *this;
}

int Regex::Search( const std::string &text )
{
    (*(std::string *)this) = text;
    return re_search(this,c_str(),(int)size(),0,(int)size(),this);
}

int Regex::Match( const std::string &text )
{
    (*(std::string *)this) = text;
    return re_match(this,c_str(),(int)size(),0,this);
}


#endif
