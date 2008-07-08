// 
// File:   rangeparser.h
// Author: ruehle
//
// Created on March 28, 2008, 10:50 AM
//

#ifndef _RANGEPARSER_H
#define	_RANGEPARSER_H

#include <list>
#include <string>

using namespace std;

class RangeParser
{
public:
    RangeParser();
    
    void Parse(string range);

private:
    struct block_t {
        block_t() {}
        block_t(const int &begin, const int &end, const int &stride)
        : _begin(begin), _end(end), _stride(stride) {}
        
        int _begin, _end, _stride;
    };

public:
    
    struct iterator {
        iterator() {}
        
        int operator*() const
            { return _current; }
        
        RangeParser::iterator &operator++();                      
        
        bool operator==(const RangeParser::iterator &);                      
        bool operator!=(const RangeParser::iterator &);                      
         
    private:
        iterator(RangeParser *, list<block_t>::iterator);
        list<block_t>::iterator _block;
        int _current;
        
        RangeParser *_parent;
        
        friend class RangeParser;
    };
    
    RangeParser::iterator begin();
    RangeParser::iterator end();
       
private:
    
    void ParseBlock(string block);
    int ToNumber(string str);
        
    list< block_t > _blocks;    
    
    bool _has_begin, _has_end;
    int _begin, _end;
    
};

inline RangeParser::iterator RangeParser::begin()
{
    return RangeParser::iterator(this, _blocks.begin());
}

inline RangeParser::iterator RangeParser::end()
{
    return RangeParser::iterator(this, _blocks.end());
}

inline RangeParser::iterator::iterator(RangeParser *parent, list<block_t>::iterator block)
    : _parent(parent), _block(block)
{
    if(block != parent->_blocks.end())
        _current = (*block)._begin;
    else
        _current = -1;
}
            
inline bool RangeParser::iterator::operator==(const RangeParser::iterator &i)
{
    return  (_block == i._block) && (_current == i._current);
}

inline bool RangeParser::iterator::operator!=(const RangeParser::iterator &i)
{
    return  !((_block == i._block) && (_current == i._current));
}


#endif	/* _RANGEPARSER_H */

