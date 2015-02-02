// vi: set expandtab ts=4 sw=4:
#ifndef atomstruct_Sequence
#define atomstruct_Sequence

#include <vector>
#include <unordered_map>
#include <string>
#include "imex.h"

namespace atomstruct {

class ATOMSTRUCT_IMEX Sequence {
public:
    typedef std::vector<char>  Contents;
protected:
    typedef std::unordered_map<std::string, char>  _1Letter_Map;
    static void  _init_rname_map();
    static _1Letter_Map  _nucleic3to1;
    static _1Letter_Map  _protein3to1;
    static _1Letter_Map  _rname3to1;

    mutable std::unordered_map<unsigned int, unsigned int>  _cache_g2ug;
    mutable std::unordered_map<unsigned int, unsigned int>  _cache_ug2g;
    mutable Contents  _cache_ungapped;
    void  _clear_cache() const
        { _cache_ungapped.clear();  _cache_g2ug.clear(); _cache_ug2g.clear(); }
    // can't inherit from vector, since we need to clear caches on changes
    Contents  _contents;
public:
    static void  assign_rname3to1(const std::string& rname, unsigned char let,
        bool protein);
    static unsigned char  nucleic3to1(const std::string &rn);
    static unsigned char  protein3to1(const std::string &rn);
    static unsigned char  rname3to1(const std::string &rn);

    template <class InputIterator> void  assign(InputIterator first,
        InputIterator last) { _clear_cache(); _contents.assign(first, last); }
    Contents::reference  at(Contents::size_type n)
        { _clear_cache(); return _contents.at(n); }
    Contents::const_iterator  begin() const { return _contents.begin(); }
    void  clear() { _clear_cache(); _contents.clear(); }
    Contents::const_iterator  end() const { return _contents.end(); }
    unsigned int  gapped_to_ungapped(unsigned int index) const;
    Contents::iterator  insert(Contents::const_iterator pos,
        Contents::size_type n, Contents::value_type val)
        { _clear_cache(); return _contents.insert(pos, n, val); }
    void  push_back(unsigned char c) { _clear_cache(); _contents.push_back(c); }
    Contents::const_reverse_iterator  rbegin() const
        { return _contents.rbegin(); }
    Contents::const_reverse_iterator  rend() const { return _contents.rend(); }
    Sequence() {}
    Sequence(const std::vector<std::string>& res_names);  // 3-letter codes
    virtual  ~Sequence() {}
    Contents::size_type  size() const { return _contents.size(); }
    void  swap(Contents& x) { _clear_cache(); _contents.swap(x); }
    const Contents&  ungapped() const;
    unsigned int  ungapped_to_gapped(unsigned int index) const;
};

}  // namespace atomstruct

#endif  // atomstruct_Sequence
