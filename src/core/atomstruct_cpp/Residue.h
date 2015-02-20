// vi: set expandtab ts=4 sw=4:
#ifndef atomstruct_Residue
#define atomstruct_Residue

#include <vector>
#include <map>
#include <string>
#include "imex.h"
#include "string_types.h"

namespace atomstruct {

class Atom;
class AtomicStructure;
class Bond;

class ATOMSTRUCT_IMEX Residue {
    friend class AtomicStructure;
public:
    typedef std::vector<Atom *>  Atoms;
    typedef std::multimap<AtomName, Atom *>  AtomsMap;
private:
    Residue(AtomicStructure *as, const std::string &name, const std::string &chain, int pos, char insert);
    char  _alt_loc;
    Atoms  _atoms;
    std::string  _chain_id;
    char  _insertion_code;
    bool  _is_helix;
    bool  _is_het;
    bool  _is_sheet;
    std::string  _name;
    int  _position;
    int  _ss_id;
    AtomicStructure *  _structure;
public:
    void  add_atom(Atom *);
    const Atoms &  atoms() const { return _atoms; }
    AtomsMap  atoms_map() const;
    std::vector<Bond*>  bonds_between(const Residue* other_res,
        bool just_first=false) const;
    const std::string &  chain_id() const { return _chain_id; }
    int  count_atom(const AtomName&) const;
    Atom *  find_atom(const AtomName&) const;
    char  insertion_code() const { return _insertion_code; }
    bool  is_het() const { return _is_het; }
    const std::string &  name() const { return _name; }
    int  position() const { return _position; }
    void  set_alt_loc(char alt_loc);
    void  set_is_helix(bool ih) { _is_helix = ih; }
    void  set_is_het(bool ih) { _is_het = ih; }
    void  set_is_sheet(bool is) { _is_sheet = is; }
    void  set_ss_id(int ssid) { _ss_id = ssid; }
    std::string  str() const;
    AtomicStructure*  structure() const { return _structure; }
    std::vector<Atom*>  template_assign(
        void (Atom::*assign_func)(const char*), const char* app,
        const char* template_dir, const char* extension) const;
};

}  // namespace atomstruct

#endif  // atomstruct_Residue
