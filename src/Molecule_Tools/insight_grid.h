#ifndef MOLECULE_TOOLS_INSIGHT_GRID_H
#define MOLECULE_TOOLS_INSIGHT_GRID_H

#include <cstdint>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/data_source/iwstring_data_source.h"

#include "Molecule_Lib/molecule.h"

namespace insight_grid {

class InsightGrid {
  private:
    float _xlen;
    float _ylen;
    float _zlen;

    double _dx;
    double _dy;
    double _dz;

    double _xmin;
    double _ymin;
    double _zmin;
    double _xmax;
    double _ymax;
    double _zmax;

    uint32_t _x;
    uint32_t _y;
    uint32_t _z;

    float * _data;
    uint32_t _npoints;

    enum class Axis {
      kXaxis,
      kYaxis,
      kZaxis,
    };

    // When forming a grid molecule in MakeGridMolecule, we can isotopically
    // label each atom in the grid with a multipe of the value of _data at
    // each point. By default this is zero, which means no isotope.
    float _isotope_multiplier;

  // private functions.

    int GetDimensions(const const_IWSubstring& buffer);
    int GetResolution(const const_IWSubstring& buffer);
    int GetGridOffsets(const const_IWSubstring& buffer);
    int MakeGridMolecule() const;

  public:
    InsightGrid();
    ~InsightGrid();

    int DebugPrint(std::ostream& output) const;

    int Build(const char* fname);
    int Build(iwstring_data_source& input);

    uint32_t npoints() const {
      return _npoints;
    }

    // Given 3 coordinates, return the index of that pont in _data;
    int IndexInData(float x, float y, float z) const;

    // Given a point in space (x,y,z) and a distance around it,
    // set all the corresponding indices in _data in `in_range`.
    // Note however that indices where _data is zero are skipped.
    int IndicesInRange(float x, float y, float z, float distance, int* in_range) const;

    // Return in `acc` the sum of all data values that are within `distance` of
    // (x,y,z).
    // Returns 1 if any values are added to `acc` - it is not reset. Designed
    // to be called across multiple atoms.
    int Sum(float x, float y, float z, float distance, Accumulator<double>& acc) const;

    // `ndx` must be an array of dimension NPoints. For each _data[ndx] that is set
    // accumulate _data[ndx] in `acc`.
    void Score(const int * ndx, Accumulator<double>& acc) const;

    void set_isotope_multiplier(float s) {
      _isotope_multiplier = s;
    }

    // Create a molecule, `grid_mol` where each grid point becomes an Atom
    // of type `element`. `grid_mol` will not have any atoms that are closer
    // than `threshold` from any atom in `m`.
    int MakeMolecule(const Element* element, float threshold,
                     const Molecule& m,
                     Molecule& grid_mol) const;
};

}  // namespace insight_grid
#endif  //  MOLECULE_TOOLS_INSIGHT_GRID_H
