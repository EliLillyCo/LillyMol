#ifndef MOLECULE_LIB_COORDINATE_BOX_H_
#define MOLECULE_LIB_COORDINATE_BOX_H_

#include <cstdint>
#include <tuple>

#include "Foundational/iwstring/iwstring.h"
#include "Molecule_Lib/space_vector.h"

namespace coordinate_box {

// The idea is that by establishing a box around a molecule, and a grid,
// then the position of any atom can be specified by a single number. As long
// as there is a means of transforming from that single number to a cell
// within the box.
// Since this is 3 dimensions, there are two dimensions for the box - it
// is unconstrainted in the z dimension.
// Set CELL_NUMBER_CHECK_IN_BOUNDS to enable bounds checking.
class CoordinateBox {
  private:
    // The length of each cell.
    double _cell_size;

    // The number of cells along each side of the box. Z is unspecified.
    int _cell[2];

    // Once a position in the box is determined, that then gets
    // translated to where the origin is.
    Space_Vector<double> _origin;

  public:
    CoordinateBox();

    // Syntax is %B{cell_size,x,y,z}. This function assumes
    // just the csv values.
    int BuildFromSmilesToken(const const_IWSubstring& token);

    // For a position in space, the corresponding cell.
    template <typename T>
    int CellNumber(const Space_Vector<T>& coords) const {
      return CellNumber(coords.x(), coords.y(), coords.z());
    }
    template <typename T>
    int CellNumber(T x, T y, T z) const;


    // For a cell number, return the position in space.
    template <typename T> 
    std::tuple<T, T, T> CoordinatesAsTuple(int cell_number) const;
    template <typename T>
    Space_Vector<T> CoordinatesAsVector(int cell_number) const;
};

// The concentric box positions points in space via two numbers,
// 1. Layer Number
// 2. Position within the layer
// Beware possible compiler dependent alignment issues.
struct LayerPosition {
  uint32_t layer;
  uint64_t position_in_layer;
  LayerPosition() {
    layer = 0;
    position_in_layer = 0;
  }
  LayerPosition(uint32_t l, uint64_t p) : layer(l), position_in_layer(p) {
  }
};

// A cell number that defines a set of coordinates in relation to a
// concentric set of cubic shells. The shells are defined by a resolution
// specified in the constructor.
// The cell numbering is not public.
// Points going through the conversion process will be truncated at
// values that are multiples of `resolution`.
// Note that to avoid possible overflow issues, prefer to use the API
// calls involving the LayerPosition struct.
// Typical usage:
// coordinate_box::ConcentricBox box;
// Space_Vector<float> original_coords(x, y, z);
// const LayerPosition lpip = box.Position(coords);
// const Space_Vector<float> recovered = box.LayerAndPositionToCoords(lpip);
// expect diff between original and recovered to be small.
class ConcentricBox {
  private:
    double _resolution;
    // Several computations need _resolution*0.5
    double _half_resolution;

  // private functions.
    template <typename T> Space_Vector<T> _right_face(uint32_t layer, uint64_t within_layer) const;
    template <typename T> Space_Vector<T> _left_face(uint32_t layer, uint64_t within_layer) const;
    template <typename T> Space_Vector<T> _somewhere_along_x(uint32_t layer, uint64_t within_layer) const;

    template <typename T> uint64_t CellNumber(const Space_Vector<T>& coords, uint32_t layer, uint32_t x, uint32_t y, uint32_t z) const;
    template <typename T> uint32_t LayerNumber(const Space_Vector<T>& coords) const;
    template <typename T> LayerPosition Position(const Space_Vector<T>& coords,
                        uint32_t layer,
                        uint64_t x, uint64_t y, uint64_t z) const;
    template <typename T> uint64_t PositionInLayer(const Space_Vector<T>& coords,
                        uint32_t layer,
                        uint64_t x, uint64_t y, uint64_t z) const;
    template <typename T> Space_Vector<T> LayerAndPositionToCoords(uint32_t layer, uint64_t within_layer) const;

  public:
    // A default resolution which should be good enough for most interatomic distances.
    explicit ConcentricBox(double res = 0.001);

    double resolution() const { return _resolution;}

    // Given a position in space, return the corresponding LayerPosition
    template <typename T>
    LayerPosition Position(const Space_Vector<T>& coords) const;
    // Given a position in space, return the corresponding cell number.
    template <typename T>
    uint64_t CellNumber(const Space_Vector<T>& coords) const;

    // For a cell number, return the position in space.
    template <typename T>
    Space_Vector<T> CellToCoordinates(const LayerPosition& layer_pos_in_layer) const;
    template <typename T>
    Space_Vector<T> CellNumberToCoordinates(uint64_t cell_number) const;
};

// Layer and positions are written as
// layer:position.
// Parse 's' for a leading `layer:position` form. If successful
// returns the number of characters consumed, and fills in `destination`.
// If unsuccessful, 0 is returned and `destination` is indeterminate.
int FromString(const const_IWSubstring& s, LayerPosition& destination);

IWString &
operator << (IWString & output, const LayerPosition& layer_position);

// Helper function for ConcentricBox, exposed just for testing.
uint32_t CellToLayer(uint64_t cell);

} // namespace coordinate_box

#endif // MOLECULE_LIB_COORDINATE_BOX_H_
