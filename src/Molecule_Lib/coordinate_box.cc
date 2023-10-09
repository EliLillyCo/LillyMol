#include <cmath>

#include "coordinate_box.h"

namespace coordinate_box {
using std::cerr;

CoordinateBox::CoordinateBox() {
  _cell_size = 0.0;
  _cell[0] = 0.0;
  _cell[1] = 0.0;
}

int
CoordinateBox::BuildFromSmilesToken(const const_IWSubstring& from_smiles) {
  resizable_array_p<const_IWSubstring> tokens;
  if (from_smiles.split(tokens, ',') != 3) {
    cerr << "CoordinateBox::BuildFromSmilesToken:must contain 3 tokens " << from_smiles << '\n';
    return 0;
  }

  for (const const_IWSubstring* token : tokens) {
    if (token->empty()) {
      cerr << "CoordinateBox::BuildFromSmilesToken:empty token(s) " << from_smiles << '\n';
      return 0;
    }
  }

  if (! tokens[0]->numeric_value(_cell_size) || _cell_size <= 0.0) {
    cerr << "CoordinateBox::BuildFromSmilesToken:invalid resolution " << from_smiles << "'\n";
    return 0;
  }

  double dims[2];
  for (int i = 0; i < 2; ++i) {
    if (! tokens[i + 1]->numeric_value(dims[i]) || dims[i] <= 0.0) {
      cerr << "CoordinateBox::BuildFromSmilesToken:invalid cell dimensions " << from_smiles << '\n';
      return 0;
    }
  }

  // Make sure that the dimensions are roughly divisible by _cell_size;
  for (int i = 0; i < 2; ++i) {
    double hopefully_integer = dims[i] / _cell_size;
    if (hopefully_integer - static_cast<int>(hopefully_integer + 0.00001) > 0.1) {
      cerr << "CoordinateBox::BuildFromSmilesToken:invalid cell size " << _cell_size << " for dimension " << dims[i] << '\n';
      return 0;
    }
    _cell[i] = static_cast<int>(hopefully_integer);
  }

  return 1;
}

// Each cell_number maps to a unique position in the box. The first task here is
// to identify the 3 positions defined by `cell_number`.
// Note that Z is unconstrained.
template <typename T> 
std::tuple<T, T, T>
CoordinateBox::CoordinatesAsTuple(int cell_number) const {
  assert(cell_number >= 0);

  const int zcell = cell_number / (_cell[0] * _cell[1]);
  cell_number = cell_number % (_cell[0] * _cell[1]);
  const int ycell = cell_number / _cell[0];
  const int xcell = cell_number % _cell[0];

  const T x = _origin.x() + _cell_size * xcell;
  const T y = _origin.y() + _cell_size * ycell;
  const T z = _origin.z() + _cell_size * zcell;
  return std::tuple<T, T, T>{x, y, z};
}
template std::tuple<float, float, float> CoordinateBox::CoordinatesAsTuple(int cell_number) const;

template <typename T>
Space_Vector<T>
CoordinateBox::CoordinatesAsVector(int cell_number) const {
  auto [x, y, z] = CoordinatesAsTuple<T>(cell_number);
  return Space_Vector<T>(x, y, z);
}

template Space_Vector<float> CoordinateBox::CoordinatesAsVector(int cell_number) const;
template Space_Vector<double> CoordinateBox::CoordinatesAsVector(int cell_number) const;

// Note that we do not check for out of range X or Y values.
// Can be enabled via an ifdef
template < typename T>
int
CoordinateBox::CellNumber(T x, T y, T z) const {
  int xcount = static_cast<int>(x / _cell_size + 0.49999);
  int ycount = static_cast<int>(y / _cell_size + 0.49999);
  int zcount = static_cast<int>(z / _cell_size + 0.49999);

#ifdef CELL_NUMBER_CHECK_IN_BOUNDS
  if (x < 0.0 || xcount >= _cell[0])
    return -1;
  if (y < 0.0 || ycount >= _cell[1])
    return -1;
  if (z < 0.0)
    return -1;
#endif

  if (xcount >= _cell[0] || ycount >= _cell[1]) {
    cerr << "CoordinateBox::CellNumber:value out of range " << x << ',' << y << ',' << z << '\n';
    return -1;
  }

  return (_cell[0] * _cell[1]) * zcount +
          _cell[0] * ycount +
          xcount;
}

template int CoordinateBox::CellNumber(float x, float y, float z) const;
template int CoordinateBox::CellNumber(double x, double y, double z) const;

uint32_t
IntCubeRoot(double d) {
  double result = pow(d, 1.0 / 3.0);
  return static_cast<uint32_t>(result);
}

uint64_t
IntCubeRoot(uint64_t d) {
  const double cube_root = pow(static_cast<double>(d + 0.001), 1.0 / 3.0);
  uint32_t result = static_cast<uint64_t>(cube_root);
  return result;
}

// layer N will have cell numbers out to N^3-1, because the centre cell
// is given number 0.
// layer 0 -> 0
// layer 1 -> [1,26]       // 3^3-1
// layer 2 -> [27,124]     // 5^3-1
// layer 3 -> [125,342]    // 7^3-1
// layer 4 -> [343,728]    // 9^3-1
uint32_t
CellToLayer(uint64_t cell) {
  if (cell <= 26) {
    if (cell == 0) {
      return 0;
    }
    return 1;
  }
  const uint64_t croot = IntCubeRoot(cell);
  return (croot - 1) / 2 + 1;
}

ConcentricBox::ConcentricBox(double res) : _resolution(res) {
  _half_resolution = _resolution * 0.5;
}

//#define DEBUG_CELL_NUMBER_TO_COORDINATES

template <typename T>
Space_Vector<T>
ConcentricBox::CellNumberToCoordinates(uint64_t cell_number) const {
  if (cell_number == 0) {
    return Space_Vector<T>(0, 0, 0);
  }

  const uint32_t layer = CellToLayer(cell_number);
#ifdef DEBUG_CELL_NUMBER_TO_COORDINATES
  cerr << "cell_number " << cell_number << " IntCubeRoot " << layer << '\n';
#endif

  const uint64_t n = 2 * layer + 1;

  uint64_t within_layer = layer == 1 ? cell_number - 1 : cell_number - (n - 2) * (n - 2) * (n - 2);;

  return LayerAndPositionToCoords<T>(layer, within_layer);
}

template <typename T>
Space_Vector<T>
ConcentricBox::LayerAndPositionToCoords(uint32_t layer,
                                        uint64_t within_layer) const {
  const uint64_t n = 2 * layer + 1;
  const uint64_t n2 = n * n;

  // Keep track of the boundaries of each face. `tot` will be the sum
  // of faces passed over so far.
  uint64_t tot = 0;
  uint64_t next_tot = tot + n2;  // Where the next face starts.

#ifdef DEBUG_CELL_NUMBER_TO_COORDINATES
  cerr << "1 cmp " << within_layer << " and " << next_tot << '\n';
#endif
  if (within_layer < next_tot) {
    return _left_face<T>(layer, within_layer - tot);
  }
  tot = next_tot;
  next_tot = tot + n2;
#ifdef DEBUG_CELL_NUMBER_TO_COORDINATES
  cerr << "2 cmp " << within_layer << " and " << next_tot << '\n';
#endif
  if (within_layer < next_tot) {
    return _right_face<T>(layer, within_layer - tot);
  }

  tot = next_tot;
  return _somewhere_along_x<T>(layer, within_layer - tot);
}

// Left face is a full N*N
template <typename T>
Space_Vector<T>
ConcentricBox::_left_face(uint32_t layer, uint64_t within_layer) const {

  const T box_start = - static_cast<T>(layer) * _resolution;
  const T x = box_start;
  // Number of cells on each side of the face
  const int n = 2 * layer + 1;
#ifdef DEBUG_LEFT_FACE
  cerr << "in _left_face layer " << layer << " within " << within_layer << " box_start " << box_start << " n = " << n << '\n';
#endif
  const T y = box_start + static_cast<T>(within_layer / n) * _resolution;
  const T z = box_start + static_cast<T>(within_layer % n) * _resolution;
#ifdef DEBUG_LEFT_FACE
  cerr << "_left_face: layer " << layer << " x " << x << " y " << y << " z " << z << '\n';
#endif
  return Space_Vector<T>(x, y, z);
}

template <typename T>
Space_Vector<T>
ConcentricBox::_right_face(uint32_t layer, uint64_t within_layer) const {
  const T box_start = - static_cast<T>(layer) * _resolution;
  const T x = - box_start;
  // Number of cells on each side of the face
  const int n = 2 * layer + 1;
  const T y = box_start + static_cast<T>(within_layer / n) * _resolution;
  const T z = box_start + static_cast<T>(within_layer % n) * _resolution;
#ifdef DEBUG_RIGHT_FACE
  cerr << "_right_face: layer " << layer << " x " << x << " y " << y << " z " << z << '\n';
#endif
  return Space_Vector<T>(x, y, z);
}

// The cell is not an xmin or xmax, but somewhere in the middle.
// It is a layer of points around the outside.
// First the lower layer in [-z, z].  (lower layer)
// Then ascending (-y, y) along +z  (front face)
// Then at constant Y [-z, z]  (top layer)
// Then (-y, y) along -z (back face).
template <typename T>
Space_Vector<T>
ConcentricBox::_somewhere_along_x(uint32_t layer, uint64_t within_layer) const {
  const T box_start = - static_cast<T>(layer - 1) * _resolution;
  const uint64_t n = 2 * layer + 1;
  const uint64_t cells_per_slice = 4 * n - 4;
  const uint64_t my_slice = within_layer / cells_per_slice;
  const uint64_t within_slice = within_layer % cells_per_slice;
  const T x = box_start + my_slice * _resolution;
#ifdef DEBUG_SOMEWHERE_ALONG_X
  cerr << "_somewhere_along_x cells_per_slice " << cells_per_slice << " my_slice " << my_slice << " x= " << x << " within_slice " << within_slice << " box_start " << box_start << '\n';
#endif
  T y;
  T z;
  if (within_slice < n) {   // Lower face
    y = box_start - _resolution;
    z = box_start - _resolution + within_slice * _resolution;
  } else if (within_slice < (n + n - 2)) {  // Front face
    const uint64_t j = within_slice - n;
    y = box_start + j * _resolution;
    z = box_start - _resolution + (n - 1) * _resolution;
  } else if (within_slice < (n + n - 2 + n)) {  // Top face
    const uint64_t j = within_slice - (n + n - 2);
    y = - box_start + _resolution;
    z = box_start - _resolution + j * _resolution;
  } else {  // back face
    const uint64_t j = within_slice - (n + n - 2 + n);
    y = box_start + j * _resolution;
    z = box_start - _resolution;
  }
  return Space_Vector<T>(x, y, z);
}

template <typename T>
Space_Vector<T>
ConcentricBox::CellToCoordinates(const LayerPosition& layer_pos_in_layer) const {
  return LayerAndPositionToCoords<T>(layer_pos_in_layer.layer, layer_pos_in_layer.position_in_layer);
}

// In order to find the layer, we need to start with the largest
// absolute value of a coordinate.
// Return a tuple of the max coordinate value, and the index of
// that max value.
enum Dimension {
  x, y, z
};

template <typename T>
std::tuple<T, Dimension>
MaxAbsCoordinate(const Space_Vector<T>& coords) {
  T result = abs(coords.x());
  Dimension dim = x;
  if (abs(coords.y()) > result) {
    result = abs(coords.y());
    dim = y;
  }
  if (abs(coords.z()) > result) {
    result = abs(coords.z());
    dim = z;
  }
  return {result, dim};
}

template <typename T>
uint32_t
ConcentricBox::LayerNumber(const Space_Vector<T>& coords) const {
  const auto [max_coord, max_dim] = MaxAbsCoordinate(coords);
  return static_cast<uint32_t>((max_coord + _half_resolution) / _resolution);
}

//#define DEBUG_CELL_NUMBER

template <typename T>
uint64_t
ConcentricBox::CellNumber(const Space_Vector<T>& coords) const {
  const uint32_t layer = LayerNumber(coords);
#ifdef DEBUG_CELL_NUMBER
  const auto [max_coord, max_dim] = MaxAbsCoordinate(coords);
  cerr << "CellNumber::max_coord " << max_coord << " layer " << layer << '\n';
#endif
  if (layer == 0) {
    return 0;
  }

  const uint64_t n = 2 * layer + 1;

  // Deliberately promote to double. Not sure if it will make any difference.
  const double box_min = - static_cast<double>(n) * _half_resolution;
#ifdef DEBUG_CELL_NUMBER
  cerr << coords <<  " max_coord " << max_coord << " along " << max_dim << " layer " << layer << " n = " << n << " box_min " << box_min << '\n';
#endif

  const uint32_t x = static_cast<uint32_t>((coords.x() - box_min) / _resolution);
  const uint32_t y = static_cast<uint32_t>((coords.y() - box_min) / _resolution);
  const uint32_t z = static_cast<uint32_t>((coords.z() - box_min) / _resolution);

#ifdef DEBUG_CELL_NUMBER
  cerr << "Int coords " << x << " " << y << " " << z << " layer " << layer << '\n';
#endif

  return CellNumber(coords, layer, x, y, z);
}

//Given `coords` which are known to be in `layer`, and integer equivalent
// coordinates (x, y, z), return the corresponding cell number.
// Logic branches depending on whether `coords` is on the left face, the
// right face, or somewhere in between.
// On the lower face there are n points.
// On the front face there are (n-2) points.
// On the top face there are n points.
// On the back face there are (n-2) points.
template <typename T>
uint64_t
ConcentricBox::CellNumber(const Space_Vector<T>& coords,
                          uint32_t layer,
                          uint32_t x, uint32_t y, uint32_t z) const {
  // Number of cells along each edge of the box in this layer.
  const uint64_t n = 2 * layer + 1;

  // The inner cell.
  const uint64_t n2 = n - 2;
  uint64_t result = layer == 1 ? 1 : n2 * n2 * n2;

#ifdef DEBUG_CELL_NUMBER
  cerr << "ConcentricBox::CellNumberX: layer " << layer << " x " << x << " y " << y << " z " << z << " n = " << n << " start " << result << '\n';
#endif

  const uint64_t position = PositionInLayer(coords, layer, x, y, z);

  return result + position;
}

template <typename T>
LayerPosition
ConcentricBox::Position(const Space_Vector<T>& coords) const {
  const uint32_t layer = LayerNumber(coords);
#ifdef DEBUG_CELL_NUMBER
  const auto [max_coord, max_dim] = MaxAbsCoordinate(coords);
  cerr << "CellNumber::max_coord " << max_coord << " layer " << layer << '\n';
#endif
  if (layer == 0) {
    LayerPosition result;
    return result;
  }

  const uint64_t n = 2 * layer + 1;

  // Deliberately promote to double. Not sure if it will make any difference.
  const double box_min = - static_cast<double>(n) * _half_resolution;
#ifdef DEBUG_CELL_NUMBER
  cerr << coords <<  " max_coord " << max_coord << " along " << max_dim << " layer " << layer << " n = " << n << " box_min " << box_min << '\n';
#endif

  const uint64_t x = static_cast<uint64_t>((coords.x() - box_min) / _resolution);
  const uint64_t y = static_cast<uint64_t>((coords.y() - box_min) / _resolution);
  const uint64_t z = static_cast<uint64_t>((coords.z() - box_min) / _resolution);

#ifdef DEBUG_CELL_NUMBER
  cerr << "Int coords " << x << " " << y << " " << z << " layer " << layer << '\n';
#endif

  const uint64_t position = PositionInLayer(coords, layer, x, y, z);

  LayerPosition result(layer, position);

  return result;
}

template <typename T>
uint64_t
ConcentricBox::PositionInLayer(const Space_Vector<T>& coords,
                        uint32_t layer,
                        uint64_t x, uint64_t y, uint64_t z) const {
  // Number of cells along each edge of the box in this layer.
  const uint64_t n = 2 * layer + 1;

  // Perform calculations in 64 bit precision to avoid overflow.
  uint64_t position = 0;

#ifdef DEBUG_CELL_NUMBER
  cerr << "ConcentricBox::CellNumberX: layer " << layer << " x " << x << " y " << y << " z " << z << " n = " << n << " start " << result << '\n';
#endif

  if (x == 0) {  // Left Face
    position += y * n + z;
  } else if (x == n - 1) {  // Right Face
    position += n * n;
    position += y * n + z;
  } else {  // Somewhere in the middle. Skip over layers.
    position += 2 * n * n;  // Skip over both faces (n*n + n*n)
    position += (x - 1) * (4 * n - 4);  // Skip over previous x layers (if any)
    if (y == 0) {
      position += z;
    } else if (y == n - 1) {
      position += n + (n - 2);  // Skip over right wall
      position += z;
    } else {  // In the middle;
      if (coords.z() < 0.0) {
        position += (n + (n - 2) + n);
      } else {
        position += n;
      }
#ifdef DEBUG_CELL_NUMBER
      cerr << "CellNumberX " << coords << " layer " << layer << " x " << x << " y " << y << " n " << n << " z " << z << '\n';
#endif
      position += (y - 1);
    }
  }

  return position;
}

constexpr char kComponentSeparator = ':';

int
FromString(const const_IWSubstring& s, LayerPosition& destination) {
  destination.layer = 0;
  destination.position_in_layer = 0;
  const int nchars = s.length();
  int i = 0;
  for ( ; i < nchars; ++i) {
    if (s[i] == kComponentSeparator) {
      i++;
      break;
    }
    int digit = s[i] - '0';
    if (digit >= 0 && digit <= 9) {
      destination.layer = destination.layer * 10 + digit;
    } else {
      return 0;
    }
  }

  if (i <= 1) {
    return 0;
  }

  bool got_position_characters = false;
  for ( ; i < nchars; ++i) {
    int digit = s[i] - '0';
    if (digit >= 0 && digit <= 9) {
      destination.position_in_layer = destination.position_in_layer * 10 + digit;
      got_position_characters = true;
    } else {
      if (got_position_characters) {
        return i;
      } else {
        return 0;
      }
    }
  }

  if (got_position_characters) {
    return nchars;
  }

  return 0;
}

IWString &
operator << (IWString & output, const LayerPosition& layer_position) {
  output << layer_position.layer << kComponentSeparator << layer_position.position_in_layer;
  return output;
}

template uint64_t ConcentricBox::CellNumber(const Space_Vector<float>& coords) const;
template uint64_t ConcentricBox::CellNumber(const Space_Vector<double>& coords) const;

template Space_Vector<float>
ConcentricBox::CellNumberToCoordinates(uint64_t cell_number) const;
template Space_Vector<double>
ConcentricBox::CellNumberToCoordinates(uint64_t cell_number) const;

template Space_Vector<float>
ConcentricBox::CellToCoordinates(const LayerPosition& layer_pos_in_layer) const;
template Space_Vector<double>
ConcentricBox::CellToCoordinates(const LayerPosition& layer_pos_in_layer) const;

template LayerPosition ConcentricBox::Position<float>(Space_Vector<float> const&) const;
template LayerPosition ConcentricBox::Position<double>(Space_Vector<double> const&) const;
}  // namespace coordinate_box
