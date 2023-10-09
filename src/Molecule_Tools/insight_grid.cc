#include <iostream>
#include <tuple>

#include "Foundational/iwmisc/misc.h"

#include "Molecule_Tools/insight_grid.h"

namespace insight_grid {

using std::cerr;

InsightGrid::InsightGrid() {
  _xlen = 0;
  _ylen = 0;
  _zlen = 0;

  _dx = 0.0;
  _dy = 0.0;
  _dz = 0.0;

  _x = 0;
  _y = 0;
  _z = 0;

  _data = nullptr;
  _npoints = 0;

  _isotope_multiplier = 0.0f;
}

InsightGrid::~InsightGrid() {
  if (_data != nullptr) {
    delete [] _data;
  }
}

int
InsightGrid::Build(const char * fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "InsightGrid::Build:cannot open '" << fname << "'\n";
    return 0;
  }

  return Build(input);
}

int
InsightGrid::Build(iwstring_data_source& input) {
  input.set_strip_trailing_blanks(1);

  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "InsightGrid::Build:empty file\n";
    return 0;
  }

  if (buffer != "GRD InsightII ASCII grid map") {
    cerr << "InsightGrid::Build:invalid header record\n";
    cerr << buffer << '\n';
    return 0;
  }

  if (! input.next_record(buffer)) {
    cerr << "InsightGrid::Build:premature EOF\n";
    return 0;
  }

  if (! input.next_record(buffer)) {
    cerr << "InsightGrid::Build:premature EOF\n";
    return 0;
  }

  if (! GetDimensions(buffer)) {
    cerr << "InsightGrid::Build:cannot parse dimensions\n";
    cerr << buffer << '\n';
    return 0;
  }

  if (! input.next_record(buffer)) {
    cerr << "InsightGrid::Build:premature EOF\n";
    return 0;
  }

  if (! GetResolution(buffer)) {
    cerr << "InsightGrid::Build:cannot parse dimensions\n";
    cerr << buffer << '\n';
  }

  _npoints = (_x + 1) * (_y + 1) * (_z + 1);
  _data = new_float(_npoints);

  if (! input.next_record(buffer)) {
    cerr << "InsightGrid::Build:premature EOF\n";
    return 0;
  }

  if (! GetGridOffsets(buffer)) {
    cerr << "InsightGrid::Build:cannot interpret grid offsets\n";
    cerr << buffer << '\n';
    return 0;
  }

  input.set_strip_leading_blanks(1);

  int nonzero_values = 0;
  Accumulator<double> acc;
  for (unsigned int i = 0; i < _npoints; ++i) {
    if (! input.next_record(buffer)) {
      cerr << "InsightGrid::Build:premature EOF\n";
      return 0;
    }

    if (! buffer.numeric_value(_data[i])) {
      cerr << "InsightGrid::BUild:invalid numeric '" << buffer << "'\n";
      return 0;
    }
    acc.extra(_data[i]);
    if (_data[i] > 0.0f) {
      ++nonzero_values;
    }
  }

  cerr << "Read " << _npoints << " data records, found " << nonzero_values << " nonzero values\n";
  cerr << "Values btw " << acc.minval() << " and " << acc.maxval() << " ave " << acc.average() << '\n';
  DebugPrint(cerr);
#ifdef ECHO_NONZERO_VALUES
  for (int i = 0; i < n; ++i) {
    if (_data[i] > 0.0f) {
      cerr << " data " << i << " value " << _data[i] << '\n';
    }
  }
#endif

  return 1;
}

int
InsightGrid::DebugPrint(std::ostream& output) const {
  output << "Grid containing " << _npoints << " points\n";
  output << "Ranges x [" << _xmin << ',' << _xmax << "] y ["
                         << _ymin << ',' << _ymax << "] z ["
                         << _zmin << ',' << _zmax << "]\n";
  output << "Sizes " << _x << ',' << _y << ',' << _z << '\n';
  output << "Granularity " << _dx << ',' << _dy << ',' << _dz << '\n';
  return 1;
}

// Parse something that looks like...
// 65.800    63.700    73.500    90.000    90.000    90.000
int
InsightGrid::GetDimensions(const const_IWSubstring& buffer) {
  int i = 0;
  const_IWSubstring token;
  for (int col = 0; buffer.nextword(token, i); ++col) {
    float value;
    if (! token.numeric_value(value)) {
      cerr << "InsightGrid::GetDimensions:invalid dimension '" << token << "'\n";
      return 0;
    }

    if (col == 0) {
      _xlen = value;
    } else if (col == 1) {
      _ylen = value;
    } else if (col == 2) {
      _zlen = value;
    }
  }

  return 1;
}

// Parse the record that looks like
//       94        91       105
int
InsightGrid::GetResolution(const const_IWSubstring& buffer) {
  int i = 0;
  const_IWSubstring token;
  for (int col = 0; buffer.nextword(token, i); ++col) {
    int value;
    if (! token.numeric_value(value) || value <= 0) {
      cerr << "InsightGrid::GetResolution:invalid value '" << token  << "'\n";
      return 0;
    }

    if (col == 0) {
      _x = value;
      _dx = _xlen / static_cast<double>(value);
    } else  if (col == 1) {
      _y = value;
      _dy = _ylen / static_cast<double>(value);
    } else  if (col == 2) {
      _z = value;
      _dz = _zlen / static_cast<double>(value);
    }
  }

  return 1;
}

// The grid offsets
//     1       -98        -4       -47        44       415       520
int
InsightGrid::GetGridOffsets(const const_IWSubstring& buffer) {
  int i = 0;
  const_IWSubstring token;
  for (int col = 0; buffer.nextword(token, i); ++col) {
    int value;
    if (! token.numeric_value(value)) {
      cerr << "InsightGrid::GetGridOffsets:invalid offset " << token << '\n';
      return 0;
    }
    if (col == 1)  {
      _xmin = value * _dx;
    } else if (col == 2) {
      _xmax = _xmin + (_x + 1) * _dx;
    } else if (col == 3) {
      _ymin = value * _dy;
    } else if (col == 4) {
      _ymax = _ymin + (_y + 1) * _dy;
    } else if (col == 5) {
      _zmin = value * _dz;
    } else if (col == 6) {
      _zmax = _zmin + (_z + 1) * _dz;
    }
  }

  return 1;
}

std::tuple<int, int>
GridPoints(float value,
           float distance,
           float minval,
           float dx,
           float maxval) {
  float lower = value - distance;
  float upper = value + distance;
  if (lower < minval) {
    lower = minval;
  } 
  if (upper > maxval) {
    upper = maxval;
  }

  int istart = static_cast<int>((lower - minval) / dx);
  int istop = static_cast<int>((upper - minval) / dx);
  return std::make_tuple(istart, istop);
}

//#define DEBUG_SUM

// The hardest part here is indexing into the _data array.
// By observation, the X value varies fastest.
int
InsightGrid::Sum(float x, float y, float z, float distance,
                Accumulator<double>& acc) const {
#ifdef DEBUG_SUM
  cerr << "Sum x " << x << " y " << y << " z " << z << '\n';
#endif
  if (x + distance < _xmin || y + distance < _ymin || z + distance < _zmin) {
#ifdef DEBUG_SUM
    cerr << "Out of range <min\n";
#endif
    return 0;
  }
  if (x - distance > _xmax || y - distance > _ymax || z - distance > _zmax) {
#ifdef DEBUG_SUM
    cerr << "Out of range >max\n";
#endif
    return 0;
  }

  auto [xstart, xstop] = GridPoints(x, distance, _xmin, _dx, _xmax);
  auto [ystart, ystop] = GridPoints(y, distance, _ymin, _dy, _ymax);
  auto [zstart, zstop] = GridPoints(z, distance, _zmin, _dz, _zmax);

#ifdef DEBUG_SUM
  cerr << "xstart " << xstart << " xstop " << xstop << " ystart " << ystart << " ystop " << ystop << " zstart " << zstart << " zstop " << zstop << '\n';
#endif

  // Avoid taking square roots.
  const double distance_squared = distance * distance;

  int rc = 0;

  for (int i = zstart; i < zstop; ++i) {
    double zz = _zmin + i * _dz;
    uint32_t zoffset = i * ((_y + 1) * (_x + 1));
//  cerr << " zz " << zz << " zoffset " << zoffset << '\n';
    for (int j = ystart; j < ystop; ++j) {
      uint32_t yoffset = zoffset + j * (_x + 1);
      double yy = _ymin + j * _dy;
      for (int k = xstart; k < xstop; ++k) {
        double xx = _xmin + k * _dx;
        double dist = ((x - xx) * (x - xx) +
                       (y - yy) * (y - yy) +
                       (z - zz) * (z - zz));
//      cerr << "at " << (yoffset + k) << " value " << _data[yoffset + k] << " distance " << sqrt(dist) << " cmp " << distance << '\n';
        if (dist > distance_squared) {
          continue;
        }
        acc.extra(_data[yoffset + k]);
        if (_data[yoffset + k] > 0.0f) {
          rc = 1;
        }
      }
    }
  }

  return rc;
}

int
InsightGrid::IndexInData(float x, float y, float z) const {
  if (x < _xmin || x > _xmax) {
    return -1;
  }

  if (y < _ymin || y > _ymax) {
    return -1;
  }

  if (z < _zmin || z > _zmax) {
    return -1;
  }

  int result = 0;

  int iz = static_cast<int>((z - _zmin) / _dz);
  result += iz * (_x + 1) * (_y + 1);
  int iy = static_cast<int>((y - _ymin) / _dy);
  result += iy * (_x + 1);
  int ix = static_cast<int>((x - _xmin) / _dx);
  return result + ix;

}

int
InsightGrid::IndicesInRange(float x, float y, float z, float distance, int* in_range) const {
  if (x + distance < _xmin || y + distance < _ymin || z + distance < _zmin) {
    return 0;
  }
  if (x - distance > _xmax || y - distance > _ymax || z - distance > _zmax) {
    return 0;
  }

  auto [xstart, xstop] = GridPoints(x, distance, _xmin, _dx, _xmax);
  auto [ystart, ystop] = GridPoints(y, distance, _ymin, _dy, _ymax);
  auto [zstart, zstop] = GridPoints(z, distance, _zmin, _dz, _zmax);

#ifdef DEBUGINDICESINRANGE
  cerr << "coords " << x << ',' << y << ',' << z << '\n';
  cerr << "xstart " << xstart << " xstop " << xstop << " ystart " << ystart << " ystop " << ystop << " zstart " << zstart << " zstop " << zstop << '\n';
#endif
  // Avoid taking square roots.
  const double distance_squared = distance * distance;

  int rc = 0;

  for (int i = zstart; i < zstop; ++i) {
    double zz = _zmin + i * _dz;
    uint32_t zoffset = i * ((_y + 1) * (_x + 1));
    for (int j = ystart; j < ystop; ++j) {
      uint32_t yoffset = zoffset + j * (_x + 1);
      double yy = _ymin + j * _dy;
      for (int k = xstart; k < xstop; ++k) {
        double xx = _xmin + k * _dx;
        double dist = ((x - xx) * (x - xx) +
                       (y - yy) * (y - yy) +
                       (z - zz) * (z - zz));
        if (dist > distance_squared) {
          continue;
        }
        if (_data[yoffset + k] == 0.0f) {
          continue;
        }
        in_range[yoffset + k] = 1;
        ++rc;
      }
    }
  }

  return rc;
}

void
InsightGrid::Score(const int * ndx, Accumulator<double>& acc) const {
  for (uint32_t i = 0; i < _npoints; ++i) {
    if (ndx[i] == 0) {
      continue;
    }
    acc.extra(_data[i]);
  }
}

int
InsightGrid::MakeMolecule(const Element* element,
                          float threshold,
                          const Molecule& m,
                          Molecule& grid_mol) const {
  std::unique_ptr<int[]> indices_to_use(new_int(_npoints));
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m.atom(i);
    IndicesInRange(a.x(), a.y(), a.z(), 2.0 * threshold, indices_to_use.get());
  }

  for (uint32_t i = 0; i < _npoints; ++i) {
    indices_to_use[i] = ! indices_to_use[i];
  }

  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m.atom(i);
    IndicesInRange(a.x(), a.y(), a.z(), threshold, indices_to_use.get());
  }

  for (uint32_t i = 0; i < _npoints; ++i) {
    indices_to_use[i] = ! indices_to_use[i];
  }

#ifdef DEBUG_MAKE_MOLECULE
  int number_in_range = 0;
  for (int i = 0; i < _npoints; ++i) {
    if (indices_to_use[i]) {
      ++number_in_range;
    }
  }
  cerr << "There are " << number_in_range << " indices far from mol\n";
#endif

  grid_mol.resize(_npoints + m.natoms());

  for (uint32_t i = 0; i <= _z; ++i) {
    const double z = _zmin + i * _dz;
    const uint32_t zoffset = i * (_y + 1) * (_x + 1);
    for (uint32_t j = 0; j <= _y; ++j) {
      const double y = _ymin + j * _dy;
      const uint32_t yoffset = zoffset + j * (_x + 1);
      for (uint32_t k = 0; k <= _x; ++k) {
        if (! indices_to_use[yoffset + k]) {
          continue;
        }

        const double x = _xmin + k * _dx;
        int atom_number = grid_mol.natoms();
        grid_mol.add(element);
        grid_mol.setxyz(atom_number, x, y, z);
        if (_isotope_multiplier > 0) {
          grid_mol.set_isotope(atom_number,
                        static_cast<int>(_isotope_multiplier * _data[yoffset + k]));
        }
      }
    }
  }

  return m.natoms();
}

}  // namespace insight_grid
