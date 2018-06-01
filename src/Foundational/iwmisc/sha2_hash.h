/*
 * Copyright (C) 2017 Jerome Migne
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef HASH_H
#define HASH_H
#include <cstddef>
#include <cstdint>
#include <vector>

class Hash {
  public:
    using Byte = std::uint8_t;

    virtual ~Hash() {}

    virtual void reset() = 0;
    virtual void update(const void*, std::size_t) = 0;
    virtual void terminate(std::vector<Byte>& out) = 0;

    void compute(const void* in, std::size_t size, std::vector<Byte>& out) {
        reset();
        update(in, size);
        terminate(out);
    }
};

#endif
