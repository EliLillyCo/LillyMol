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
#ifndef SHA2_H
#define SHA2_H
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <array>
#include "sha2_hash.h"

#ifdef SHA2_CPP
#define SHA2_TEMPLATE template
#else
#define SHA2_TEMPLATE extern template
#endif

namespace sha2 {
constexpr bool cpu_is_big_endian = false;

template<typename Word> using State = std::array<Word, 8>;
template<typename Word>
using Block_bytes = std::array<Hash::Byte, 16 * sizeof(Word)>;

extern const State<std::uint32_t> h0_224;
extern const State<std::uint32_t> h0_256;
extern const State<std::uint64_t> h0_384;
extern const State<std::uint64_t> h0_512;
extern const State<std::uint64_t> h0_512_224;
extern const State<std::uint64_t> h0_512_256;

template<typename Word>
void copy_bigendian(const Word* src, std::size_t words_len, Word* dst) {
    if (cpu_is_big_endian) {
        std::copy_n(src, words_len, dst);
        return;
    }
    std::transform(src, src + words_len, dst, [](Word w) {
        for (char *low = reinterpret_cast<char*>(&w),
                  *high = low + (sizeof(w) - 1);
             low < high; low++, high--)
        {
            std::swap(*low, *high);
        }
        return w;
    });
}

void process_block(State<std::uint32_t>&, const Hash::Byte* block);
void process_block(State<std::uint64_t>&, const Hash::Byte* block);
void padding(State<std::uint32_t>&, Block_bytes<std::uint32_t>& buffer,
             std::uint64_t byte_count);
void padding(State<std::uint64_t>&, Block_bytes<std::uint64_t>& buffer,
             std::uint64_t byte_count);
}

template<typename Word>
struct Sha2_ctx {
    using State = sha2::State<Word>;
    using Block = sha2::Block_bytes<Word>;

    State h;
    Block buffer;
    std::uint64_t byte_count;

    void init(const State& h0) {
        h = h0;
        byte_count = 0;
    }

    void update(const void* data, std::size_t len) {
        auto* p = static_cast<const Hash::Byte*>(data);
        auto* p_end = p + len;
        if (std::size_t i = byte_count & (buffer.size() - 1)) {
            auto n = std::min(buffer.size() - i, len);
            std::copy_n(p, n, buffer.begin() + i);
            p += n;
            if (i + n == buffer.size()) {
                sha2::process_block(h, buffer.data());
            }
        }
        for (auto* q = p; (q += buffer.size()) <= p_end; p = q) {
            sha2::process_block(h, p);
        }
        std::copy(p, p_end, buffer.begin());
        byte_count += len;
    }

    void terminate() {
        sha2::padding(h, buffer, byte_count);
        sha2::process_block(h, buffer.data());
    }
};

template<typename Word, std::size_t digest_bits_len,
         const sha2::State<Word>& h0>
class Sha2 : public Hash {
    Sha2_ctx<Word> ctx;
  public:
    static constexpr std::size_t digest_size = digest_bits_len / 8;

    virtual void reset() override {ctx.init(h0);}

    Sha2() {reset();}

    virtual void update(const void* data, std::size_t len) override {
        ctx.update(data, len);
    }

    virtual void terminate(std::vector<Byte>& out) override {
        ctx.terminate();
        constexpr std::size_t words_len = (digest_size + sizeof(Word) / 2)
                                        / sizeof(Word);
        out.resize(words_len * sizeof(Word));
        sha2::copy_bigendian(ctx.h.data(), words_len,
                             reinterpret_cast<Word*>(out.data()));
        out.resize(digest_size);
    }
};

SHA2_TEMPLATE class Sha2_ctx<std::uint32_t>;
SHA2_TEMPLATE class Sha2_ctx<std::uint64_t>;
SHA2_TEMPLATE class Sha2<std::uint32_t, 224, sha2::h0_224>;
SHA2_TEMPLATE class Sha2<std::uint32_t, 256, sha2::h0_256>;
SHA2_TEMPLATE class Sha2<std::uint64_t, 384, sha2::h0_384>;
SHA2_TEMPLATE class Sha2<std::uint64_t, 512, sha2::h0_512>;
SHA2_TEMPLATE class Sha2<std::uint64_t, 224, sha2::h0_512_224>;
SHA2_TEMPLATE class Sha2<std::uint64_t, 256, sha2::h0_512_256>;
#undef SHA2_TEMPLATE

using Sha224     = Sha2<std::uint32_t, 224, sha2::h0_224>;
using Sha256     = Sha2<std::uint32_t, 256, sha2::h0_256>;
using Sha384     = Sha2<std::uint64_t, 384, sha2::h0_384>;
using Sha512     = Sha2<std::uint64_t, 512, sha2::h0_512>;
using Sha512_224 = Sha2<std::uint64_t, 224, sha2::h0_512_224>;
using Sha512_256 = Sha2<std::uint64_t, 256, sha2::h0_512_256>;

#endif
