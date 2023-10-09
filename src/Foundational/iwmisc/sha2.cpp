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
#define SHA2_CPP
#include "sha2.h"
#include <cstring>
#include <algorithm>


namespace sha2 {

const State<uint32_t> h0_224 = {
    0xc1059ed8, 0x367cd507, 0x3070dd17, 0xf70e5939,
    0xffc00b31, 0x68581511, 0x64f98fa7, 0xbefa4fa4
};
const State<uint32_t> h0_256 = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
};

const State<uint64_t> h0_384 = {
    0xcbbb9d5dc1059ed8, 0x629a292a367cd507,
    0x9159015a3070dd17, 0x152fecd8f70e5939,
    0x67332667ffc00b31, 0x8eb44a8768581511,
    0xdb0c2e0d64f98fa7, 0x47b5481dbefa4fa4
};

const State<uint64_t> h0_512 = {
    0x6a09e667f3bcc908, 0xbb67ae8584caa73b,
    0x3c6ef372fe94f82b, 0xa54ff53a5f1d36f1,
    0x510e527fade682d1, 0x9b05688c2b3e6c1f,
    0x1f83d9abfb41bd6b, 0x5be0cd19137e2179
};

const State<uint64_t> h0_512_224 = {
    0x8c3d37c819544da2, 0x73e1996689dcd4d6,
    0x1dfab7ae32ff9c82, 0x679dd514582f9fcf,
    0x0f6d2b697bd44da8, 0x77e36f7304c48942,
    0x3f9d85a86a1d36c8, 0x1112e6ad91d692a1
};

const State<uint64_t> h0_512_256 = {
    0x22312194fc2bf72c, 0x9f555fa3c84c64c2,
    0x2393b86b6f53b151, 0x963877195940eabd,
    0x96283ee2a88effe3, 0xbe5e1e2553863992,
    0x2b0199fc2c85b8aa, 0x0eb72ddc81c52ca2
};

} // namespace sha2

namespace {

constexpr array<uint32_t, 64> K_256 = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

constexpr array<uint64_t, 80> K_512 = {
    0x428a2f98d728ae22, 0x7137449123ef65cd,
    0xb5c0fbcfec4d3b2f, 0xe9b5dba58189dbbc,
    0x3956c25bf348b538, 0x59f111f1b605d019,
    0x923f82a4af194f9b, 0xab1c5ed5da6d8118,
    0xd807aa98a3030242, 0x12835b0145706fbe,
    0x243185be4ee4b28c, 0x550c7dc3d5ffb4e2,
    0x72be5d74f27b896f, 0x80deb1fe3b1696b1,
    0x9bdc06a725c71235, 0xc19bf174cf692694,
    0xe49b69c19ef14ad2, 0xefbe4786384f25e3,
    0x0fc19dc68b8cd5b5, 0x240ca1cc77ac9c65,
    0x2de92c6f592b0275, 0x4a7484aa6ea6e483,
    0x5cb0a9dcbd41fbd4, 0x76f988da831153b5,
    0x983e5152ee66dfab, 0xa831c66d2db43210,
    0xb00327c898fb213f, 0xbf597fc7beef0ee4,
    0xc6e00bf33da88fc2, 0xd5a79147930aa725,
    0x06ca6351e003826f, 0x142929670a0e6e70,
    0x27b70a8546d22ffc, 0x2e1b21385c26c926,
    0x4d2c6dfc5ac42aed, 0x53380d139d95b3df,
    0x650a73548baf63de, 0x766a0abb3c77b2a8,
    0x81c2c92e47edaee6, 0x92722c851482353b,
    0xa2bfe8a14cf10364, 0xa81a664bbc423001,
    0xc24b8b70d0f89791, 0xc76c51a30654be30,
    0xd192e819d6ef5218, 0xd69906245565a910,
    0xf40e35855771202a, 0x106aa07032bbd1b8,
    0x19a4c116b8d2d0c8, 0x1e376c085141ab53,
    0x2748774cdf8eeb99, 0x34b0bcb5e19b48a8,
    0x391c0cb3c5c95a63, 0x4ed8aa4ae3418acb,
    0x5b9cca4f7763e373, 0x682e6ff3d6b2b8a3,
    0x748f82ee5defb2fc, 0x78a5636f43172f60,
    0x84c87814a1f0ab72, 0x8cc702081a6439ec,
    0x90befffa23631e28, 0xa4506cebde82bde9,
    0xbef9a3f7b2c67915, 0xc67178f2e372532b,
    0xca273eceea26619c, 0xd186b8c721c0c207,
    0xeada7dd6cde0eb1e, 0xf57d4f7fee6ed178,
    0x06f067aa72176fba, 0x0a637dc5a2c898a6,
    0x113f9804bef90dae, 0x1b710b35131c471b,
    0x28db77f523047d84, 0x32caab7b40c72493,
    0x3c9ebe0a15c9bebc, 0x431d67c49c100d4c,
    0x4cc5d4becb3e42b6, 0x597f299cfc657e2a,
    0x5fcb6fab3ad6faec, 0x6c44198c4a475817
};

template<typename Uword>
inline Uword shr(Uword x, size_t n) {return x >> n;}

template<typename Uword>
inline Uword rotr(Uword x, size_t n) {
    return x >> n | x << (8 * sizeof(Uword) - n);
}

template<typename Word>
inline Word ch(Word x, Word y, Word z) {
    return (x & y) ^ (~x & z);
}

template<typename Word>
inline Word maj(Word x, Word y, Word z) {
    return (x & y) ^ (x & z) ^ (y & z);
}

struct Sigma_256 {
    using Word = uint32_t;
    static Word big_0(Word x) {
        return rotr(x,  2) ^ rotr(x, 13) ^ rotr(x, 22);
    }
    static Word big_1(Word x) {
        return rotr(x,  6) ^ rotr(x, 11) ^ rotr(x, 25);
    }
    static Word little_0(Word x) {
        return rotr(x,  7) ^ rotr(x, 18) ^  shr(x,  3);
    }
    static Word little_1(Word x) {
        return rotr(x, 17) ^ rotr(x, 19) ^  shr(x, 10);
    }
};

struct Sigma_512 {
    using Word = uint64_t;
    static Word big_0(Word x) {
        return rotr(x, 28) ^ rotr(x, 34) ^ rotr(x, 39);
    }
    static Word big_1(Word x) {
        return rotr(x, 14) ^ rotr(x, 18) ^ rotr(x, 41);
    }
    static Word little_0(Word x) {
        return rotr(x,  1) ^ rotr(x,  8) ^  shr(x,  7);
    }
    static Word little_1(Word x) {
        return rotr(x, 19) ^ rotr(x, 61) ^  shr(x,  6);
    }
};

template<typename Sigma, typename Array>
void do_process_block(sha2::State<typename Array::value_type>& H,
                      const Hash::Byte* block, const Array& K)
{
    using Word = typename Array::value_type;
    constexpr auto rounds_count = K.size();

    Word w[rounds_count];

    sha2::copy_bigendian(reinterpret_cast<const Word*>(block), 16, w);

    for (size_t t = 16; t < rounds_count; t++) {
        w[t] = Sigma::little_1(w[t-2])  + w[t-7]
             + Sigma::little_0(w[t-15]) + w[t-16];
    }

    Word a = H[0],
         b = H[1],
         c = H[2],
         d = H[3],
         e = H[4],
         f = H[5],
         g = H[6],
         h = H[7];

    for (size_t t = 0; t < rounds_count; t++) {
        Word T1 = h + Sigma::big_1(e) + ch(e, f, g) + K[t] + w[t];
        Word T2 = Sigma::big_0(a) + maj(a, b, c);
        h = g;
        g = f;
        f = e;
        e = d + T1;
        d = c;
        c = b;
        b = a;
        a = T1 + T2;
    }

    H[0] += a;
    H[1] += b;
    H[2] += c;
    H[3] += d;
    H[4] += e;
    H[5] += f;
    H[6] += g;
    H[7] += h;
}

template<typename Word>
void do_padding(sha2::State<Word>& h, sha2::Block_bytes<Word>& buffer,
                uint64_t byte_count)
{
    size_t i = byte_count & (buffer.size() - 1);
    buffer[i++] = 0x80;
    if (i >= buffer.size() - 2 * sizeof(Word)) {
        fill(buffer.begin() + i, buffer.end(), 0);
        sha2::process_block(h, buffer.data());
        i = 0;
    }
    auto end_zeros = buffer.data() + (buffer.size() - sizeof(uint64_t));
    fill(buffer.data() + i, end_zeros, 0);
    uint64_t len = byte_count * 8;
    sha2::copy_bigendian(&len, 1, reinterpret_cast<uint64_t*>(end_zeros));
}

} // namespace

namespace sha2 {

void process_block(State<uint32_t>& h, const Hash::Byte* block)
{
    do_process_block<Sigma_256>(h, block, K_256);
}

void process_block(State<uint64_t>& h, const Hash::Byte* block)
{
    do_process_block<Sigma_512>(h, block, K_512);
}

void padding(State<uint32_t>& h, Block_bytes<uint32_t>& buffer,
             uint64_t byte_count)
{
    do_padding(h, buffer, byte_count);
}

void padding(State<uint64_t>& h, Block_bytes<uint64_t>& buffer,
             uint64_t byte_count)
{
    do_padding(h, buffer, byte_count);
}

} // namespace sha2
