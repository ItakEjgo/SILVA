#include <parlay/sequence.h>
#pragma once
#pragma once

#include <bits/stdc++.h>
#include <parlay/internal/binary_search.h>
#include "../hilbert.h"

namespace geobase
{
    using namespace std;
    using FT = double;
    // using FT = float;
    constexpr FT FT_INF_MIN = numeric_limits<FT>::min();
    constexpr FT FT_INF_MAX = numeric_limits<FT>::max();
    constexpr FT FT_EPS = numeric_limits<FT>::epsilon();

    struct break_down
    {
        FT sort_time = 0;
        FT leaf_time = 0;
        FT inte_time = 0;
        FT split_time = 0;
        FT slice_time = 0;
        void clear()
        {
            sort_time = 0;
            leaf_time = 0;
            inte_time = 0;
            split_time = 0;
            slice_time = 0;
        }
    };

    inline int dcmp(const FT &x)
    {
        if (fabs(x) < FT_EPS)
            return 0;
        return x < 0 ? -1 : 1;
    }

    bool less_msb(unsigned int x, unsigned int y)
    {
        return x < y && x < (x ^ y);
    }

    struct Point
    {
        size_t id;
        FT x, y;
        unsigned long long morton_id;
        Point() {}
        Point(FT _x, FT _y) : x(_x), y(_y) {}
        Point(size_t _id, FT _x, FT _y) : id(_id), x(_x), y(_y)
        {
            // morton_id = mortonIndex();
            // morton_id = interleave_bits();
        }

        bool operator==(const Point &p) const
        {
            return !(id - p.id) && !dcmp(x - p.x) && !dcmp(y - p.y);
        }

        friend std::ostream &operator<<(std::ostream &os, const Point &p)
        {
            os << fixed << setprecision(6) << p.id << ": (" << p.x << ", " << p.y << ")";
            // os << "(" << p.x << ", " << p.y << ")";
            return os;
        }

        // bool operator < (const Point &b) const{
        //     return (morton_id < b.morton_id) ||
        //         (morton_id == b.morton_id && id < b.id);
        // }

        /* return Z value of this point */
        unsigned long long interleave_bits() const
        {
            // Pun the x and y coordinates as integers: Just re-interpret the bits.
            //
            auto ix = static_cast<unsigned int>(x);
            auto iy = static_cast<unsigned int>(y);
            // cout << ix << ", " << iy << endl;
            // cout << bitset<32>(ix) << endl;
            // cout << bitset<32>(iy) << endl;

            auto ret = 0ull;
            for (auto i = 0; i < 32; i++)
            {
                ret |= ((ix & (1ull << i)) << (i + 1)) | ((iy & (1ull << i)) << i);
            }
            // cout << bitset<64>(ret) << endl;
            return ret;
        }

        unsigned long long overlap_bits() const
        {
            auto ix = static_cast<unsigned long long>(x);
            auto iy = static_cast<unsigned long long>(y);
            unsigned long long p[] = {ix, iy};
            return hilbert_c2i(2, 32, p);
        }

        long long mortonIndex() const
        {
            // Pun the x and y coordinates as integers: Just re-interpret the bits.
            //
            auto ix = static_cast<unsigned int>(x);
            auto iy = static_cast<unsigned int>(y);
            // cout << ix << " " << iy << endl;

            // Since we're assuming 2s complement arithmetic (99.99% of hardware today),
            // we'll need to convert these raw integer-punned floats into
            // their corresponding integer "indices".

            // Smear their sign bits into these for twiddling below.
            //
            const auto ixs = static_cast<int>(ix) >> 31;
            const auto iys = static_cast<int>(iy) >> 31;

            // This is a combination of a fast absolute value and a bias.
            //
            // We need to adjust the values so -FLT_MAX is close to 0.
            //
            ix = (((ix & 0x7FFFFFFFL) ^ ixs) - ixs) + 0x7FFFFFFFL;
            iy = (((iy & 0x7FFFFFFFL) ^ iys) - iys) + 0x7FFFFFFFL;

            // Now we have -FLT_MAX close to 0, and FLT_MAX close to UINT_MAX,
            // with everything else in-between.
            //
            // To make this easy, we'll work with x and y as 64-bit integers.
            //
            long long xx = ix;
            long long yy = iy;

            // Dilate and combine as usual...

            xx = (xx | (xx << 16)) & 0x0000ffff0000ffffLL;
            yy = (yy | (yy << 16)) & 0x0000ffff0000ffffLL;

            xx = (xx | (xx << 8)) & 0x00ff00ff00ff00ffLL;
            yy = (yy | (yy << 8)) & 0x00ff00ff00ff00ffLL;

            xx = (xx | (xx << 4)) & 0x0f0f0f0f0f0f0f0fLL;
            yy = (yy | (yy << 4)) & 0x0f0f0f0f0f0f0f0fLL;

            xx = (xx | (xx << 2)) & 0x3333333333333333LL;
            yy = (yy | (yy << 2)) & 0x3333333333333333LL;

            xx = (xx | (xx << 1)) & 0x5555555555555555LL;
            yy = (yy | (yy << 1)) & 0x5555555555555555LL;

            return xx | (yy << 1);
        }
    };

    typedef pair<Point, Point> Bounding_Box;

    struct diff_type{
        size_t add_cnt, remove_cnt;
        parlay::sequence<Point> add, remove;
        diff_type(){
            add_cnt = 0;
            remove_cnt = 0;
        }
        diff_type(size_t add_sz, size_t remove_sz){
            add_cnt = 0;
            remove_cnt = 0;
            add = parlay::sequence<Point>::uninitialized(add_sz);
            remove = parlay::sequence<Point>::uninitialized(remove_sz);
        }
        void add_point(Point &p, bool reverse = false){
            if (!reverse){
                parlay::assign_uninitialized(add[add_cnt++], p);
            }
            else{
                parlay::assign_uninitialized(remove[remove_cnt++], p);
            }
        }
        void remove_point(Point &p, bool reverse = false){
            if (!reverse){
                parlay::assign_uninitialized(remove[remove_cnt++], p);
            }
            else{
                parlay::assign_uninitialized(add[add_cnt++], p);
            }
        }
        void compact(){
            add.resize(add_cnt);
            remove.resize(remove_cnt);
        }
        void reset(){
            add_cnt = 0;
            remove_cnt = 0;
        }
        void reset(size_t add_sz, size_t remove_sz){
            add_cnt = 0;
            remove_cnt = 0;
            add = parlay::sequence<Point>::uninitialized(add_sz);
            remove = parlay::sequence<Point>::uninitialized(remove_sz);
        }
    };

}
