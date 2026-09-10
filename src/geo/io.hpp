#include <parlay/parallel.h>
#include <parlay/sequence.h>
#pragma once
#include "point.hpp"
#include "operations.hpp"
namespace geobase {
    template <class T>
    auto read_pts(T &P, ifstream &fin, bool real_data = false)
    {
        if (!real_data)
        {
            size_t n, d;
            fin >> n >> d;
            P.resize(n);
            size_t id;
            FT x, y;
            FT x_min(FT_INF_MAX), x_max(FT_INF_MIN), y_min(FT_INF_MAX), y_max(FT_INF_MIN);
            for (size_t i = 0; i < n; i++)
            {
                // fin >> id >> x >> y;
                fin >> x >> y;
                x_max = max(x_max, x);
                x_min = min(x_min, x);
                y_max = max(y_max, y);
                y_min = min(y_min, y);
                id = i;
                auto cur_p = Point(id, x, y);
                P[i] = cur_p;
            }
            return Bounding_Box({Point(x_min, y_min), Point(x_max, y_max)});
        }
        else
        {
            size_t id;
            FT x, y;
            FT x_min(FT_INF_MAX), x_max(FT_INF_MIN), y_min(FT_INF_MAX), y_max(FT_INF_MIN);
            P.clear();
            while (fin >> id >> x >> y)
            {
                x *= 1000000;
                y *= 1000000;
                if (x < 0 || y < 0)
                {
                    continue; // ignore outliers
                }
                x_max = max(x_max, x);
                x_min = min(x_min, x);
                y_max = max(y_max, y);
                y_min = min(y_min, y);
                auto cur_p = Point(id, x, y);
                P.emplace_back(cur_p);
            }
            return Bounding_Box({Point(x_min, y_min), Point(x_max, y_max)});
        }
    }

    template <class T>
    void print_binary(T x)
    {
        cout << bitset<sizeof(x) * 8>(x) << endl;
    }

    template <class MBR>
    void print_mbr(MBR &mbr)
    {
        cout << fixed << setprecision(6) << "[(" << mbr.first.x << ", " << mbr.first.y << "), (" << mbr.second.x << ", " << mbr.second.y << ")" << "]" << endl;
    }

    template <class T>
    auto generate_range_query(T &P, size_t n)
    {
        auto id1 = rand() % n;
        auto id2 = rand() % n;
        FT xmin = min(P[id1].x, P[id2].x),
           xmax = max(P[id1].x, P[id2].x),
           ymin = min(P[id1].y, P[id2].y),
           ymax = max(P[id1].y, P[id2].y);
        return make_pair(Point(xmin, ymin), Point(xmax, ymax));
    }

    template <class In>
    auto read_range_query(In qry_in)
    {
        ifstream fin(qry_in);
        size_t n, d;
        fin >> n >> d;
        parlay::sequence<Bounding_Box> ret(n);
        for (size_t i = 0; i < n; i++)
        {
            fin >> ret[i].first.x >> ret[i].first.y >> ret[i].second.x >> ret[i].second.y;
        }
        return ret;
    }

    template <class In>
    auto read_range_query(In qry_in, size_t q_type, size_t &maxSize, bool is_real = false)
    {
        ifstream fin(qry_in);
        if (q_type == 8)
        { // range report, need maxSize
            fin >> maxSize;
        }
        size_t n, d;
        fin >> n >> d;
        parlay::sequence<Bounding_Box> ret(n);
        parlay::sequence<size_t> cnt(n);
        for (size_t i = 0; i < n; i++)
        {
            fin >> cnt[i] >> ret[i].first.x >> ret[i].first.y >> ret[i].second.x >> ret[i].second.y;
        }
        return make_tuple(cnt, ret);
    }

    // data type/helper functions for nearest neighbor search
    }
