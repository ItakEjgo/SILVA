#pragma once
#include <cstddef>
#include "geo/point.hpp"
#include "geo/operations.hpp"
#include "geo/io.hpp"

namespace mvq {

struct Config {
    size_t leaf_size = 32;
    size_t maxSize = 100;
    geobase::Bounding_Box largest_mbr;

    static Config& get() {
        static Config instance;
        return instance;
    }
};

struct ProfilerStats {
    geobase::break_down zd_build_break_down;
    geobase::break_down cpam_build_break_down;
    double zd_leaf_copy_time = 0.0;
    double zd_inte_copy_time = 0.0;
    double leaf_time = 0.0;
    double inte_time = 0.0;
    size_t visited_leaf = 0;
    size_t visited_inte = 0;

    static ProfilerStats& get() {
        static ProfilerStats instance;
        return instance;
    }
};

} // namespace mvq