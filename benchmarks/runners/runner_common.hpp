#pragma once
#include <vector>
#include <atomic>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <parlay/sequence.h>
#include <silva/geo/point.hpp>
#include <silva/geo/operations.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> BoostPoint;
typedef std::pair<BoostPoint, size_t> Value;

extern std::atomic<size_t> boost_live_mem;
