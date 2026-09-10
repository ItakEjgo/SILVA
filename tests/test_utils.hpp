#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>

#include <cpam/cpam.h>
#include <parlay/primitives.h>
#include <silva/geo/point.hpp>
#include <silva/geo/operations.hpp>
#include <silva/geo/io.hpp>
#include <parlay/internal/get_time.h>
#include <parlay/hash_table.h>
#include <silva/index/mvq.hpp>
#include <silva/index/cpamz.hpp>
#include <helper/time_loop.h>

#include <silva/core/hilbert.h>
#include <silva/index/cpambb.hpp>

#define TEST	//	print for correctness check

using namespace std;
using namespace geobase;

#include <silva/core/global_config.hpp>


