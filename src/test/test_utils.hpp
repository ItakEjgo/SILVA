#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>

#include <cpam/cpam.h>
#include <parlay/primitives.h>
#include "../geo/point.hpp"
#include "../geo/operations.hpp"
#include "../geo/io.hpp"
#include <parlay/internal/get_time.h>
#include <parlay/hash_table.h>
#include "../mvq.hpp"
#include "../cpamz.hpp"
#include "helper/time_loop.h"

#include "../hilbert.h"
#include "../cpambb.hpp"

#define TEST	//	print for correctness check

using namespace std;
using namespace geobase;

#include "../global_config.hpp"


