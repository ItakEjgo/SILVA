#pragma once

#include <bits/stdc++.h>

#include <cpam/cpam.h>
#include <parlay/primitives.h>
#include "geo/point.hpp"
#include "geo/operations.hpp"
#include "geo/io.hpp"
#include <parlay/internal/get_time.h>
#include "mvq.hpp"


namespace CPAMZ{

	using namespace std;

	using key_type = pair<unsigned long long, unsigned long long>;	// morton_id, id
	using geobase::Point;
	using val_type = Point;

	struct entry {
	using key_t = key_type;
	using val_t = val_type;
	using aug_t = key_type;

	static inline bool comp(key_t a, key_t b) { return a < b; }
	static aug_t get_empty() { return {0, 0}; }
	static aug_t from_entry(key_t k, val_t v) { return get_empty(); }
	static aug_t combine(aug_t a, aug_t b) { return std::max(a, b); }
	};

	//	cpam leaf size is [B, 2B]
	using zmap = cpam::pam_map<entry, 16>;
	using par = std::tuple<entry::key_t, entry::val_t>;

	template<class T, class MBR>
	auto filter_range(T &tree, MBR &query_mbr, bool use_hilbert = false){
		auto BL = query_mbr.first;
		auto UR = query_mbr.second;
		unsigned long long Z_min = 0, Z_max = 0;
		if (!use_hilbert){
			auto BR = geobase::Point(UR.x, BL.y), UL = geobase::Point(BL.x, UR.y);
			Z_min = std::min(BL.interleave_bits(), std::min(std::min(UR.interleave_bits(), UL.interleave_bits()), BR.interleave_bits()));
			Z_max = std::max(BL.interleave_bits(), std::max(std::max(UR.interleave_bits(), UL.interleave_bits()), BR.interleave_bits()));
		}
		else{
			unsigned long long p1[] = {static_cast<unsigned long long>(BL.x), static_cast<unsigned long long>(BL.y)};	//	BL
			unsigned long long p2[] = {static_cast<unsigned long long>(UR.x), static_cast<unsigned long long>(UR.y)};	//	UR
			unsigned long long pointlo[] = {p1[0], p1[1]};
			unsigned long long work[] = {p2[0], p2[1]};
  			hilbert_box_pt(2, sizeof(unsigned long long), 8 * sizeof(unsigned long long), 1, pointlo, work);
			Z_min = hilbert_c2i(2, 32, pointlo);
			unsigned long long pointhi[] = {p2[0], p2[1]};
			work[0] = p1[0], work[1] = p1[1];
  			hilbert_box_pt(2, sizeof(unsigned long long), 8 * sizeof(unsigned long long), 0, work, pointhi);
			Z_max = hilbert_c2i(2, 32, pointhi);
		}
		entry::key_t small(Z_min, 0);
		entry::key_t large(Z_max, numeric_limits<unsigned long long>::max());
		auto f2 = [&](par cur){ 
			return BL.x <= get<1>(cur).x && get<1>(cur).x <= UR.x &&
				BL.y <= get<1>(cur).y && get<1>(cur).y <= UR.y;
		};
		return zmap::filter(zmap::range(tree, small, large), f2);
	}

	template<typename M>
	auto plain_map_spatial_diff(M &lhs, M &rhs, geobase::Bounding_Box &query){
		auto l_pts = zmap::values(filter_range(lhs, query));
		auto r_pts = zmap::values(filter_range(rhs, query));
		// cout << "debug: " << l_pts.size() << ", " << r_pts.size() << endl;
		auto [add, remove] = merge_pts(l_pts, r_pts);
		return make_tuple(add, remove);
	}

	template<typename M>
	auto map_spatial_diff(M &lhs, M &rhs, geobase::Bounding_Box &query){
		auto filtered_lhs = filter_range(lhs, query);
		auto filtered_rhs = filter_range(rhs, query);
		auto add = zmap::values(zmap::map_difference(filtered_rhs, filtered_lhs));
		auto remove = zmap::values(zmap::map_difference(filtered_lhs, filtered_rhs));
		return make_tuple(add, remove);
	}

	template<typename M>
	auto map_diff(M &lhs, M &rhs, geobase::Bounding_Box &query){
		auto add = zmap::values(zmap::map_difference(rhs, lhs));
		auto remove = zmap::values(zmap::map_difference(lhs, rhs));
		return make_tuple(add, remove);
	}

	template<class T, class MBR>
	auto range_report(T &CPAMZ, MBR &query_mbr, bool use_hilbert = false){
		auto BL = query_mbr.first;
		auto UR = query_mbr.second;
		unsigned long long Z_min = 0, Z_max = 0;
		if (!use_hilbert){
			auto BR = geobase::Point(UR.x, BL.y), UL = geobase::Point(BL.x, UR.y);
			Z_min = std::min(BL.interleave_bits(), std::min(std::min(UR.interleave_bits(), UL.interleave_bits()), BR.interleave_bits()));
			Z_max = std::max(BL.interleave_bits(), std::max(std::max(UR.interleave_bits(), UL.interleave_bits()), BR.interleave_bits()));
		}
		else{
			unsigned long long p1[] = {static_cast<unsigned long long>(BL.x), static_cast<unsigned long long>(BL.y)};	//	BL
			unsigned long long p2[] = {static_cast<unsigned long long>(UR.x), static_cast<unsigned long long>(UR.y)};	//	UR
			unsigned long long pointlo[] = {p1[0], p1[1]};
			unsigned long long work[] = {p2[0], p2[1]};
  			hilbert_box_pt(2, sizeof(unsigned long long), 8 * sizeof(unsigned long long), 1, pointlo, work);
			Z_min = hilbert_c2i(2, 32, pointlo);
			unsigned long long pointhi[] = {p2[0], p2[1]};
			work[0] = p1[0], work[1] = p1[1];
  			hilbert_box_pt(2, sizeof(unsigned long long), 8 * sizeof(unsigned long long), 0, work, pointhi);
			Z_max = hilbert_c2i(2, 32, pointhi);
		}
		// cout << "BL: " << BL.x << " " << BL.y << " " << BL.overlap_bits() << endl;
		// cout << "UR: " << UR.x << " " << UR.y << " " << UR.overlap_bits() << endl;
		// cout << Z_min << " " << Z_max << endl;
		entry::key_t small(Z_min, 0);
		entry::key_t large(Z_max, numeric_limits<unsigned long long>::max());
		// T r2_map = zmap::range(CPAMZ, Z_min, Z_max);
		T r2_map = zmap::range(CPAMZ, small, large);
		// cout << r2_map.size() << endl;
		auto f2 = [&](par cur){ 
			return BL.x <= get<1>(cur).x && get<1>(cur).x <= UR.x &&
				BL.y <= get<1>(cur).y && get<1>(cur).y <= UR.y;
		};
		r2_map = zmap::filter(r2_map, f2);
		auto ret = zmap::values(r2_map); 
		return ret;
	}

	template<class T, class MBR>
	auto zMAP_range_report(T &zMAP, MBR query_mbr){
		auto BL = query_mbr.first;
		auto UR = query_mbr.second;
		auto BR = geobase::Point(UR.x, BL.y), UL = geobase::Point(BL.x, UR.y);
		auto Z_min = std::min(BL.mortonIndex(), std::min(std::min(UR.mortonIndex(), UL.mortonIndex()), BR.mortonIndex()));
		auto Z_max = std::max(BL.mortonIndex(), std::max(std::max(UR.mortonIndex(), UL.mortonIndex()), BR.mortonIndex()));

		parlay::sequence<Point> ret = {};
		for (auto low_it = zMAP.lower_bound(Z_min); low_it != zMAP.upper_bound(Z_max); low_it++){
			if (BL.x <= get<1>(*low_it).x && get<1>(*low_it).x <= UR.x &&
				BL.y <= get<1>(*low_it).y && get<1>(*low_it).y <= UR.y ){
				ret.emplace_back(get<1>(*low_it));
			}
		}
		return ret;
	}
	
	template<class PT>
	auto zMAP_init(PT &P){
		size_t n = P.size();
		map<long long, geobase::Point> m1;

		for (size_t i = 0; i < n; i++){
			m1[P[i].morton_id] = P[i];
		}
		return m1;
	}

	template<typename M>
	auto map_diff(M &lhs, M &rhs){
		auto add = zmap::values(zmap::map_difference(rhs, lhs));
		auto remove = zmap::values(zmap::map_difference(lhs, rhs));
		return make_tuple(add, remove);
	}
	
	template<class PT>
	auto CPAMZ_init(PT &P, bool use_hilbert = false){
		size_t n = P.size();
		
		parlay::parallel_for(0, n, [&](int i){
			P[i].morton_id = use_hilbert ? P[i].overlap_bits() : P[i].interleave_bits();
		});

		parlay::sequence<par> entries(n);
		parlay::parallel_for(0, n, [&](int i){
			entries[i] = {{P[i].morton_id, P[i].id}, P[i]};
			// entries[i] = {P[i]->id, P[i]};
		});
		zmap m1(entries);
		return m1;
	}	

	template<typename PT, typename M>
	auto CPAMZ_insert(PT &P, M &mmp, bool use_hilbert = false){
		size_t n = P.size();

		parlay::parallel_for(0, n, [&](int i){
			P[i].morton_id = use_hilbert ? P[i].overlap_bits() : P[i].interleave_bits();
		});
		parlay::sequence<par> insert_pts(n);
		parlay::parallel_for(0, n, [&](int i){
			insert_pts[i] = {{P[i].morton_id, P[i].id}, P[i]};
			// insert_pts[i] = {P[i]->id, P[i]};
		});
		auto m2 = zmap::multi_insert(mmp, insert_pts);

		return m2;
	}

	template<typename PT, typename M>
	auto CPAMZ_delete(PT &P, M &mmp, bool use_hilbert = false){
		size_t n = P.size();

		parlay::parallel_for(0, n, [&](int i){
			P[i].morton_id = use_hilbert ? P[i].overlap_bits() : P[i].interleave_bits();
		});
		parlay::sequence<par> delete_pts(n);
		// parlay::sequence<pair<unsigned long long, long long> > delete_pts(n);
		parlay::parallel_for(0, n, [&](int i){
			delete_pts[i] = {{P[i].morton_id, P[i].id}, P[i]};
			// delete_pts[i] = {P[i]->morton_id, P[i]->id};
			// insert_pts[i] = {P[i]->id, P[i]};
		});
		auto m2 = zmap::multi_delete(mmp, delete_pts);

		return m2;
	}

	template<typename PT, typename M>
	auto CPAMZ_commit(M &mmp, PT &P_insert, PT &P_delete){
		// parlay::internal::timer t("CPAMZ breakdown", true);
		auto new_ver = CPAMZ_delete(P_delete, mmp);	//	new	version
		// t.next("delete time");
		new_ver = CPAMZ_insert(P_insert, new_ver); 
		// t.next("insert time");
		return new_ver;
	}

	template<typename M>
	auto CPAMZ_merge(M &base, M &v1, M &v2){
		parlay::internal::timer t("merge", false);
		auto [add1, remove1] = map_diff(base, v1);
		auto [add2, remove2] = map_diff(base, v2); 
		t.next("diff time");

		auto [insert1, delete1, update1] = geobase::filter_diff_results(add1, remove1);
		auto [insert2, delete2, update2] = geobase::filter_diff_results(add2, remove2);

		auto [no_conflict_insert, conflict_insert] = geobase::merge_by_id_with_conflict(insert1, insert2);	// insert conflict
		auto [no_conflict_update, conflict_update] = geobase::merge_by_id_with_conflict(update1, update2);	// update conflict
		//	delate = delete + update
		auto [no_conflict_delete, conflict_delete] = geobase::merge_by_id_with_conflict(delete1, update2);	//	delete and update conflict
		auto [no_conflict_delete1, conflict_delete1] = geobase::merge_by_id_with_conflict(delete2, update1);	//	delete and update conflict, another case

		no_conflict_delete.append(no_conflict_delete1);
		conflict_delete.append(conflict_delete1);
		t.next("conflict time");

		auto new_ver = CPAMZ_commit(base, no_conflict_insert, no_conflict_delete);
		t.next("commit time");

		return make_tuple(new_ver, conflict_insert, conflict_update, conflict_delete);
	}

	//	return size of interior nodes and sizeof leaf nodes size, respectively
	auto size_in_bytes(){
		size_t inte_used = zmap::GC::used_node();
		size_t internal_nodes_space = sizeof(typename zmap::GC::regular_node) * inte_used;
		auto [used, unused] = parlay::internal::get_default_allocator().stats();
		return make_tuple(internal_nodes_space, used);
	}


}