#pragma once
#include "test_utils.hpp"
namespace ZDTest{

	template<typename PT>
	void spatial_join_test(PT &P1, PT &P2, FT &point_dis){
		auto P_set = get_sorted_points(P1);
		mvq::Tree zdtree(mvq::Config::get().leaf_size);
		zdtree.build(P_set);	// initial version for P1

		P_set = get_sorted_points(P2);
		mvq::Tree zdtree2(mvq::Config::get().leaf_size);
		zdtree2.build(P_set);

		parlay::sequence<pair<geobase::Point, geobase::Point> > join_res = {};

		auto avg_time = time_loop(
			3, 1.0,
			[&](){
				join_res.clear();
			},
			[&](){
				zdtree.two_version_spatial_join(zdtree.root, zdtree2.root, point_dis, join_res, mvq::Config::get().largest_mbr, mvq::Config::get().largest_mbr);
			},
			[&]{}
		);
		cout << fixed << setprecision(6) << "[zdtree] spatial-join time (avg): " << avg_time << endl;

		#ifdef TEST
			auto cmp = [](const auto &lhs, const auto &rhs){
				if (lhs.first.id != rhs.first.id) return lhs.first.id < rhs.first.id;
				return lhs.second.id < rhs.second.id;
			};
			join_res = parlay::sort(join_res, cmp);
			parlay::sequence<pair<geobase::Point, geobase::Point> > bf_res = {};
			
			for (auto &pt1: P1){
				for (auto &pt2: P2){
					if (dcmp(geobase::point_point_sqrdis(pt1, pt2) - point_dis * point_dis) <= 0){
						bf_res.emplace_back(pt1, pt2);
					}
				}
			}

			bf_res = parlay::sort(bf_res, cmp);
			if (join_res != bf_res){
				cout << "[ERROR] Incorrect Join Results" << endl;
				cout << "Join Res: " << join_res.size() << endl;
				for (auto &par: join_res){
					cout << fixed << setprecision(6) << sqrt(point_point_sqrdis(par.first, par.second)) << " | " << par.first << ", " << par.second << endl;
				}
				cout << "BF Res: " << bf_res.size() << endl;
				for (auto &par: bf_res){
					cout << fixed << setprecision(6) << sqrt(point_point_sqrdis(par.first, par.second)) << " | " << par.first << ", " << par.second << endl;
				}
			}
			else{
				cout << "[INFO] Correct Join Restuls" << endl;
				cout << "join/bf size: " << join_res.size() << ", " << bf_res.size() << endl;
			}
		#endif
	}

	template<typename PT, typename RQ>
    void plain_spatial_diff_test_latency(PT &P, RQ &range_queries, parlay::sequence<size_t> &batch_sizes, size_t &insert_ratio, bool dual_traverse = false){
        /*  build tree */
		auto P_set = get_sorted_points(P);
		mvq::Tree zdtree(mvq::Config::get().leaf_size);
		zdtree.build(P_set);	// initial version

		auto max_batch_size = batch_sizes[batch_sizes.size() - 1];
        /* get insert, delete points */
        auto P_test = geobase::shuffle_point(P, max_batch_size);
        auto [P_insert_set, P_delete_set] = geobase::split_insert_delete(P_test, insert_ratio, P.size());

        /* get insert, delete points */
		for (auto &batch_size: batch_sizes){
			cout << "[INFO] Batch Size: " << batch_size << endl;
			auto insert_num = batch_size / 10 * insert_ratio;
			auto delete_num = batch_size / 10 * (10 - insert_ratio);

			auto P_insert = P_insert_set.substr(0, insert_num);
			auto P_delete = P_delete_set.substr(0, delete_num);

			auto P_newver = geobase::collect_newver_point(P, P_insert, P_delete);
	
			P_set = get_sorted_points(P_newver);
			mvq::Tree newtree(mvq::Config::get().leaf_size);
			newtree.build(P_set);
	
			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());
			// parlay::sequence<size_t> pts1Cnt(range_queries.size());
			// parlay::sequence<size_t> pts2Cnt(range_queries.size());
			// map<size_t, size_t> pts_map, diff_map;

			auto pts1 = parlay::sequence<Point>::uninitialized(2 * mvq::Config::get().maxSize); 
			auto pts2 = parlay::sequence<Point>::uninitialized(2 * mvq::Config::get().maxSize); 
			for (size_t i = 0; i < range_queries.size(); i++){
				auto avg_time = time_loop(
					3, 1.0,
					[&](){},
					[&](){
						if (!dual_traverse){
							size_t cnt1 = 0, cnt2 = 0;
							zdtree.range_report(range_queries[i], mvq::Config::get().largest_mbr, cnt1, pts1);
							pts1.resize(cnt1);
							newtree.range_report(range_queries[i], mvq::Config::get().largest_mbr, cnt2, pts2);
							pts2.resize(cnt2);
							diff_type ret_diff(mvq::Config::get().maxSize, mvq::Config::get().maxSize);
							merge_pts(pts1, pts2, ret_diff);
							ret_diff.compact();
							addCnt[i] = ret_diff.add.size();
							removeCnt[i] = ret_diff.remove.size();
						}
						else{
							diff_type ret_diff(mvq::Config::get().maxSize, mvq::Config::get().maxSize);
							zdtree.spatial_two_version_diff(zdtree.root, newtree.root, range_queries[i], mvq::Config::get().largest_mbr, ret_diff);
							ret_diff.compact();
							addCnt[i] = ret_diff.add.size();
							removeCnt[i] = ret_diff.remove.size();
						}
					},
					[&]{}
				);
				cout << fixed << setprecision(6) << i << " " << avg_time << endl;
			}
			
	
			
			#ifdef TEST
				string file_name = "output/zd_spatial_diff_plain-" + to_string(batch_size); 
				if (dual_traverse){
					file_name = file_name + "-dual";
				}
				ofstream spatialDiffOut(file_name);
				for (size_t i = 0; i < range_queries.size(); i++){
					spatialDiffOut<< addCnt[i] << " " << removeCnt[i] << endl;
				}
			#endif
		}
    }


	template<typename PT, typename RQ>
    void spatial_diff_test_latency(PT &P, RQ &range_queries, parlay::sequence<size_t> &batch_sizes, size_t &insert_ratio){
        /*  build tree */
		auto P_set = get_sorted_points(P);
		mvq::Tree zdtree(mvq::Config::get().leaf_size);
		zdtree.build(P_set);	// initial version

		auto max_batch_size = batch_sizes[batch_sizes.size() - 1];
        /* get insert, delete points */
        auto P_test = geobase::shuffle_point(P, max_batch_size);
        auto [P_insert_set, P_delete_set] = geobase::split_insert_delete(P_test, insert_ratio, P.size());
		
		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) batch_size = P.size();

			cout << "[INFO] Batch Size: " << batch_size << endl;
			auto insert_num = batch_size / 10 * insert_ratio;
			auto delete_num = batch_size / 10 * (10 - insert_ratio);

			auto P_insert = P_insert_set.substr(0, insert_num);
			auto P_delete = P_delete_set.substr(0, delete_num);
			
			// get version 1 by deletion
			auto P_delete_sorted = get_sorted_points(P_delete);
			auto new_ver = zdtree.multi_version_batch_delete_sorted(P_delete_sorted, zdtree.root);	//	new version
			
			// get version 2 by insertion
			auto P_insert_sorted = get_sorted_points(P_insert);
			auto new_ver2 = zdtree.multi_version_batch_insert_sorted(P_insert_sorted, new_ver);
			
			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());

			diff_type ret_diff(mvq::Config::get().maxSize, mvq::Config::get().maxSize);

			for (size_t i = 0; i < range_queries.size(); i++){
				auto avg_time = time_loop(
					3, 1.0,
					[&](){},
					[&](){
						ret_diff.reset(mvq::Config::get().maxSize, mvq::Config::get().maxSize);
						zdtree.spatial_two_version_diff(zdtree.root, new_ver2, range_queries[i], mvq::Config::get().largest_mbr, ret_diff);
						ret_diff.add.resize(ret_diff.add_cnt);
						ret_diff.remove.resize(ret_diff.remove_cnt);
						addCnt[i] = ret_diff.add.size();
						removeCnt[i] = ret_diff.remove.size();
					},
					[&]{}
				);
				cout << fixed << setprecision(6) << i << " " << avg_time << endl;
			}

			
			#ifdef TEST
				string file_name = "output/zd_spatial_diff-" + to_string(batch_size); 
				ofstream spatialDiffOut(file_name);
				for (size_t i = 0; i < range_queries.size(); i++){
					spatialDiffOut<< addCnt[i] << " " << removeCnt[i] << endl;
				}
			#endif

		}
    }

	/*	50% insertion, 50% deletion	*/
	template<typename PT, typename RQ>
	auto spatial_diff_test_fix_size(PT P,  RQ &range_queries, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto P_set = get_sorted_points(P);
		mvq::Tree zdtree(mvq::Config::get().leaf_size);
		zdtree.build(P_set);	// initial version

		parlay::sequence<size_t> addCnt(range_queries.size());
		parlay::sequence<size_t> removeCnt(range_queries.size());

		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) break;
			cout << "[ratio (insert/delete)]: " << batch_size << "/" << 10 - batch_size << endl;

			auto insert_num = P.size() / 100 * batch_size;	//	0.5 insert
			auto delete_num = P.size() / 100 * (10 - batch_size);	//	0.5 delete

			auto P_insert = P.substr(0, insert_num);
			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += P.size();
			});

			auto P_delete = P.substr(P.size() - delete_num, delete_num);

			cout << "# of insertion/deletion: " << P_insert.size() << ", " << P_delete.size();

			// get version 1 by deletion
			auto P_delete_sorted = get_sorted_points(P_delete);
			auto new_ver = zdtree.multi_version_batch_delete_sorted(P_delete_sorted, zdtree.root);	//	new version

			// get version 2 by insertion
			auto P_insert_sorted = get_sorted_points(P_insert);
			auto new_ver2 = zdtree.multi_version_batch_insert_sorted(P_insert_sorted, new_ver);

	    	auto diff_avg = time_loop(
		    	3, 1.0, [&]() {},
		    	[&]() {
					parlay::parallel_for(0, range_queries.size(), [&](int i){
						diff_type ret_diff(mvq::Config::get().maxSize, mvq::Config::get().maxSize);
						zdtree.spatial_two_version_diff(zdtree.root, new_ver2, range_queries[i], mvq::Config::get().largest_mbr, ret_diff);
						addCnt[i] = ret_diff.add.size();
						removeCnt[i] = ret_diff.remove.size();
					});
		    	},
	    		[&](){} 
			);
			cout << fixed << setprecision(6) << "[zdtree]: spatial diff time (avg): " << diff_avg << endl;
		}
	}

	/*	50% insertion, 50% deletion	*/
	template<typename PT, typename RQ>
	auto spatial_diff_test_fix_ratio(PT P,  RQ &range_queries, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto P_set = get_sorted_points(P);
		mvq::Tree zdtree(mvq::Config::get().leaf_size);
		zdtree.build(P_set);	// initial version

		parlay::sequence<size_t> addCnt(range_queries.size());
		parlay::sequence<size_t> removeCnt(range_queries.size());

		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) break;
			cout << "[batch-size]: " << batch_size << endl;

			auto insert_num = batch_size / 2;	//	0.5 insert
			auto delete_num = batch_size - insert_num;	//	0.5 delete

			auto P_insert = P.substr(0, insert_num);
			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += P.size();
			});

			auto P_delete = P.substr(P.size() - delete_num, delete_num);

			cout << "# of insertion/deletion: " << P_insert.size() << ", " << P_delete.size();

			// get version 1 by deletion
			auto P_delete_sorted = get_sorted_points(P_delete);
			auto new_ver = zdtree.multi_version_batch_delete_sorted(P_delete_sorted, zdtree.root);	//	new version

			// get version 2 by insertion
			auto P_insert_sorted = get_sorted_points(P_insert);
			auto new_ver2 = zdtree.multi_version_batch_insert_sorted(P_insert_sorted, new_ver);

	    	auto diff_avg = time_loop(
		    	3, 1.0, [&]() {},
		    	[&]() {
					parlay::parallel_for(0, range_queries.size(), [&](int i){
						diff_type ret_diff(mvq::Config::get().maxSize, mvq::Config::get().maxSize);
						zdtree.spatial_two_version_diff(zdtree.root, new_ver2, range_queries[i], mvq::Config::get().largest_mbr, ret_diff);
						addCnt[i] = ret_diff.add.size();
						removeCnt[i] = ret_diff.remove.size();
					});
		    	},
	    		[&](){} 
			);
			cout << fixed << setprecision(6) << "[zdtree]: spatial diff time (avg): " << diff_avg << endl;
		}
	}



	template<typename PT, typename RQ>
	auto spatial_diff_test(PT P,  RQ &range_queries, parlay::sequence<size_t> &batch_sizes, bool &early_end, bool use_hilbert = false){
		auto P_set = get_sorted_points(P);
		mvq::Tree zdtree(mvq::Config::get().leaf_size);
		zdtree.build(P_set);	// initial version

		parlay::sequence<size_t> addCnt(range_queries.size());
		parlay::sequence<size_t> removeCnt(range_queries.size());

		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) batch_size = P.size();
			cout << "[batch-size]: " << batch_size << endl;

			auto P_insert = P.substr(0, batch_size);
			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += P.size();
			});

			auto P_delete = P.substr(P.size() - batch_size, batch_size);

			// get version 1 by deletion
			auto P_delete_sorted = get_sorted_points(P_delete);
			auto new_ver = zdtree.multi_version_batch_delete_sorted(P_delete_sorted, zdtree.root);	//	new version

			// get version 2 by insertion
			auto P_insert_sorted = get_sorted_points(P_insert);
			auto new_ver2 = zdtree.multi_version_batch_insert_sorted(P_insert_sorted, zdtree.root);

			// cout << P_insert.size() << ", " << P_delete.size() << endl;
			// cout << zdtree.collect_records(new_ver).size() << ", " << zdtree.collect_records(new_ver2).size() << endl;



			// cout << "[INFO] add set size, P_insert size: " << add.size() << ", " << P_insert.size() << endl;
			// cout << "[INFO] remove set size, P_delete size: " << remove.size() << ", " << P_delete.size() << endl;
			if (!early_end){
				decltype(zdtree.root) commit_ver;
	    		auto commit_avg = time_loop(
		    		3, 1.0, [&]() {
						commit_ver.reset();
					},
		    		[&]() {
						commit_ver = zdtree.commit(zdtree.root, P_insert, P_delete);
		    		},
	    			[&](){} 
				);

				// cout << "init size = " << zdtree.collect_records(zdtree.root).size() << endl;
				// cout << "commit size = " << zdtree.collect_records(commit_ver).size() << endl;
				// cout << "commit finished." << endl;

				decltype(zdtree.root) merge_ver;
				parlay::sequence<Point> conflict_insert, conflict_update, conflict_delete;

				auto merge_avg = time_loop(
		    		3, 1.0, [&]() {
						merge_ver.reset();
						conflict_insert.clear();
						conflict_update.clear();
						conflict_delete.clear();
					},
		    		[&]() {
						tie(merge_ver, conflict_insert, conflict_update, conflict_delete) = zdtree.merge(zdtree.root, new_ver, new_ver2);
		    		},
	    			[&](){} 
				);

				cout << "[INFO] commit, merge size: " << zdtree.collect_records(commit_ver).size() << ", " << zdtree.collect_records(merge_ver).size() << endl;
				cout << fixed << setprecision(6) << "[zdtree]: spatial commit time (avg): " << commit_avg << endl;
				cout << fixed << setprecision(6) << "[zdtree]: spatial merge time (avg): " << merge_avg << endl;
			}

		}

		// return ret;
	}

	template<typename PT>
	void diff_test(PT P, int batch_percent = 10, bool use_hilbert = false){
		auto P_set = get_sorted_points(P);
		mvq::Tree zdtree(mvq::Config::get().leaf_size);
		zdtree.build(P_set);	// initial version

		auto batch_size = P.size() * batch_percent / 100;	//	insertion 10%

		auto P_insert = P.substr(0, batch_size);
		parlay::parallel_for(0, P_insert.size(), [&](size_t j){
			P_insert[j].id += P.size();
		});
		
		auto P_insert_sorted = get_sorted_points(P_insert);
		auto new_ver = zdtree.multi_version_batch_insert_sorted(P_insert_sorted, zdtree.root);

		auto delete_size = std::min(P.size(), (size_t)2 * batch_size);
		auto P_delete = P.substr(0, delete_size);
		// auto P_delete = P.substr(0, 10);
		auto P_delete_sorted = get_sorted_points(P_delete);
		new_ver = zdtree.multi_version_batch_delete_sorted(P_delete_sorted, new_ver);	//	new version

		auto add_sz = 0, remove_sz = 0;
	    auto zd_diff_avg = time_loop(
		    3, 1.0, [&]() {},
		    [&]() {
				auto [add, remove] = zdtree.diff(zdtree.root, new_ver, 64);
				add_sz = add.size();
				remove_sz = remove.size();
		    },
	    [&](){} );

		cout << "add size: " << add_sz << endl;
		cout << "remove size: " << remove_sz << endl;
		cout << fixed << setprecision(6) << "[zdtree] diff time (avg): " << zd_diff_avg << endl;
	}
	
	

	template<typename PT>
	void multi_version_query_test(PT P, string query_dir, int batch_percent = 10, int version_num = 6){
		// build zdtree initial version
		auto P_set = get_sorted_points(P);
		mvq::Tree zdtree(mvq::Config::get().leaf_size);
		zdtree.build(P_set);
		cout << "[INFO] Tree build finished." << endl;
		zdtree.multi_version_roots.emplace_back(zdtree.root);

		auto num_insert_version = version_num / 2;
		auto num_delete_version = version_num - num_insert_version;
		auto batch_size = P.size() * batch_percent / 100;
		
		shared_ptr<mvq::BaseNode> new_ver = zdtree.root;
		// insert 3 versions
		for (auto i = 0; i != num_insert_version; i++){
			auto P_insert = P.substr(i * batch_size, batch_size);

			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += (i + 1) * P.size();
			});

			auto P_insert_sorted = get_sorted_points(P_insert);
			new_ver = zdtree.multi_version_batch_insert_sorted(P_insert_sorted, new_ver);
			zdtree.multi_version_roots.emplace_back(new_ver);
		}
		cout << "[INFO] Tree insertion finished." << endl;

		//	delete 3 versions
		for (auto i = 0; i != num_delete_version; i++){
			auto P_delete = P.substr(i * batch_size, batch_size);

			parlay::parallel_for(0, P_delete.size(), [&](size_t j){
				P_delete[j].id += (i + 1) * P.size();
			});

			auto P_delete_sorted = get_sorted_points(P_delete);
			new_ver = zdtree.multi_version_batch_delete_sorted(P_delete_sorted, new_ver);
			zdtree.multi_version_roots.emplace_back(new_ver);
		}

		cout << "[INFO] Tree deletion finished." << endl;


		// zdtree.print_leaf(zdtree.multi_version_roots[0]);
		// zdtree.print_leaf(zdtree.multi_version_roots[2]);
		// return;

		// range count query test
		for (auto i = 0; i != 3; i++){	//	0, 1, 2 represent small, median, and large regions
			auto cur_query_dir = query_dir + "/1.in-" + to_string(i) + ".qry";
			auto [cnt, range_count_querys] = read_range_query(cur_query_dir, 8, mvq::Config::get().maxSize);
			// range_count_querys = range_count_querys.substr(0, 1000);
			// cout << "test i" << endl;

			// cout << zdtree.multi_version_roots.size() << endl;
			for (size_t j = 0; j < zdtree.multi_version_roots.size(); j++){
				parlay::sequence<size_t> rangeCnt(range_count_querys.size());
				auto rangeCnt_avg = time_loop(
					3, 1, [&]() {},
					// 1, 0, [&]() {},
					[&]() {
						parlay::parallel_for(
							0, range_count_querys.size(),
							[&]( size_t k ) {
								rangeCnt[k] = zdtree.range_count(zdtree.multi_version_roots[j], range_count_querys[k], mvq::Config::get().largest_mbr);
						});
					},
				[&](){} );

				cout << fixed << setprecision(6) << "zdtree range count time (avg) for region " << i << " on version " << j << ": " << rangeCnt_avg << endl;
				auto output_name = "mvq_range_count-" + to_string(i) + "-on-" + to_string(j) + ".txt";
				ofstream regionCntOut(output_name);
				for (size_t k = 0; k < range_count_querys.size(); k++){
					regionCntOut << rangeCnt[k] << endl;
				}
			}
		}

		// range report test
		for (auto i = 0; i != 2; i++){
			auto cur_query_dir = query_dir + "/1.in-" + to_string(i) + ".qry";
			auto [cnt, range_report_querys] = read_range_query(cur_query_dir, 8, mvq::Config::get().maxSize);

			for (size_t j = 0; j < zdtree.multi_version_roots.size(); j++){
				parlay::sequence<parlay::sequence<Point> > rangeReport(range_report_querys.size());
				parlay::sequence<size_t> rangeReportCnt(range_report_querys.size(), 0);
				for (size_t k = 0; k < range_report_querys.size(); k++){
					rangeReport[k].resize(mvq::Config::get().maxSize);	
				}
				auto rangeReport_avg = time_loop(
					3, 1.0, [&]() {},
					[&]() {
						parlay::parallel_for(
							0, range_report_querys.size(),
							[&]( size_t k ) {
								rangeReportCnt[k] = 0;
								zdtree.range_report(zdtree.multi_version_roots[j], range_report_querys[k], mvq::Config::get().largest_mbr, rangeReportCnt[k], rangeReport[k]);
						});
					},
				[&](){} );
				cout << fixed << setprecision(6) << "zdtree range report time (avg) for region " << i << " on version " << j << ": " << rangeReport_avg << endl;
				auto output_name = "mvq_range_report-" + to_string(i) + "-on-" + to_string(j) + ".txt";
				ofstream regionCntOut(output_name);
				for (size_t k = 0; k < range_report_querys.size(); k++){
					regionCntOut << rangeReportCnt[k] << endl;
				}
			}
		}
		
		
	}

	template<typename PT>
	void multi_version_test(PT P, string dir, int start_year = 14, int version_num = 5){
		auto cur_year = start_year;
		
		parlay::sequence<geobase::Point> P_delete[version_num], P_insert[version_num], P_update[version_num], P_updove[version_num];

		for (auto i = 0; i != version_num; i++){
			// cout << "[INFO] Year: " << cur_year << "-" << cur_year + 1 << " status:" << endl;
			auto delete_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-delete.txt";
			auto insert_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-insert.txt";
			auto update_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-update.txt";
			auto updove_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-update_remove.txt";
			// cout << delete_file_name << endl << insert_file_name << endl;
			ifstream fin_delete(delete_file_name);
			ifstream fin_insert(insert_file_name);
			ifstream fin_update(update_file_name);
			ifstream fin_updove(updove_file_name);

			auto delete_mbr = read_pts(P_delete[i], fin_delete, 1);
			auto insert_mbr = read_pts(P_insert[i], fin_insert, 1);
			auto update_mbr = read_pts(P_update[i], fin_update, 1);
			auto updove_mbr = read_pts(P_updove[i], fin_updove, 1);

			// cout << "# of delete pts: " << P_delete[i].size() << endl;
			// cout << "# of insert pts: " << P_insert[i].size() << endl;

			if (P_update[i].size() != P_updove[i].size()){
				cout << "[ERROR]: inconsistent # of update pts!" << endl;
			}
			// cout << "# of update pts: " << P_update[i].size() << endl;

			P_delete[i].append(P_updove[i]);
			P_insert[i].append(P_update[i]);

			cur_year += 1;
			delete_mbr = insert_mbr; // useless, just remove warning
			update_mbr = updove_mbr;
		}

		mvq::Tree zdtree(mvq::Config::get().leaf_size);

		/* Build initial version */
		auto build_avg = time_loop(
            3, 1.0, [&]() {
				zdtree.clear();
			},
            [&]() {
                auto P_set = get_sorted_points(P);
                zdtree.build(P_set);
            },
   	    	[&](){
			});

		zdtree.multi_version_roots.emplace_back(zdtree.root);	//	store the initial version 0
		auto tot = zdtree.num_of_nodes();
		auto stat = zdtree.get_tree_statistics();

		cout << "[zdtree init build time]: " << fixed << setprecision(6) << build_avg << " Seconds" << endl;
		cout << "[zdtree num of total tree nodes]: " << tot << endl;
		cout << "[zdtree node nums]: " << stat.num_inte_nodes << " interior nodes, " << stat.num_leaf_nodes << " leaf nodes" << endl;
		cout << "[zdtree memory usage]: " << stat.get_total() << " MB" << endl;
		// cout << "[memory usage for inte nodes]: " << 1.0 * stat.mem_inte_nodes / 1024.0 / 1024.0 << " MB" << endl 
		//      << "[memory usage for leaf nodes]: " << 1.0 * stat.mem_leaf_nodes / 1024.0 / 1024.0 << " MB"  << endl;
		
		shared_ptr<mvq::BaseNode> new_ver = zdtree.root;

		for (auto i = 0; i < version_num; i++){
			cout << "dealing with version " << i + 1 << ":" << endl;

			auto commit_avg = time_loop(
				3, 1.0, 
				[&]() {
					new_ver.reset();
				},
				[&]() {
					new_ver = zdtree.commit(zdtree.multi_version_roots[i], P_insert[i], P_delete[i]);
				},
				[&](){});

			zdtree.multi_version_roots.emplace_back(new_ver);
			cout << "[new ver commit time]: " << fixed << setprecision(6) << commit_avg << " Seconds" << endl;

			tot = zdtree.num_of_nodes();
			auto cur_stat = zdtree.get_tree_statistics();

			cout << "[num of total tree nodes]: " << tot << endl;
			cout << "[zdtree memory usage]: " << cur_stat.get_total() - stat.get_total() << " MB" << endl;
			cout << "[zdtree tot inte memory usage]: " << 1.0 * cur_stat.mem_inte_nodes / 1024.0 / 1024.0 << " MB" << endl;
			cout << "[zdtree tot leaf memory usage]: " << 1.0 * cur_stat.mem_leaf_nodes / 1024.0 / 1024.0 << " MB" << endl;
			cout << "[zdtree increased nodes]: " << cur_stat.num_inte_nodes - stat.num_inte_nodes << " interior nodes, " <<
													cur_stat.num_leaf_nodes - stat.num_leaf_nodes << " leaf nodes" << endl;
			stat = cur_stat;

			// cout << "[memory usage for inte nodes]: " << 1.0 * stat.mem_inte_nodes / 1024.0 / 1024.0 << " MB" << endl 
			//  	<< "[memory usage for leaf nodes]: " << 1.0 * stat.mem_leaf_nodes / 1024.0 / 1024.0 << " MB"  << endl;	
		}
	}

    template<class PT>
    void build_test(PT P){
		// auto node_num = 0;
		// unordered_map<int, int> leaf_map = {};
		mvq::Tree zdtree(mvq::Config::get().leaf_size);

        auto zdtree_build_avg = time_loop(
            3, 1.0, [&]() {
				zdtree.clear();
			},
            [&]() {
                auto P_set = get_sorted_points(P);
                zdtree.build(P_set);
				zdtree.multi_version_roots.emplace_back(zdtree.root);
				// node_num = zdtree.num_of_nodes();
				// leaf_map = zdtree.leaf_size_status();
            },
        [&](){
		});
		auto stat = zdtree.get_tree_statistics();
		// auto node_num = zdtree.num_of_nodes();

		// cout << "total # of nodes: " << node_num << endl;
		cout << "[zdtree stats]: " << endl <<
			"[# of inte nodes]: " << stat.num_inte_nodes << endl <<
			"[# of leaf nodes]: " << stat.num_leaf_nodes << endl <<
			"[memory usage for inte nodes]: " << 1.0 * stat.mem_inte_nodes / 1024.0 / 1024.0 << " MB" << endl <<
			"[memory usage for leaf nodes]: " << 1.0 * stat.mem_leaf_nodes / 1024.0 / 1024.0 << " MB"  << endl;
        cout << fixed << setprecision(6) << "[zdtree]: build time (avg): " << zdtree_build_avg << endl;
    }






        
    //     shared_ptr<mvq::BaseNode> new_ver = nullptr;

    //     P_set = get_sorted_address(P_delete);   //  delete P_delete
    //     new_ver = zdtree.multi_version_batch_delete_sorted(P_set, zdtree.root);





    //     set<int> delete_id1 = {0, 1, 2, 3, 4, 5}, update_id1 = {6, 7, 8, 9, 10}, insert_id1 = {n, n + 1, n + 2, n + 3, n + 4, n + 5};
    //     set<int> delete_id2 = {0, 1, 2, 6, 7, 8}, update_id2 = {3, 4, 5, 9, 11}, insert_id2 = {n, n + 1, n + 2, n + 6, n + 7, n + 8};

    //     parlay::sequence<Point> P_insert1 = {}, P_delete1 = {}, P_update1 = {};
    //     parlay::sequence<Point> P_insert2 = {}, P_delete2 = {}, P_update2 = {};


	// 	for (auto &id: delete_id1) P_delete1.emplace_back(P[id]);
	// 	for (auto &id: delete_id2) P_delete2.emplace_back(P[id]);

	// 	for (auto &id: update_id1) P_update1.emplace_back(P[id]);
	// 	for (auto &id: update_id2) P_update2.emplace_back(P[id]);

    //     shared_ptr<mvq::BaseNode> new_ver = nullptr;
    //     shared_ptr<mvq::BaseNode> ver1 = nullptr;
    //     shared_ptr<mvq::BaseNode> ver2 = nullptr;



	// 	for (auto &pt: P_update1) {cout << pt.x << ", " << pt.y << " ";} cout << endl;
	// 	for (auto &pt: P_update2) {cout << pt.x << ", " << pt.y << " ";} cout << endl;

	// 	auto [conflict_insert, conflict_update, conflict_delate] = zdtree.merge(zdtree.root, ver1, ver2);

	// 	cout << "insert conflict: "; for (auto &pt: conflict_insert) {cout << pt->id << " ";} cout << endl;
	// 	cout << "update conflict: "; for (auto &pt: conflict_update) {cout << pt->id << " ";} cout << endl;
	// 	cout << "delate conflict: "; for (auto &pt: conflict_delate) {cout << pt->id << " ";} cout << endl;
		
    // }




		

	// 	auto P_set = get_sorted_address(P);
	// 	zdtree.build(P_set);	//	original tree;
	// 	seq_zdtree.build(P_set);
	// 	cout << "original hash: " << zdtree.tree_hash(zdtree.root) << ", " << seq_zdtree.tree_hash(seq_zdtree.root) << endl;




	// 	cout << "after hash (par v.s. seq): " << zdtree.tree_hash(new_ver) << " " << seq_zdtree.tree_hash(seq_new_ver) << endl;
	// 	cout << fixed << setprecision(6) << "insert time (par v.s. seq): " << zdtree_time_cost << " " << seq_zdtree_time_cost << endl;


	// 	seq_zdtree.multi_version_roots.clear();
	// 	seq_zdtree.multi_version_roots.emplace_back(seq_zdtree.root);
	// 	seq_zdtree.multi_version_roots.emplace_back(seq_new_ver);
	// 	cout << "zdtree tot node count (two-version par v.s. seq): " << zdtree.node_count() << " " << seq_zdtree.node_count() << endl;
		
	// 	cout << endl;
	// 	// P_set = get_sorted_address(P);



    template<typename PT>
    void batch_insert_test(PT P, parlay::sequence<size_t> &batch_sizes){
	    auto n = P.size();

	    mvq::Tree zdtree(mvq::Config::get().leaf_size);
	    auto P_set = get_sorted_points(P); // build original tree

	    zdtree.build(P_set);

		auto rand_p = shuffle_point(P);

		shared_ptr<mvq::BaseNode> new_ver = nullptr;
		
		for (auto &num_processed: batch_sizes){
			if (num_processed > P.size()) num_processed = P.size();
	    	auto P2 = rand_p.substr(0, num_processed);	// first x%
	    	parlay::parallel_for (0, P2.size(), [&](int i){
		    	P2[i].id = n + i;
	    	});

			bool print_flag = true;
	    	auto zdtree_insert_avg = time_loop(
		    	3, 1.0, [&]() {
					new_ver.reset();
				},
		    	[&]() {
					// parlay::internal::timer t("debug", true);
			    	P_set = get_sorted_points(P2);
					// t.next("sort time");
			    	new_ver = zdtree.multi_version_batch_insert_sorted(P_set, zdtree.root);
					// t.next("insert time");
		    	},
	    	[&](){
				if (print_flag){
					cout << "# of points: " << zdtree.collect_records(new_ver).size() << endl;
					print_flag = false;
				}
			} );

			cout << "[batch_size]: " << num_processed << endl;		
	    	cout << fixed << setprecision(6) << "[zdtree]: batch insert time (avg): " << zdtree_insert_avg << endl;
		}
    }

    template<typename PT>
    void batch_delete_test(PT P, parlay::sequence<size_t> &batch_sizes){
	    mvq::Tree zdtree(mvq::Config::get().leaf_size);
	    auto P_set = get_sorted_points(P); // build original tree
	    zdtree.build(P_set);
	    shared_ptr<mvq::BaseNode> new_ver = nullptr;

		auto rand_p = shuffle_point(P);
		// auto rand_p = P;
		

		for (auto &num_processed: batch_sizes){
			if (num_processed > P.size()) num_processed = P.size();
	    	auto P2 = rand_p.substr(0, num_processed);	// first x%
			// bool print_flag = true;
	    	auto zdtree_delete_avg = time_loop(
		    	3, 1.0, [&]() {
					new_ver.reset();
				},
		    	[&]() {
					// parlay::internal::timer t("debug", true);
			    	// P_set = get_sorted_address(P2);
			    	P_set = get_sorted_points(P2);
					// t.next("sort time");
			    	new_ver = zdtree.multi_version_batch_delete_sorted(P_set, zdtree.root);
					// t.next("delete time");
		    	},
	    	[&](){
			} );
			cout << "# of points: " << zdtree.collect_records(new_ver).size() << endl;
			cout << "[batch_size]: " << num_processed << endl;
	    	// cout << fixed << setprecision(6) << "[leaf copy time]: " << zd_leaf_copy_time << endl;
	    	// cout << fixed << setprecision(6) << "[inte copy time]: " << zd_inte_copy_time << endl;
	    	cout << fixed << setprecision(6) << "[zdtree]: batch delete time (avg): " << zdtree_delete_avg << endl;
		}
    }


	
	//     mvq::Tree zdtree(mvq::Config::get().leaf_size);
	//     // auto P_set = get_sorted_address(P); // build original tree
	//     auto P_set = get_sorted_points(P); // build original tree




	//     // P_set = get_sorted_address(P3);
	//     P_set = get_sorted_points(P3);
	//     zdtree.build(P_set); // build for 1+x%
	//     parlay::hashtable<parlay::hash_numeric<int> > table(P2.size(), parlay::hash_numeric<int>{});




	//     cout << fixed << setprecision(6) << "[zdtree]: delete time (avg): " << zdtree_delete_avg << endl;
    // }

    template<typename PT>
    void batch_insert_test(PT P, PT P_insert, PT P_delete){
	    mvq::Tree zdtree(mvq::Config::get().leaf_size);
	    // auto P_set = get_sorted_address(P); // build original tree
	    auto P_set = get_sorted_points(P); // build original tree
	    zdtree.build(P_set);
		zdtree.multi_version_roots.emplace_back(zdtree.root);

		auto tot_node_count = zdtree.num_of_nodes();
		cout << "[initial node count]: " << tot_node_count << endl;

	    shared_ptr<mvq::BaseNode> new_ver = nullptr;
	    auto zdtree_insert_avg = time_loop(
		    3, 1.0, [&]() {},
		    [&]() {
			    P_set = get_sorted_points(P_insert);
			    new_ver = zdtree.multi_version_batch_insert_sorted(P_set, zdtree.root);
		    },
	    [&](){} );

		zdtree.multi_version_roots.emplace_back(new_ver);

	    cout << fixed << setprecision(6) << "[zdtree]: insert time (avg): " << zdtree_insert_avg << endl;

		tot_node_count = zdtree.num_of_nodes();
		cout << "[after_insertion]: " << tot_node_count << endl;
		
	    // zdtree.clear();
		// P_set = get_sorted_points(P);
	    // zdtree.build(P_set); 

	    parlay::hashtable<parlay::hash_numeric<int> > table(P_delete.size(), parlay::hash_numeric<int>{});

		auto inserted_ver = new_ver;
	    auto zdtree_delete_avg = time_loop(
		    3, 1.0, [&]() { },
		    [&]() {
			    P_set = get_sorted_points(P_delete);
			    // new_ver = zdtree.multi_version_batch_delete_sorted(P_set, zdtree.root);
			    new_ver = zdtree.multi_version_batch_delete_sorted(P_set, inserted_ver);
		    },
	    [&](){} );

		zdtree.multi_version_roots.emplace_back(new_ver);		


	    cout << fixed << setprecision(6) << "[zdtree]: delete time (avg): " << zdtree_delete_avg << endl;

		tot_node_count = zdtree.num_of_nodes();
		cout << "[after_deletion]: " << tot_node_count << endl;

		zdtree.multi_version_roots[1].reset();
		tot_node_count = zdtree.num_of_nodes();
		cout << "[after_remove_intermediate]: " << tot_node_count << endl;

    }


	template<class PT, class RQ>
	void range_count_test(PT P, RQ querys, parlay::sequence<size_t> &cnt){
	    size_t n = P.size();
	    auto P2 = P.substr(0, P.size() / 2);   // insert a half
	    for (size_t i = 0; i < P2.size(); i++){
		    P2[i].id = n + i;
	    }
	    mvq::Tree zdtree(mvq::Config::get().leaf_size);

	    auto P_set = get_sorted_points(P);
	    zdtree.build(P_set);	// 	build for P

		parlay::sequence<size_t> rangeCnt(querys.size());
		for (size_t i = 0; i < querys.size(); i++){
			auto avg_time = time_loop(
				3, 1.0, 
				[&]() {},
				[&]() {					
					rangeCnt[i] = zdtree.range_count(querys[i], mvq::Config::get().largest_mbr);
				},
				[&](){} );
			if (rangeCnt[i] != cnt[i]){
				cout << "[ERROR] Incorrect" << endl;
			}
			else{
				cout << fixed << setprecision(6) << rangeCnt[i] << " " << avg_time << endl;
			}
		}
		


		// ofstream regionCntOut("zd_range_count.txt");
		
	}

	template<class PT, class RQ>
	void range_report_test(PT P, RQ querys, parlay::sequence<size_t> &cnt){
	    mvq::Tree zdtree(mvq::Config::get().leaf_size);
	    auto P_set = get_sorted_points(P);
	    zdtree.build(P_set);	// 	build for P

		parlay::sequence<parlay::sequence<Point> > rangeReport(querys.size());
		parlay::sequence<size_t> rangeReportCnt(querys.size(), 0);

		for (size_t i = 0; i < querys.size(); i++){
			rangeReport[i].resize(cnt[i]);	
		}

		for (size_t i = 0; i < querys.size(); i++){
			// print_mbr(querys[i]);
			auto avg_time = time_loop(
				3, 1.0, 
				[&]() {
					rangeReportCnt[i] = 0;
				},
				[&]() {					
					zdtree.range_report(querys[i], mvq::Config::get().largest_mbr, rangeReportCnt[i], rangeReport[i]);
				},
				[&](){} );
			if (rangeReportCnt[i] != cnt[i]){
				cout << "[ERROR] Incorrect" << endl;
				cout << rangeReportCnt[i] << " " << cnt[i] << endl;
			}
			else{
				cout << fixed << setprecision(6) << rangeReportCnt[i] << " " << avg_time << endl;
			}
		}


		// cout << fixed << setprecision(6) << "[zdtree] range report time (avg): " << rangeReport_avg << endl;


	}


	template<class PT>
	void knn_test(PT P, size_t k = 10, size_t q_num = 50000){
	    mvq::Tree zdtree(mvq::Config::get().leaf_size);
	    // auto P_set = get_sorted_address(P);
	    auto P_set = get_sorted_points(P);
	    zdtree.build(P_set);	// 	build for P

		parlay::sequence<FT> knn_sqrdis(q_num);

		auto knnReport_avg = time_loop(
			3, 1.0, [&]() {},
			[&]() {
				for (size_t i = 0; i < q_num; i++){
					knn_sqrdis[i] = zdtree.knn_report(k, P[i], mvq::Config::get().largest_mbr).top().second;
				}
			},
		[&](){} );
		cout << fixed << setprecision(6) << "[zdtree] knn report time (avg): " << knnReport_avg << endl;

		ofstream knnReportOut("zd-knn.res");
		for (size_t i = 0; i < q_num; i++){
			knnReportOut << knn_sqrdis[i] << endl;
			// size_t id = i < 5 ? i : n / 2 + i;
			// knnReportOut << knn_sqrdis[i] << ", brute force: " << geobase::knn_bf(k, *P_set[i], P) << endl;
			// knnReportOut << knn_sqrdis[id] << ", brute force: " << knn_bf(k, P_set[id], P_set) << endl;
		}
		// cout << "knn finished." << endl;
	}


    template<class PT, class RQ>
    void zdtree_test(PT P, RQ querys, int q_type){
	    size_t n = P.size();
	    auto P2 = P.substr(0, P.size() / 2);   // insert a half
	    for (size_t i = 0; i < P2.size(); i++){
		    P2[i].id = n + i;
	    }
	    mvq::Tree zdtree(mvq::Config::get().leaf_size);
	    auto P_set = get_sorted_address(P);
	    zdtree.build(P_set);	// 	build for P

	    // auto phash = zdtree.tree_hash(zdtree.root);	// hash value for construct P directly
	    // cout << "build zdtree finished." << endl;

	    /* insertion test */

	    /* in order traverse */
	    // zdtree.in_order_traverse(zdtree.root, 0, 0, 32, true);
	    // zdtree.in_order_traverse(zdtree.root);
	    // cout << "traverse finished." << endl;

	    /* range count test */
	    // if (q_type == 4){
		    parlay::sequence<size_t> rangeCnt(querys.size());
		    auto rangeCnt_avg = time_loop(
			    3, 1.0, [&]() {},
        	    [&]() {
				    parlay::parallel_for(
    	        	    0, querys.size(),
					    [&]( size_t i ) {
						    // rangeCnt[i] = zdtree.range_count(querys[i]);
						    rangeCnt[i] = zdtree.range_count(querys[i], mvq::Config::get().largest_mbr);
            	    });
			    },
		    [&](){} );
		    cout << fixed << setprecision(6) << "zdtree range count time (avg): " << rangeCnt_avg << endl;
		    ofstream regionCntOut("zd_range_count.txt");
		    for (size_t i = 0; i < querys.size(); i++){
			    regionCntOut << rangeCnt[i] << endl;
		    }
			// cout << "range count finished." << endl;
	    // }

	    /* range report test */
	    // if (q_type == 8){
		    parlay::sequence<parlay::sequence<Point*> > rangeReport(querys.size());
		    parlay::sequence<size_t> rangeReportCnt(querys.size(), 0);
		    for (size_t i = 0; i < querys.size(); i++){
			    rangeReport[i].resize(mvq::Config::get().maxSize);	
		    }
		    auto rangeReport_avg = time_loop(
			    3, 1.0, [&]() {},
        	    [&]() {
				    parlay::parallel_for(
    	        	    0, querys.size(),
					    [&]( size_t i ) {
						    rangeReportCnt[i] = 0;
						    // zdtree.range_report(querys[i], rangeReportCnt[i], rangeReport[i]);
						    zdtree.range_report(querys[i], mvq::Config::get().largest_mbr, rangeReportCnt[i], rangeReport[i]);
            	    });
			    },
		    [&](){} );
		    cout << fixed << setprecision(6) << "zdtree range report time (avg): " << rangeReport_avg << endl;
		    ofstream regionReportOut("zd_range_report.txt");
		    for (size_t i = 0; i < querys.size(); i++){
			    regionReportOut << rangeReportCnt[i] << endl;
		    }
			// cout << "range report finshied." << endl;
	    // }

	    /* knn report test */
	    // if (q_type == 1){
		    size_t k = 10;
		    P_set.resize(n);
		    parlay::parallel_for(0, n, [&](size_t i){
			    P_set[i] = i < P.size() ? &P[i] : &P2[i - P.size()];
		    });
		    size_t knn_test_size = n;
		    parlay::sequence<FT> knn_sqrdis(knn_test_size);
		    auto knnReport_avg = time_loop(
			    3, 1.0, [&]() {},
        	    [&]() {
				    parlay::parallel_for(
    	        	    0, knn_test_size,
					    [&]( size_t i ) {
						    // knn_sqrdis[i] = zdtree.knn_report(k, &P[i]).top().second;
						    // knn_sqrdis[i] = zdtree.knn_report(k, &P[i], mvq::Config::get().largest_mbr).top().second;
						    knn_sqrdis[i] = zdtree.knn_report(k, P_set[i], mvq::Config::get().largest_mbr).top().second;
            	    });
			    },
		    [&](){} );
		    cout << fixed << setprecision(6) << "zdtree knn report time (avg): " << knnReport_avg << endl;
		    ofstream knnReportOut("zd_knn_report.txt");
		    for (size_t i = 0; i < 10; i++){
			    size_t id = i < 5 ? i : n / 2 + i;
			    // knnReportOut << knn_sqrdis[i] << ", brute force: " << geobase::knn_bf(k, *P_set[i], P) << endl;
			    knnReportOut << knn_sqrdis[id] << ", brute force: " << knn_bf(k, *P_set[id], P_set) << endl;
		    }
			cout << "knn finished." << endl;
	    // }
    }

    template<class PT, class RQ>
    void CPAMZ_test(PT &P, RQ querys){
	    auto CPAMZ = CPAMZ::CPAMZ_init(P);
    #ifdef DEBUG
	    cout << "CPAMZ size = " << CPAMZ.size() << endl;
    #endif
	    parlay::internal::timer t;
	    parlay::sequence<Point> range_query_res;
	    for (auto query_mbr: querys){
		    range_query_res = CPAMZ::range_report(CPAMZ, query_mbr);
	    }
	    cout << "Z-CPAM query time (avg): " << t.next_time() << endl;
    #ifdef DEBUG
	    // auto range_query_res = CPAMZ::range_report(CPAMZ[0], query_mbr);
	    cout << "[Z-CPAM RESULT] range query results: " << range_query_res.size() << "| ";
	    for (auto p: range_query_res){
		    cout << "(" << p.x << ", " << p.y << ")" << " ";
	    }
	    cout << endl;	
    #endif
    }

    template<class PT, class RQ>
    void zMAP_test(PT &P, RQ querys){
	    auto zMAP = CPAMZ::zMAP_init(P);
    #ifdef DEBUG
	    cout << "zMAP size = " << zMAP.size() << endl;
    #endif
	    parlay::internal::timer t;
	    parlay::sequence<Point> range_query_res;
	    for (auto query_mbr: querys){
		    range_query_res = CPAMZ::zMAP_range_report(zMAP, query_mbr);
	    }
	    cout << "zMAP query time (avg): " << t.next_time() << endl;
	    // auto range_query_res = CPAMZ::range_report(CPAMZ[0], query_mbr);
    #ifdef DEBUG
	    cout << "[Z-Map RESULT] range query results: " << range_query_res.size() << "| ";
	    for (auto p: range_query_res){
		    cout << "(" << p.x << ", " << p.y << ")" << " ";
	    }
	    cout << endl;	
    #endif
    }
}