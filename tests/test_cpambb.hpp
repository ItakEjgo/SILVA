#pragma once
#include "test_utils.hpp"
namespace CPAMBB{


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

			if (P_update[i].size() != P_updove[i].size()){
				cout << "[ERROR]: inconsistent # of update pts!" << endl;
			}

			P_delete[i].append(P_updove[i]);
			P_insert[i].append(P_update[i]);

			cur_year += 1;
			delete_mbr = insert_mbr; // useless, just remove warning
			update_mbr = updove_mbr;
		}

		vector<CPAMBB::zmap> all_versions;
		CPAMBB::zmap tree;

		/* Build initial version */
		auto build_avg = time_loop(
            3, 1.0, [&]() {
				tree.clear();
			},
            [&]() {
				tree = CPAMBB::map_init(P);	// initi
            },
   	    	[&](){
			});

		cout << "[cpambb init build time]: " << fixed << setprecision(6) << build_avg << " Seconds" << endl;

		all_versions.emplace_back(tree);

		auto f_noop = [&](const auto &et){	return 0; };

		std::unordered_map<size_t, bool> mmp = {}, num_mmp = {};
		double cur_mem = 0.0, prev_mem = 0.0; 
		size_t cur_inte_num = 0, cur_leaf_num = 0, cur_leaf_sz = 0;
		size_t pre_inte_num = 0, pre_leaf_num = 0, pre_leaf_sz = 0;
		cur_leaf_sz += pre_leaf_sz;
		
		prev_mem = 1.0 * tree.size_in_bytes(f_noop, mmp);
		tie(pre_inte_num, pre_leaf_num, pre_leaf_sz) = tree.node_stats(num_mmp);

		cout << "[init-version memory]: " << prev_mem / 1024.0 / 1024.0 << " MB" << endl;
		cout << "[init-version node nums]: " << pre_inte_num << " interior nodes, " << pre_leaf_num << " leaf nodes" << endl;
		
		vector<zmap> new_ver(version_num);

		for (auto i = 0; i < version_num; i++){
			cout << "dealing with version " << i + 1 << ":" << endl;

			auto commit_avg = time_loop(
				3, 1.0, 
				[&]() {
					new_ver[i].clear();
				},
				[&]() {
					new_ver[i] = CPAMBB::map_commit(all_versions[i], P_insert[i], P_delete[i]);
				},
				[&](){});

			all_versions.emplace_back(new_ver[i]);
			cur_mem = 0;
			mmp.clear();
			num_mmp.clear();
			cur_inte_num = cur_leaf_num = 0;

			for (size_t j = 0; j < all_versions.size(); j++){
				cur_mem += 1.0 * all_versions[j].size_in_bytes(f_noop, mmp); 	// accumulate all version memories, shared pointers only count once
				auto [tmp_inte_num, tmp_leaf_num, tmp_leaf_sz] = all_versions[j].node_stats(num_mmp);
				cur_inte_num += tmp_inte_num, cur_leaf_num += tmp_leaf_num;
			}
			cout << "[new ver commit time]: " << fixed << setprecision(6) << commit_avg << " Seconds" << endl;
			cout << "[cpambb memory usage]: " << (cur_mem - prev_mem) / 1024.0 / 1024.0  << " MB" << endl;
			cout << "[cpambb node nums]: " << cur_inte_num - pre_inte_num << " interior nodes, " << cur_leaf_num - pre_leaf_num << " leaf nodes" << endl;
			prev_mem = cur_mem;
			pre_inte_num = cur_inte_num;
			pre_leaf_num = cur_leaf_num;
		}
	}

	template<typename PT, typename RQ>
	void plain_spatial_diff_test_latency(PT &P, RQ &range_queries, parlay::sequence<size_t> &batch_sizes, size_t &insert_ratio){
		/*  build tree */
		auto cpambb0 = CPAMBB::map_init(P);	//	initial version
		auto max_batch_size = batch_sizes[batch_sizes.size() - 1];
		/* get insert, delete points */
		auto P_test = geobase::shuffle_point(P, max_batch_size);
		auto [P_insert_set, P_delete_set] = geobase::split_insert_delete(P_test, insert_ratio, P.size());

		for (auto &batch_size: batch_sizes){
			cout << "[INFO] Batch Size: " << batch_size << endl;
			auto insert_num = batch_size / 10 * insert_ratio;
			auto delete_num = batch_size / 10 * (10 - insert_ratio);

			auto P_insert = P_insert_set.substr(0, insert_num);
			auto P_delete = P_delete_set.substr(0, delete_num);

			auto P_newver = geobase::collect_newver_point(P, P_insert, P_delete);
			auto cpambb1 = CPAMBB::map_init(P_newver);

			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());

			auto l_pts = parlay::sequence<Point>::uninitialized(2 * mvq::Config::get().maxSize);
			auto r_pts = parlay::sequence<Point>::uninitialized(2 * mvq::Config::get().maxSize);

			for (size_t i = 0; i < range_queries.size(); i++){
				auto avg_time = time_loop(
					3, 1.0,
					[&](){},
					[&](){
						diff_type ret_diff(mvq::Config::get().maxSize, mvq::Config::get().maxSize);
						CPAMBB::plain_map_spatial_diff(cpambb0, cpambb1, range_queries[i], ret_diff, l_pts, r_pts);
						ret_diff.compact();
						addCnt[i] = ret_diff.add.size();
						removeCnt[i] = ret_diff.remove.size();
					},
					[&]{}
				);
				cout << fixed << setprecision(6) << i << " " << avg_time << endl;
			}

			#ifdef TEST
				string file_name = "output/cpambb_spatial_diff_plain-" + to_string(batch_size); 
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
		auto cpambb0 = CPAMBB::map_init(P);	//	initial version
		auto max_batch_size = batch_sizes[batch_sizes.size() - 1];

        /* get insert, delete points */
        auto P_test = geobase::shuffle_point(P, max_batch_size);
        auto [P_insert_set, P_delete_set] = geobase::split_insert_delete(P_test, insert_ratio, P.size());

		for (auto &batch_size: batch_sizes){
			cout << "[INFO] Batch Size: " << batch_size << endl;
			auto insert_num = batch_size / 10 * insert_ratio;
			auto delete_num = batch_size / 10 * (10 - insert_ratio);

			auto P_insert = P_insert_set.substr(0, insert_num);
			auto P_delete = P_delete_set.substr(0, delete_num);

			auto cpambb1 = CPAMBB::map_delete(P_delete, cpambb0); 
			auto cpambb2 = CPAMBB::map_insert(P_insert, cpambb1);	//	new	version
        
        	parlay::sequence<size_t> addCnt(range_queries.size());
        	parlay::sequence<size_t> removeCnt(range_queries.size());

			auto l_pts = parlay::sequence<Point>::uninitialized(2 * mvq::Config::get().maxSize);
			auto r_pts = parlay::sequence<Point>::uninitialized(2 * mvq::Config::get().maxSize);
			for (size_t i = 0; i < range_queries.size(); i++){
				auto avg_time = time_loop(
					3, 1.0,
					[&](){},
					[&](){
						// for (size_t i = 0; i < range_queries.size(); i++){
							diff_type ret_diff(mvq::Config::get().maxSize, mvq::Config::get().maxSize);
							CPAMBB::plain_map_spatial_diff(cpambb0, cpambb2, range_queries[i], ret_diff, l_pts, r_pts);
							ret_diff.compact();
							addCnt[i] = ret_diff.add.size();
							removeCnt[i] = ret_diff.remove.size();
						// }
					},
					[&]{}
				);
				cout << fixed << setprecision(6) << i << " " << avg_time << endl;
			}
			
        
			#ifdef TEST
        		string file_name = "output/cpambb_spatial_diff-" + to_string(batch_size); 
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
		auto cpambb0 = CPAMBB::map_init(P);	//	initial version

		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) break;
			cout << "[ratio (insert/delete)]: " << batch_size << "/" << 10 - batch_size << endl;	

			auto insert_num = P.size() / 100 * batch_size;
			auto delete_num = P.size() / 100 * (10 - batch_size);

			auto P_insert = P.substr(0, insert_num);
			auto P_delete = P.substr(P.size() - delete_num, delete_num);
	
			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += P.size();
			});

			cout << "# of insertion/deletion: " << P_insert.size() << ", " << P_delete.size();
	
			auto cpambb1 = CPAMBB::map_delete(P_delete, cpambb0); 
			auto cpambb2 = CPAMBB::map_insert(P_insert, cpambb1);	//	new	version
	
			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());

			auto diff_avg = time_loop(
				3, 1.0, [&]() {},
				[&]() {
					parlay::parallel_for(0, range_queries.size(), [&](int i){
						diff_type ret_diff(mvq::Config::get().maxSize, mvq::Config::get().maxSize);
						CPAMBB::map_spatial_diff(cpambb0, cpambb2, range_queries[i], ret_diff);
						ret_diff.compact();
						addCnt[i] = ret_diff.add.size();
						removeCnt[i] = ret_diff.remove.size();
					});
				},
			[&](){} );
			cout << fixed << setprecision(6) << "[CPAMBB] spatial-diff time (avg): " << diff_avg << endl;
		}
	}
	

	/*	50% insertion, 50% deletion	*/
	template<typename PT, typename RQ>
	auto spatial_diff_test_fix_ratio(PT P,  RQ &range_queries, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto cpambb0 = CPAMBB::map_init(P);	//	initial version

		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) break;
			cout << "[batch-size]: " << batch_size << endl;
	
			auto P_insert = P.substr(0, batch_size / 2);
			auto P_delete = P.substr(P.size() - batch_size / 2, batch_size / 2);
	
			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += P.size();
			});

			cout << "# of insertion/deletion: " << P_insert.size() << ", " << P_delete.size();
	
			auto cpambb1 = CPAMBB::map_delete(P_delete, cpambb0); 
			auto cpambb2 = CPAMBB::map_insert(P_insert, cpambb1);	//	new	version
	
			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());

			auto diff_avg = time_loop(
				3, 1.0, [&]() {},
				[&]() {
					parlay::parallel_for(0, range_queries.size(), [&](int i){
						diff_type ret_diff(mvq::Config::get().maxSize, mvq::Config::get().maxSize);
						CPAMBB::map_spatial_diff(cpambb0, cpambb2, range_queries[i], ret_diff);
						ret_diff.compact();
						addCnt[i] = ret_diff.add.size();
						removeCnt[i] = ret_diff.remove.size();
					});
				},
			[&](){} );
			cout << fixed << setprecision(6) << "[CPAMBB] spatial-diff time (avg): " << diff_avg << endl;
		}
	}

	template<typename PT, typename RQ>
	auto spatial_diff_test(PT P, RQ &range_queries, parlay::sequence<size_t> &batch_sizes, bool &early_end, bool use_hilbert = false){
		auto cpambb0 = CPAMBB::map_init(P);	//	initial version

		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) batch_size = P.size();
			cout << "[batch-size]: " << batch_size << endl;

			auto P_insert = P.substr(0, batch_size);
			auto P_delete = P.substr(P.size() - batch_size, batch_size);

			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += P.size();
			});

			auto cpambb1 = CPAMBB::map_delete(P_delete, cpambb0); 
			auto cpambb2 = CPAMBB::map_insert(P_insert, cpambb0);	//	new	version

			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());
			// auto addCnt = 0;
			// auto removeCnt = 0;
			// parlay::sequence<Point> ret;



			if (!early_end){
				decltype(cpambb0) commit_ver;
	    		auto commit_avg = time_loop(
		    		3, 1.0, [&]() {
						commit_ver.clear();
					},
		    		[&]() {
						commit_ver = CPAMBB::map_commit(cpambb0, P_insert, P_delete);
		    		},
	    			[&](){} 
				);

				parlay::sequence<Point> conflict_insert, conflict_update, conflict_delete;
				decltype(cpambb0) merge_ver;
				auto merge_avg = time_loop(
		    		3, 1.0, [&]() {
						merge_ver.clear();
						conflict_insert.clear();
						conflict_update.clear();
						conflict_delete.clear();
					},
		    		[&]() {
						tie(merge_ver, conflict_insert, conflict_update, conflict_delete) = CPAMBB::map_merge(cpambb0, cpambb1, cpambb2);
		    		},
	    			[&](){} 
				);
				cout << "[INFO] commit, merge size: " << commit_ver.size() << ", " << merge_ver.size() << endl; 
				cout << fixed << setprecision(6) << "[CPAMBB]: spatial-commit time (avg): " << commit_avg << endl;
				cout << fixed << setprecision(6) << "[CPAMBB]: spatial-merge time (avg): " << merge_avg << endl;
			}

		}

		// return ret;
	}


	template<typename PT>
	void multi_version_query_test(PT P, string query_dir, int batch_percent = 10, int version_num = 6){
		// build zdtree initial version
		auto CPAMZ = CPAMBB::map_init(P);
		cout << "build finished" << endl;

		auto num_insert_version = version_num / 2;
		auto num_delete_version = version_num - num_insert_version;
		auto batch_size = P.size() * batch_percent / 100;
		

		parlay::sequence<decltype(CPAMZ)> all_versions(7);
		all_versions[0] = CPAMZ;
		// insert 3 versions
		for (auto i = 0; i != num_insert_version; i++){
			auto P_insert = P.substr(i * batch_size, batch_size);

			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += (i + 1) * P.size();
			});

			all_versions[i + 1] = CPAMBB::map_insert(P_insert, all_versions[i]); 
		}
		cout << "insert finished" << endl;	
		//	delete 3 versions
		for (auto i = 0; i != num_delete_version; i++){
			auto P_delete = P.substr(i * batch_size, batch_size);

			parlay::parallel_for(0, P_delete.size(), [&](size_t j){
				P_delete[j].id += (i + 1) * P.size();
			});
			
			all_versions[i + 4] = CPAMBB::map_delete(P_delete, all_versions[i + 3]); 
			// all_versions.push_back(new_version);
		}
		cout << "delete finished" << endl;	
		
		// range count query test
		for (auto i = 0; i != 1; i++){	//	0, 1, 2 represent small, median, and large regions
			auto cur_query_dir = query_dir + "/1.in-" + to_string(i) + ".qry";
			auto [cnt, range_count_querys] = read_range_query(cur_query_dir, 8, mvq::Config::get().maxSize);

			for (size_t j = 0; j < all_versions.size(); j++){

				parlay::sequence<size_t> rangeCnt(range_count_querys.size());
				auto rangeCnt_avg = time_loop(
					3, 1.0, [&]() {},
					[&]() {
						parlay::parallel_for(
							0, range_count_querys.size(),
							[&]( size_t k ) {
								rangeCnt[k] = CPAMBB::range_count(all_versions[j], range_count_querys[k]); 
						});
					},
				[&](){} );

				cout << fixed << setprecision(6) << "CPAM-BB range count time (avg) for region " << i << " on version " << j << ": " << rangeCnt_avg << endl;
				auto output_name = "CPAM-BB_range_count-" + to_string(i) + "-on-" + to_string(j) + ".txt";
				ofstream regionCntOut(output_name);
				for (size_t k = 0; k < range_count_querys.size(); k++){
					regionCntOut << rangeCnt[k] << endl;
				}
			}
		}

		// range report query test
		for (auto i = 0; i != 1; i++){	//	0, 1, 2 represent small, median, and large regions
			auto cur_query_dir = query_dir + "/1.in-" + to_string(i) + ".qry";
			auto [cnt, range_report_querys] = read_range_query(cur_query_dir, 8, mvq::Config::get().maxSize);

			for (size_t j = 0; j < all_versions.size(); j++){

				parlay::sequence<size_t> rangeCnt(range_report_querys.size());
				auto rangeReport_avg = time_loop(
					3, 1.0, [&]() {},
					[&]() {
						parlay::parallel_for(
							0, range_report_querys.size(),
							[&]( size_t k ) {
		    					rangeCnt[k] = CPAMBB::range_report(all_versions[j], range_report_querys[k]).size();
						});
					},
				[&](){} );

				cout << fixed << setprecision(6) << "CPAM-BB range report time (avg) for region " << i << " on version " << j << ": " << rangeReport_avg << endl;
				auto output_name = "CPAM-BB_range_report-" + to_string(i) + "-on-" + to_string(j) + ".txt";
				ofstream regionCntOut(output_name);
				for (size_t k = 0; k < range_report_querys.size(); k++){
					regionCntOut << rangeCnt[k] << endl;
				}
			}
		}
	}

	template<typename PT>
	void diff_test(PT P, int batch_percent = 10, bool use_hilbert = false){
		auto cpambb0 = CPAMBB::map_init(P);	//	initial version

		auto batch_size = P.size() * batch_percent / 100;	//	insertion 10%
		auto P_insert = P.substr(0, batch_size);
		parlay::parallel_for(0, P_insert.size(), [&](size_t j){
			P_insert[j].id += P.size();
		});

		auto P_delete = P.substr(0, 2 * batch_size);

		auto cpambb1 = CPAMBB::map_insert(P_insert, cpambb0);	//	new	version
		cpambb1 = CPAMBB::map_delete(P_delete, cpambb1); 

		auto add_sz = 0, remove_sz = 0;
	    auto cpam_diff_avg = time_loop(
		    3, 1.0, [&]() {},
		    [&]() {
				auto [add, remove] = CPAMBB::map_diff(cpambb0, cpambb1);
				add_sz = add.size();
				remove_sz = remove.size();
		    },
	    [&](){} );
		

		cout << "add size: " << add_sz << endl;
		cout << "remove size: " << remove_sz << endl;
		if (use_hilbert) cout << "[Hilbert-CPAMBB]: ";
		else cout << "[Zorder-CPAMBB]: ";
		cout << fixed << setprecision(6) << "CPAMBB diff time (avg): " << cpam_diff_avg << endl;
	}


	template<typename PT>
	void build_test(PT P, bool use_hilbert = false){
		CPAMBB::zmap tree;

		auto cpam_build_avg = time_loop(
			3, 1.0, [&](){
				tree.clear();
			},
			[&](){
				tree = CPAMBB::map_init(P);
			},
		[&](){} );
		// auto [mem_inte_nodes, mem_leaf_nodes] = CPAMBB::size_in_bytes();
		auto [num_inte_nodes, num_leaf_nodes, leaf_size] = tree.node_stats();
		// cout << "leaf sz = " << 1.0 * leaf_size / 1024.0 / 1024.0 << " MB" << endl;

		auto f_noop = [&](const auto &et){
			return 0;
		};

		cout << "[cpambb memory usage]: " << endl <<
			"[# of inte nodes]: " << num_inte_nodes << endl << 
			"[# of leaf nodes]: " << num_leaf_nodes << endl <<
			"[tree size]: " << 1.0 * tree.size_in_bytes(f_noop) / 1024.0 / 1024.0 << " MB" << endl;
			// "[memory usage for inte nodes]: " << 1.0 * mem_inte_nodes / 1024.0 / 1024.0 << " MB" << endl <<
			// "[memory usage for leaf nodes]: " << 1.0 * mem_leaf_nodes / 1024.0 / 1024.0 << " MB"  << endl;

		// cout << "cpambb print stats: " << endl;
		// CPAMBB::print_stats();

		if (use_hilbert) cout << "[Hilbert-CPAMBB]: ";
		else cout << "[Zorder-CPAMBB]: ";
		cout << fixed << setprecision(6) << "build time (avg): " << cpam_build_avg << endl;
	}

	template<class PT, class RQ>
    void range_count_test(PT P, RQ querys, parlay::sequence<size_t> &cnt, bool use_hilbert = false){
	    auto tree = CPAMBB::map_init(P, use_hilbert);

		parlay::sequence<size_t> rangeCnt(querys.size());


		/* Latency Test */
		for (size_t i = 0; i < querys.size(); i++){
			auto avg_time = time_loop(
				3, 1.0, 
				[&]() {},
				[&]() {					
					rangeCnt[i] = CPAMBB::range_count(tree, querys[i], use_hilbert);
				},
				[&](){} );
			if (rangeCnt[i] != cnt[i]){
				cout << "[ERROR] Incorrect" << endl;
			}
			else{
				cout << fixed << setprecision(6) << rangeCnt[i] << " " << avg_time << endl;
			}
		}
	
    }


	template<class PT, class RQ>
    void range_report_test(PT P, RQ querys, parlay::sequence<size_t> &cnt, bool use_hilbert = false, size_t par_for_granularity = 100){
	    auto tree = CPAMBB::map_init(P, use_hilbert);

		// parlay::sequence<Bounding_Box> q2(querys.size());
		parlay::sequence<size_t> rangeCnt(querys.size());
		parlay::sequence<parlay::sequence<Point> > rangeReport(querys.size());
		for (size_t i = 0; i < querys.size(); i++){
			rangeReport[i].resize(cnt[i]);
		}

		// size_t tot_inte = 0, tot_leaf = 0; double tot_time = 0;
		// for (size_t i = 0; i < querys.size(); i++){
		auto avg_time = time_loop(
			3, 1.0, 
			[&]() {},
			[&]() {					
				// for (size_t i = 0; i < querys.size(); i++){
				parlay::parallel_for(0, querys.size(), [&](size_t i){
					rangeCnt[i] = CPAMBB::range_report(tree, querys[i], rangeReport[i], use_hilbert);
				});
				// }
			},
			[&](){} 
		);

		if (use_hilbert) cout << "[Hilbert-CPAMBB]: ";
		else cout << "[CPAMBB]: ";
		cout << fixed << setprecision(6) << "range report time (avg): " << avg_time << endl;


		/* range report sample test. */


		// }
		// cout << "[INFO] totoal visited inte: " << tot_inte << endl;
		// cout << "[INFO] totoal visited leaf: " << tot_leaf << endl;

		// size_t tot_inte, tot_leaf;


		// cout << "[tot visted inte]: " << tot_inte << endl;
		// cout << "[tot visted leaf]: " << tot_leaf << endl;

		// if (use_hilbert) cout << "[Hilbert-CPAMBB]: ";
		// else cout << "[Zorder-CPAMBB]: ";
		// cout << fixed << setprecision(6) << "range report time (avg): " << rangeReport_avg << endl;


    }

	template<class PT>
    void knn_test(PT P, size_t k = 10, size_t q_num = 50000, bool use_hilbert = false){
	    auto tree = CPAMBB::map_init(P, use_hilbert);
		
		auto knn_sqrdis = parlay::sequence<size_t>::uninitialized(q_num);
		
		auto avg_time = time_loop(
			3, 1.0, 
			[&]() {},
			[&]() {					
				for (size_t i = 0; i < q_num; i++){
					knn_sqrdis[i] = CPAMBB::knn(tree, P[i], k).top().second;
				}
			},
			[&](){} );
		cout << fixed << setprecision(6) << "[CPAMBB] KNN Latency: " << avg_time << endl;



		// if (use_hilbert) cout << "[Hilbert-CPAMBB]: ";
		// else cout << "[Zorder-CPAMBB]: ";
		// cout << fixed << setprecision(6) << "knn report time (avg): " << rangeReport_avg << endl;

		// /* Correctness Check */

    }

	

	template<typename PT>
	void batch_insert_test(PT P, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto n = P.size();
		auto m1 = CPAMBB::map_init(P, use_hilbert);	//	build original tree
		decltype(m1) m2;

		auto rand_p = shuffle_point(P);

		for (auto &num_processed: batch_sizes){
			if (num_processed > P.size()) num_processed = P.size();
			auto P2 = rand_p.substr(0, num_processed);

	    	parlay::parallel_for (0, P2.size(), [&](int i){
		    	P2[i].id = n + i;
	    	});
			bool print_flag = true;
	    	auto cpam_insert_avg = time_loop(
		    	3, 1.0, [&]() {
					m2.clear();
				},
		    	[&]() {
					m2 = CPAMBB::map_insert(P2, m1, use_hilbert);
		    	},
	    	[&](){
				if (print_flag){
					cout << "# of points: " << m2.size() << endl;
					print_flag = false;
				}			
			});

			cout << "[batch_size]: " << num_processed << endl;
			if (use_hilbert) cout << "[Hilbert-CPAMBB]: ";
			else cout << "[Zorder-CPAMBB]: ";
			cout << fixed << setprecision(6) << "batch insert time (avg): " << cpam_insert_avg << endl;
		}

	}

	template<typename PT>
	void batch_delete_test(PT P, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto m1 = CPAMBB::map_init(P, use_hilbert);	//	build original tree
		decltype(m1) m2;
		auto rand_p = shuffle_point(P);

		for (auto &num_processed: batch_sizes){
			if (num_processed > P.size()) num_processed = P.size();
			auto P2 = rand_p.substr(0, num_processed);
			bool print_flag = true;

	    	auto cpam_insert_avg = time_loop(
		    	3, 1.0, [&]() {
					m2.clear();
				},
		    	[&]() {
					m2 = CPAMBB::map_delete(P2, m1, use_hilbert);
		    	},
	    	[&](){
				if (print_flag){
					cout << "# of points: " << m2.size() << endl;
					print_flag = false;
				}
			} );

			cout << "[batch_size]: " << num_processed << endl;
			if (use_hilbert) cout << "[Hilbert-CPAMBB]: ";
			else cout << "[Zorder-CPAMBB]: ";
			cout << fixed << setprecision(6) << "batch delete time (avg): " << cpam_insert_avg << endl;
		}
	}


}

//	use CPAM (without bounding box) to build
