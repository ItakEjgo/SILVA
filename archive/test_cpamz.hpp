#pragma once
#include "test_utils.hpp"
namespace CPAMZ{

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

		vector<CPAMZ::zmap> all_versions;
		CPAMZ::zmap tree;

		/* Build initial version */
		auto build_avg = time_loop(
            3, 1.0, [&]() {
				tree.clear();
			},
            [&]() {
				tree = CPAMZ::CPAMZ_init(P);	// initi
            },
   	    	[&](){
			});

		cout << "[cpamz init build time]: " << fixed << setprecision(6) << build_avg << " Seconds" << endl;

		all_versions.emplace_back(tree);

		auto f_noop = [&](const auto &et){	return 0; };

		std::unordered_map<size_t, bool> mmp = {};
		double cur_mem = 0.0, prev_mem = 0.0; 
		prev_mem = 1.0 * tree.size_in_bytes(f_noop, mmp);

		cout << "[init-version memory]: " << prev_mem / 1024.0 / 1024.0 << " MB" << endl;
		
		vector<CPAMZ::zmap> new_ver(version_num);

		for (auto i = 0; i < version_num; i++){
			cout << "dealing with version " << i + 1 << ":" << endl;

			auto commit_avg = time_loop(
				3, 1.0, 
				[&]() {
					new_ver[i].clear();
				},
				[&]() {
					new_ver[i] = CPAMZ::CPAMZ_commit(all_versions[i], P_insert[i], P_delete[i]);
				},
				[&](){});

			all_versions.emplace_back(new_ver[i]);
			cur_mem = 0;
			mmp.clear();
			for (size_t j = 0; j < all_versions.size(); j++){
				cur_mem += 1.0 * all_versions[j].size_in_bytes(f_noop, mmp); 	// accumulate all version memories, shared pointers only count once
			}
			cout << "[new ver commit time]: " << fixed << setprecision(6) << commit_avg << " Seconds" << endl;
			cout << "[cpamz memory usage]: " << (cur_mem - prev_mem) / 1024.0 / 1024.0  << " MB" << endl;
			prev_mem = cur_mem;
		}
	}


	template<typename PT, typename RQ>
    void plain_spatial_diff_test_latency(PT &P, RQ &range_queries, parlay::sequence<size_t> &batch_sizes, size_t &insert_ratio){
        /*  build tree */
		auto cpamz0 = CPAMZ::CPAMZ_init(P);	//	initial version
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

			auto cpamz1 = CPAMZ::CPAMZ_init(P_newver);
        
			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());
	
			auto avg_time = time_loop(
				3, 1.0,
				[&](){},
				[&](){
					for (size_t i = 0; i < range_queries.size(); i++){
						auto [add, remove] = CPAMZ::plain_map_spatial_diff(cpamz0, cpamz1, range_queries[i]);
						addCnt[i] = add.size();
						removeCnt[i] = remove.size();
					}
				},
				[&]{}
			);
			cout << fixed << setprecision(6) << "[cpamz-plain] spatial-diff time (avg): " << avg_time << endl;
			
			#ifdef TEST
				string file_name = "output/cpamz_spatial_diff_plain-" + to_string(batch_size); 
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
		auto cpamz0 = CPAMZ::CPAMZ_init(P);	//	initial version
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

			auto cpamz1 = CPAMZ::CPAMZ_delete(P_delete, cpamz0);	//	new	version
			auto cpamz2 = CPAMZ::CPAMZ_insert(P_insert, cpamz1); 
			
			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());

			auto avg_time = time_loop(
				3, 1.0,
				[&](){},
				[&](){
					for (size_t i = 0; i < range_queries.size(); i++){
						auto [add, remove] = CPAMZ::map_spatial_diff(cpamz0, cpamz2, range_queries[i]);
						addCnt[i] = add.size();
						removeCnt[i] = remove.size();
					}
				},
				[&]{}
			);
			cout << fixed << setprecision(6) << "[cpamz] spatial-diff time (avg): " << avg_time << endl;
			
			#ifdef TEST
				string file_name = "output/cpamz_spatial_diff-" + to_string(batch_size); 
				ofstream spatialDiffOut(file_name);
				for (size_t i = 0; i < range_queries.size(); i++){
					spatialDiffOut<< addCnt[i] << " " << removeCnt[i] << endl;
				}
			#endif
		}
    }

	/*	different ratio: 1-9, 2-8, 3-7, ..., 8-2, 9-1 */
	template<typename PT, typename RQ>
	auto spatial_diff_test_fix_size(PT P,  RQ &range_queries, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto cpamz0 = CPAMZ::CPAMZ_init(P);	//	initial version
	
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

			auto cpamz1 = CPAMZ::CPAMZ_delete(P_delete, cpamz0);	//	new	version
			auto cpamz2 = CPAMZ::CPAMZ_insert(P_insert, cpamz1); 

			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());

				auto diff_avg = time_loop(
				3, 1.0, [&]() {},
				[&]() {
					parlay::parallel_for(0, range_queries.size(), [&](int i){
						auto [add, remove] = CPAMZ::map_spatial_diff(cpamz0, cpamz2, range_queries[i]);
						addCnt[i] = add.size();
						removeCnt[i] = remove.size();
					});
				},
			[&](){} );
			cout << fixed << setprecision(6) << "[CPAMZ]: spatial-diff time (avg): " << diff_avg << endl;
		}
	}

	/*	50% insertion, 50% deletion	*/
	template<typename PT, typename RQ>
	auto spatial_diff_test_fix_ratio(PT P,  RQ &range_queries, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto cpamz0 = CPAMZ::CPAMZ_init(P);	//	initial version
	
		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) break;
			cout << "[batch-size]: " << batch_size << endl;
			
			auto P_insert = P.substr(0, batch_size / 2);
			auto P_delete = P.substr(P.size() - batch_size / 2, batch_size / 2);
			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += P.size();
			});

			cout << "# of insertion/deletion: " << P_insert.size() << ", " << P_delete.size();

			auto cpamz1 = CPAMZ::CPAMZ_delete(P_delete, cpamz0);	//	new	version
			auto cpamz2 = CPAMZ::CPAMZ_insert(P_insert, cpamz1); 

			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());

		 	auto diff_avg = time_loop(
		    	3, 1.0, [&]() {},
		    	[&]() {
					parlay::parallel_for(0, range_queries.size(), [&](int i){
						auto [add, remove] = CPAMZ::map_spatial_diff(cpamz0, cpamz2, range_queries[i]);
						addCnt[i] = add.size();
						removeCnt[i] = remove.size();
					});
		    	},
	    	[&](){} );
			cout << fixed << setprecision(6) << "[CPAMZ]: spatial-diff time (avg): " << diff_avg << endl;
		}
	}
	


	template<typename PT, typename RQ>
	auto spatial_diff_test(PT P, RQ &range_queries, parlay::sequence<size_t> &batch_sizes, bool &early_end, bool use_hilbert = false){
		auto cpamz0 = CPAMZ::CPAMZ_init(P);	//	initial version
	
		// parlay::sequence<Point> ret;

		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) batch_size = P.size();
			cout << "[batch-size]: " << batch_size << endl;
			
			auto P_insert = P.substr(0, batch_size);
			auto P_delete = P.substr(P.size() - batch_size, batch_size);
			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += P.size();
			});

			auto cpamz1 = CPAMZ::CPAMZ_delete(P_delete, cpamz0);	//	new	version
			auto cpamz2 = CPAMZ::CPAMZ_insert(P_insert, cpamz0); 

			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());



			if (!early_end){
				decltype(cpamz0) commit_ver;
				auto commit_avg = time_loop(
					3, 1.0, [&]() {
						commit_ver.clear();
					},
					[&]() {
						commit_ver = CPAMZ::CPAMZ_commit(cpamz0, P_insert, P_delete);
					},
					[&](){} 
				);

				decltype(cpamz0) merge_ver;
				parlay::sequence<Point> conflict_insert, conflict_update, conflict_delete;
				// cout << cpamz0.size() << ", " << cpamz1.size() << ", " << cpamz2.size() << endl;
				auto merge_avg = time_loop(
					3, 1.0, [&]() {
						merge_ver.clear();
						conflict_insert.clear();
						conflict_update.clear();
						conflict_delete.clear();
					},
					[&]() {
						tie(merge_ver, conflict_insert, conflict_update, conflict_delete)  = CPAMZ::CPAMZ_merge(cpamz0, cpamz1, cpamz2);
					},
					[&](){} 
				);
				cout << "[INFO] commit, merge size: " << commit_ver.size() << ", " << merge_ver.size() << endl; 
				cout << fixed << setprecision(6) << "[CPAMZ]: spatial-commit time (avg): " << commit_avg << endl;
				cout << fixed << setprecision(6) << "[CPAMZ]: spatial-merge time (avg): " << merge_avg << endl;
			}
			
			// if (use_hilbert) cout << "[Hilbert-CPAMZ]: ";
			// else cout << "[Zorder-CPAMZ]: ";

		}
		// return ret;
	}


	template<typename PT>
	void diff_test(PT P, int batch_percent = 10, bool use_hilbert = false){
		auto tree = CPAMZ::CPAMZ_init(P);	//	initial version

		auto batch_size = P.size() * batch_percent / 100;	//	insertion 10%
		auto P_insert = P.substr(0, batch_size);
		parlay::parallel_for(0, P_insert.size(), [&](size_t j){
			P_insert[j].id += P.size();
		});

		auto P_delete = P.substr(0, 2 * batch_size);

		auto new_ver = CPAMZ::CPAMZ_insert(P_insert, tree);	//	new	version
		new_ver = CPAMZ::CPAMZ_delete(P_delete, new_ver); 

		auto add_sz = 0, remove_sz = 0;
	    auto cpam_diff_avg = time_loop(
		    3, 1.0, [&]() {},
		    [&]() {
				auto [add, remove] = CPAMZ::map_diff(tree, new_ver);
				add_sz = add.size();
				remove_sz = remove.size();
		    },
	    [&](){} );
		

		cout << "add size: " << add_sz << endl;
		cout << "remove size: " << remove_sz << endl;

		if (use_hilbert) cout << "[Hilbert-CPAMZ]: ";
		else cout << "[Zorder-CPAMZ]: ";
		cout << fixed << setprecision(6) << "spatial diff time (avg): " << cpam_diff_avg << endl;
	}

	template<typename PT>
	void build_test(PT P, bool use_hilbert = false){
		CPAMZ::zmap tree;
		auto cpam_build_avg = time_loop(
			3, 1.0, [&](){
				tree.clear();
			},
			[&](){
				tree = CPAMZ::CPAMZ_init(P);
			},
		[&](){} );

		// auto [mem_inte_nodes, mem_leaf_nodes] = CPAMZ::size_in_bytes();
		auto [num_inte_nodes, num_leaf_nodes, leaf_size] = tree.node_stats();

		auto f_noop = [&](const auto &et){
			return 0;
		};

		cout << "[cpamz memory usage]: " << endl <<
			"[# of inte nodes]: " << num_inte_nodes << endl << 
			"[# of leaf nodes]: " << num_leaf_nodes << endl <<
			"[tree size]: " << 1.0 * tree.size_in_bytes(f_noop) / 1024.0 / 1024.0 << " MB" << endl;
			// "[memory usage for inte nodes]: " << 1.0 * mem_inte_nodes / 1024.0 / 1024.0 << " MB" << endl <<
			// "[memory usage for leaf nodes]: " << 1.0 * mem_leaf_nodes / 1024.0 / 1024.0 << " MB"  << endl;

		if (use_hilbert) cout << "[Hilbert-CPAMZ]: ";
		else cout << "[Zorder-CPAMZ]: ";
		cout << fixed << setprecision(6) << "build time (avg): " << cpam_build_avg << endl;
	}
	
	template<typename PT>
	void batch_insert_test(PT P, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto n = P.size();
		auto m1 = CPAMZ::CPAMZ_init(P, use_hilbert);	//	build original tree
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
					m2 = CPAMZ::CPAMZ_insert(P2, m1, use_hilbert);
		    	},
	    	[&](){
				if (print_flag){
					cout << "# of points: " << m2.size() << endl;
					print_flag = false;
				}
			} );

			// auto [mem_inte_nodes, mem_leaf_nodes] = CPAMZ::size_in_bytes();
			// auto [num_inte_nodes, num_leaf_nodes, leaf_size] = tree.node_stats();


			cout << "[batch_size]: " << num_processed << endl;
			if (use_hilbert) cout << "[Hilbert-CPAMZ]: ";
			else cout << "[Zorder-CPAMZ]: ";
			cout << fixed << setprecision(6) << "batch insert time (avg): " << cpam_insert_avg << endl;
		}
	}

	template<typename PT>
	void batch_delete_test(PT P, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto m1 = CPAMZ::CPAMZ_init(P, use_hilbert);	//	build original tree
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
					m2 = CPAMZ::CPAMZ_delete(P2, m1, use_hilbert);
		    	},
	    	[&](){
				if (print_flag){
					cout << "# of points: " << m2.size() << endl;
					print_flag = false;
				}			
			} );
			
			cout << "[batch_size]: " << num_processed << endl;
			if (use_hilbert) cout << "[Hilbert-CPAMZ]: ";
			else cout << "[Zorder-CPAMZ]: ";
			cout << fixed << setprecision(6) << "batch delete time (avg): " << cpam_insert_avg << endl;
		}
	}

	template<class PT, class RQ>
    void range_report_test(PT &P, RQ querys, parlay::sequence<size_t> &cnt, bool use_hilbert = false){
	    auto CPAMZ = CPAMZ::CPAMZ_init(P, use_hilbert);
		// for (auto pt: P){
		// 	cout << pt.morton_id << endl;
		// }

		parlay::sequence<size_t> rangeCnt(querys.size());

		for (size_t i = 0; i < querys.size(); i++){
			auto avg_time = time_loop(
				3, 1.0, 
				[&]() {},
				[&]() {					
					rangeCnt[i] = CPAMZ::range_report(CPAMZ, querys[i], use_hilbert).size();
				},
				[&](){} );
			if (rangeCnt[i] != cnt[i]){
				cout << "[ERROR] Incorrect" << endl;
			}
			else{
				cout << fixed << setprecision(6) << rangeCnt[i] << " " << avg_time << endl;
			}
		}



		// if (use_hilbert) cout << "[Hilbert-CPAMZ]: ";
		// else cout << "[Zorder-CPAMZ]: ";
		// cout << fixed << setprecision(6) << "range report time (avg): " << rangeReport_avg << endl;


    }

	template<class PT, class RQ>
    void range_count_test(PT &P, RQ querys, parlay::sequence<size_t> &cnt, bool use_hilbert = false){
	    range_report_test(P, querys, cnt);
    }

}

