#include <atomic>
#pragma once
// --- From mvzd/src/mvq/node.hpp ---
#include "global_config.hpp"
#include "geo/point.hpp"
#include "geo/operations.hpp"
#include "geo/io.hpp"
#include <parlay/internal/get_time.h>
#include "pam/utils.h"

#define USE_PT


namespace mvq{
    inline std::atomic<size_t> global_live_mem(0);
    inline std::atomic<size_t> query_nodes_touched(0);


	using namespace std;
	using namespace geobase;
	using parlay::sequence;
	using parlay::par_do;
	using parlay::par_do_if;

	struct tree_stat{
		size_t num_inte_nodes;
		size_t num_leaf_nodes;
		size_t mem_inte_nodes;
		size_t mem_leaf_nodes;
		tree_stat():num_inte_nodes(0), num_leaf_nodes(0), mem_inte_nodes(0), mem_leaf_nodes(0){}
		auto get_total(){
			return 1.0 * (mem_inte_nodes + mem_leaf_nodes) / 1024.0 / 1024.0;
		}
	};

	struct BaseNode{
		#ifdef USE_MBR
		Bounding_Box mbr;
		BaseNode(): mbr({Point(-1, -1), Point(-1, -1)}){}
		#endif
		BaseNode(){}

		virtual~BaseNode() = default;
		virtual bool is_leaf(){ return false; }
		virtual size_t get_num_points(){ return 0; }
	};

	struct InteNode: BaseNode{
		shared_ptr<BaseNode> l_son, r_son;
		size_t num_pts;

		InteNode(): l_son(nullptr), r_son(nullptr), num_pts(0){ mvq::global_live_mem.fetch_add(sizeof(InteNode), std::memory_order_relaxed); }
		virtual ~InteNode() { mvq::global_live_mem.fetch_sub(sizeof(InteNode), std::memory_order_relaxed); }
		// InteNode(InteNode &x): l_son(x->l_son), r_son(x->r_son), num_pts(x->num_pts){}
		
		virtual bool is_leaf(){ return false; }
		virtual size_t get_num_points(){ return num_pts; }
	};

	struct LeafNode: BaseNode{
		// sequence<Point> records;
		sequence<Point> records = sequence<Point>::uninitialized(32);
		size_t tracked_mem = 0;
		void update_mem() {
			size_t current = sizeof(LeafNode) + records.capacity() * sizeof(Point);
			if (current > tracked_mem) {
				mvq::global_live_mem.fetch_add(current - tracked_mem, std::memory_order_relaxed);
			} else if (current < tracked_mem) {
				mvq::global_live_mem.fetch_sub(tracked_mem - current, std::memory_order_relaxed);
			}
			tracked_mem = current;
		}
		virtual ~LeafNode() {
			mvq::global_live_mem.fetch_sub(tracked_mem, std::memory_order_relaxed);
		}

		template<typename Records>
		LeafNode(Records &r){
			if (r.size() > 32){
				records = sequence<Point>::uninitialized(r.size());
			}
			size_t i = 0;	
			for (auto &pt: r){
				parlay::assign_uninitialized(records[i++], pt);
				// records[i] = r[i];
			}
			records.resize(r.size());
			update_mem();
		}

		template<typename Records, typename Func>
		LeafNode(Records &r, Func &f){
			if (r.size() > 32){
				records = sequence<Point>::uninitialized(r.size());
			}
			size_t i = 0;	
			for (auto &pt: r){
				if (f){
					parlay::assign_uninitialized(records[i++], pt);
				}
				// records[i] = r[i];
			}
			records.resize(i);
			update_mem();
		}

		// LeafNode(LeafNode &x): records(x->records){}

		virtual bool is_leaf(){ return true; }
		virtual size_t get_num_points(){ return records.size(); }

		void print_records(){
			cout << records.size() << endl;
			for (size_t i = 0; i < records.size(); i++){
				cout << "(" << records[i].x << ", " << records[i].y << ")" << endl; 
			}
		}
	};



// --- From mvzd/src/mvq/tree.hpp ---
class Tree{
	public:
		size_t visited_leaf = 0;
		size_t visited_inte = 0;
		size_t granularity_cutoff = 1000;
		size_t inte_node_cnt = 0;
		size_t leaf_node_cnt = 0;
		size_t leaf_size = 1;
		size_t node_cnt = 0;
		size_t leaf_cnt = 0;
		size_t record_cnt = 0;
		size_t case0 = 0;
		size_t case1 = 0;
		size_t case2 = 0;
		size_t case3 = 0;
		size_t case4 = 0;
		FT case0_time = 0;
		FT case1_time = 0;
		FT case2_time = 0;
		FT case3_time = 0;
		FT case4_time = 0;
		shared_ptr<BaseNode> root;
		vector<shared_ptr<BaseNode> > multi_version_roots = {};

		Tree(size_t _leaf_sz);

		shared_ptr<BaseNode> node_copy(shared_ptr<BaseNode> &x);

		template<typename Func>
		shared_ptr<BaseNode> node_copy(shared_ptr<BaseNode> &x, Func &f);

		// entrance of building zdtree & tree construction
		shared_ptr<BaseNode> build(std::span<Point> P, size_t b);
		void build(std::span<Point> P);


		void gc_nocheck(shared_ptr<BaseNode> &x);
		void garbage_collect(shared_ptr<BaseNode> &base, shared_ptr<BaseNode> &x);
        void garbage_collect(shared_ptr<BaseNode> &base, shared_ptr<BaseNode> &new_ver, shared_ptr<BaseNode> &x);

		void clear();

		int tree_hash(shared_ptr<BaseNode> &x);
        int tree_hash_non_associative(shared_ptr<BaseNode> &x);

		size_t num_of_nodes();	// total number of nodes in mvq index
		auto leaf_size_status();
		auto leaf_size_status(shared_ptr<BaseNode> &rt);

		auto get_tree_statistics();
		template<typename NodeMap>
		void get_tree_statistics(shared_ptr<BaseNode> &x, tree_stat &stat, NodeMap &mmp);


		auto print_leaf(shared_ptr<BaseNode> &x);

		void merge_nodes(shared_ptr<BaseNode> &lhs, shared_ptr<BaseNode> &rhs, shared_ptr<InteNode> &cur); 	// merge two sons to current node.
		void delete_merge_nodes(shared_ptr<BaseNode> &L, shared_ptr<BaseNode> &R, InteNode* cur_node);

		shared_ptr<InteNode> create_internal(shared_ptr<BaseNode> &L, shared_ptr<BaseNode> &R);	//	create an internal node, do not store pointers to original records.
		shared_ptr<LeafNode> create_leaf(std::span<Point> P, size_t b); // create a leaf, store all (pointers of) records.

		// 	in-place insertion
		void batch_insert_sorted(std::span<Point> P);
		void batch_insert_sorted_node(shared_ptr<BaseNode> &x, std::span<Point> P, size_t b);
		// 	multi-version (persistent) insertion
		shared_ptr<BaseNode> multi_version_batch_insert_sorted(std::span<Point> P, shared_ptr<BaseNode> &old_root);
		void multi_version_batch_insert_sorted_node(shared_ptr<BaseNode> &x, std::span<Point> P, size_t b);

		//	in-place deletion
		void batch_delete_sorted(std::span<Point> P);
		void batch_delete_sorted_node(shared_ptr<BaseNode> &x, std::span<Point> P, size_t b);
		//	multi-version (persistent) deletion
		shared_ptr<BaseNode> multi_version_batch_delete_sorted(std::span<Point> P, shared_ptr<BaseNode> &old_root);
		void multi_version_batch_delete_sorted_node(shared_ptr<BaseNode> &x, std::span<Point> P, size_t b);



		// collect all records in version x	
		template<typename Out>
		auto collect_records(shared_ptr<BaseNode> &x, Out &out, size_t &cnt);
		auto collect_records(shared_ptr<BaseNode> &x);
		auto collect_records(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, 
							Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter);

		// auto diff(shared_ptr<BaseNode> &tree1, shared_ptr<BaseNode> &tree2);
		template<typename DIFF>
		auto spatial_diff(shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, 
							FT x_prefix, FT y_prefix, size_t b, bool x_splitter, DIFF &ret_diff);

		template<typename JOIN_RES>
		void spatial_join(shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2, FT &point_dis, JOIN_RES &join_res,
								Bounding_Box &mbr1, FT x_prefix1, FT y_prefix1, size_t b1, bool x_splitter1,
								Bounding_Box &mbr2, FT x_prefix2, FT y_prefix2, size_t b2, bool x_splitter2);
		
		template<typename JOIN_RES>
		void two_version_spatial_join(shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2, FT &point_dis, JOIN_RES &join_res,
									Bounding_Box mbr1, Bounding_Box mbr2);

		template<typename DIFF>
		auto spatial_two_version_diff(shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, 
										DIFF &ret_diff);

		auto diff(shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2, size_t b);
		auto leaf_leaf_diff(shared_ptr<BaseNode> &leaf1, shared_ptr<BaseNode> &leaf2);

		template<typename DIFF>
		auto leaf_leaf_diff(shared_ptr<BaseNode> &leaf1, shared_ptr<BaseNode> &leaf2, Bounding_Box &query_mbr, DIFF &ret_diff);

		auto leaf_inte_diff(std::span<Point> P, shared_ptr<BaseNode> &inte, size_t b);

		template<typename DIFF>
		auto leaf_inte_diff(std::span<Point> P, shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr,
							 FT x_prefix, FT y_prefix, size_t b, bool x_splitter, DIFF &ret_diff, bool reverse = false);

		auto two_version_diff(shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2);

		//	commit, based on op do multi-version insertion, deletion, and update
		auto commit(shared_ptr<BaseNode> &old_version, sequence<Point> &P_insert, sequence<Point> &P_delete);
		// merge, works like git merge, return a conflict set
		auto merge(shared_ptr<BaseNode> &base, shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2);
		
		// intervals are left close right open: i.e., [L, R)
		// range report
		template <class Out> 
		void range_report_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter, size_t &cnt, Out &out);
		template <class Out>
		void range_report(Bounding_Box &query_mbr, Bounding_Box &cur_mbr, size_t &cnt, Out &out);
		template <class Out>
		void range_report(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, size_t &cnt, Out &out);

		// range count
		size_t range_count_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter);
		size_t range_count(Bounding_Box &query_mbr, Bounding_Box &cur_mbr);
		size_t range_count(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr);

		// k nearest neighbor report
		template<class T> 
		void knn_report_node(shared_ptr<BaseNode> &x, size_t &k, Point query_point, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter, T &nn_res);
		auto knn_report(size_t &k, Point query_point, Bounding_Box &cur_mbr);
    };

}

// --- From mvzd/src/mvq/tree_impl.hpp ---
namespace mvq {


	// 	sequence<Point> insert_points = {}, delete_points = {}, update_points = {};
	// 	size_t i = 0, j = 0;
	// 	while (i < sorted_add.size() && j < sorted_remove.size()){
	// 		if (sorted_add[i].id < sorted_remove[j].id){	//	we should insert sorted_add[i]
	// 			insert_points.emplace_back(sorted_add[i++]);
	// 		}
	// 		else if (sorted_add[i].id == sorted_remove[j].id){	//	exist in both add and remove, indicate it is an update point
	// 			update_points.emplace_back(sorted_add[i++]);
	// 			j++;
	// 		}
	// 		else{	//	we should remove sorted_remove[j]
	// 			delete_points.emplace_back(sorted_remove[j++]);
	// 		}
	// 	}
	// 	while (i < sorted_add.size()) insert_points.emplace_back(sorted_add[i++]);
	// 	while (j < sorted_remove.size()) delete_points.emplace_back(sorted_remove[j++]);
	// 	return make_tuple(insert_points, delete_points, update_points);
	// }
	auto Tree::collect_records(shared_ptr<BaseNode> &x){
		sequence<Point> ret = {};
		if (!x){
			return ret;
		}
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			ret = cur_leaf->records;
			return ret;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		sequence<Point> R;
		auto collect_left = [&](){ret = collect_records(cur_inte->l_son); }; 
		auto collect_right = [&](){R = collect_records(cur_inte->r_son); };
		par_do_if(x->get_num_points() >= granularity_cutoff,
			collect_left,
			collect_right
		);
		ret.append(R);
		return ret;
	}
	
	template<typename Out>
	auto Tree::collect_records(shared_ptr<BaseNode> &x, Out &out, size_t &cnt){
		if (!x){
			return;
		}
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (auto &p: cur_leaf->records){
				out[cnt++] = p;
			}
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		size_t R_cnt = cnt;
		if (cur_inte->l_son) R_cnt += cur_inte->l_son->get_num_points();

		auto collect_left = [&](){ collect_records(cur_inte->l_son, out, cnt); }; 
		auto collect_right = [&](){ collect_records(cur_inte->r_son, out, R_cnt); };
		par_do_if(x->get_num_points() >= granularity_cutoff,
			collect_left,
			collect_right
		);
		return;
	}

	auto Tree::collect_records(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, 
								Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter){
		sequence<Point> ret = {};
		if (!x){
			return ret;
		}
		auto flag = mbr_mbr_relation(cur_mbr, query_mbr);
		if (flag < 0) return ret;

		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			if (flag > 0) return cur_leaf->records;
			else{
				auto f = [&](auto x){	//	filter a sequence of points in the query mbr
					return point_in_mbr(x, query_mbr);
				};
				return parlay::filter(cur_leaf->records, f);
			}
		}

		auto cur_inte = static_cast<InteNode*>(x.get());

		size_t shift_b = (b + 1) / 2;
		FT split_value = 1.0 * (1u << (shift_b - 1));

		auto L_box = cur_mbr;
		auto R_box = cur_mbr;
		auto rx_prefix = x_prefix;
		auto ry_prefix = y_prefix;
		if (x_splitter){
			L_box.second.x = min(x_prefix + split_value - FT_EPS, L_box.second.x);
			R_box.first.x = max(x_prefix + split_value, R_box.first.x);
			rx_prefix += split_value;
		}
		else{
			L_box.second.y = min(y_prefix + split_value - FT_EPS, L_box.second.y);
			R_box.first.y = max(y_prefix + split_value, R_box.first.y);
			ry_prefix += split_value;
		}

		sequence<Point> R;
		auto collect_left = [&](){ret = collect_records(cur_inte->l_son, query_mbr, L_box, x_prefix, y_prefix, b - 1, !x_splitter); }; 
		auto collect_right = [&](){R = collect_records(cur_inte->r_son, query_mbr, R_box, rx_prefix, ry_prefix, b - 1, !x_splitter); };
		par_do_if(x->get_num_points() >= granularity_cutoff,
			collect_left,
			collect_right
		);
		ret.append(R);
		return ret;
	}


	Tree::Tree(size_t _leaf_sz){
		leaf_size = _leaf_sz;
	}

	void Tree::clear(){
		root.reset();
		multi_version_roots.clear();
	}

	void Tree::delete_merge_nodes(shared_ptr<BaseNode> &L, shared_ptr<BaseNode> &R, InteNode* cur_node){
		// deal with MBR, covered points of parent
		auto L_num_pts = (L == nullptr) ? 0 : L->get_num_points();
		auto R_num_pts = (R == nullptr) ? 0 : R->get_num_points();

		cur_node->num_pts = L_num_pts + R_num_pts;
		cur_node->l_son = move(L);
		cur_node->r_son = move(R);
	}

	shared_ptr<InteNode> Tree::create_internal(shared_ptr<BaseNode> &L, shared_ptr<BaseNode> &R){
  		// parlay::internal::timer t;
		// t.start();
		shared_ptr<InteNode> cur_node(new InteNode());
		
		// augmented changes happen here
		merge_nodes(L, R, cur_node);
		// t.next_time();
		// zd_build_break_down.inte_time += t.total_time();
		return cur_node;		
	}

	shared_ptr<LeafNode> Tree::create_leaf(std::span<Point> P, size_t b){
		shared_ptr<LeafNode> cur_node(new LeafNode(P));
		return cur_node;
	}

	shared_ptr<BaseNode> Tree::build(std::span<Point> P, size_t b){
		if (!b || (P.size() <= leaf_size)){
			return create_leaf(P, b);
			// return nullptr;	//	measure leaf time
		}

		// parlay::internal::timer t;
		// t.start();
		auto splitter = split_by_bit(P, 0, P.size(), b);
		// t.next_time();
		// zd_build_break_down.split_time += t.total_time();
		shared_ptr<BaseNode> L = nullptr;
		shared_ptr<BaseNode> R = nullptr;
		auto build_left = [&](){ if (0 < splitter) L = build(P.subspan(0, splitter), b - 1); }; 
		auto build_right = [&](){ if (splitter < P.size()) R = build(P.subspan(splitter), b - 1); }; 
		par_do_if(P.size() >= granularity_cutoff,
			build_left, 
			build_right);
		return create_internal(L, R);
	}

	void Tree::build(std::span<Point> P){
		if (P.empty()) return;
		root = build(P, 64);
		// multi_version_roots.emplace_back(root);
	}

	void Tree::gc_nocheck(shared_ptr<BaseNode> &x){
		if (!x) return;
		if (x->is_leaf()){
			x.reset();
			return;
		}

		auto cur_inte = static_cast<InteNode*>(x.get());
		if (x.use_count() == 1){
			auto gc_left = [&](){ gc_nocheck(cur_inte->l_son); }; 
			auto gc_right = [&](){gc_nocheck(cur_inte->r_son); };
			par_do_if(x->get_num_points() > granularity_cutoff, gc_left, gc_right);
		}
		x.reset();
	}

	void Tree::garbage_collect(shared_ptr<BaseNode> &base, shared_ptr<BaseNode> &x){
		if (!x || x == base) return;
		if (x->is_leaf()){
			x.reset();
			return;
		}
		
		if (!base || base->is_leaf()){	//	base is null or leaf, just garbage_collect all nodes in x.
			gc_nocheck(x);
			return;
		}

		auto base_inte =static_cast<InteNode*>(base.get()); 
		auto cur_inte = static_cast<InteNode*>(x.get());

		if (x->get_num_points() > granularity_cutoff){
			auto gc_left = [&](){ garbage_collect(base_inte->l_son, cur_inte->l_son); };
			auto gc_right = [&](){	garbage_collect(base_inte->r_son, cur_inte->r_son); };
			par_do(gc_left, gc_right); 
		}
		x.reset();
	}


	void Tree::garbage_collect(shared_ptr<BaseNode> &base, shared_ptr<BaseNode> &new_ver, shared_ptr<BaseNode> &x){
		if (!x || x == base || x == new_ver) {
			cout << "[DEBUG] reused node found" << endl;
			return;
		}

		if (x->is_leaf()){
			cout << "[DEBUG] freed node found" << endl;
			x.reset();
			return;
		}
		
		if (!base || base->is_leaf()){	//	base is null or leaf, garbage collect based on x and new_ver.
			cout << "[DEBUG] go to new_ver-x" << endl;
			garbage_collect(new_ver, x);
			return;
		}
		else if (!new_ver || new_ver->is_leaf()){	//	new_ver is null or leaf, garbage collect based on x and base.
			cout << "[DEBUG] go to base-x" << endl;
			garbage_collect(base, x);
			return;
		}
		
		auto base_inte =static_cast<InteNode*>(base.get()); 
		auto cur_inte = static_cast<InteNode*>(x.get());
		auto newver_inte = static_cast<InteNode*>(new_ver.get());

		if (x->get_num_points() > granularity_cutoff){
			auto gc_left = [&](){ garbage_collect(base_inte->l_son, cur_inte->l_son, newver_inte->l_son); };
			auto gc_right = [&](){	garbage_collect(base_inte->r_son, cur_inte->r_son, newver_inte->r_son); };
			par_do(gc_left, gc_right); 
		}
		x.reset();
	}

	void Tree::batch_insert_sorted_node(shared_ptr<BaseNode> &x, std::span<Point> P, size_t b){
		if (x == nullptr){
			x = build(P, b);
			return;
		}
		auto less = [&](auto lhs, auto rhs){
			return lhs.morton_id < rhs.morton_id || 
					(lhs.morton_id == rhs.morton_id && lhs.id < rhs.id);
			// auto msd = 0;
			// if (geobase::less_msb(static_cast<unsigned int>(lhs.x) ^ static_cast<unsigned int>(rhs.x), static_cast<unsigned int>(lhs.y) ^ static_cast<unsigned int>(rhs.y))) 
			// 	msd = 1;
			// return !msd ? lhs.x < rhs.x : lhs.y < rhs.y;
		};
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			if (!b || cur_leaf->records.size() + P.size() <= leaf_size){	// current leaf is not full
				cur_leaf->records = parlay::merge(cur_leaf->records, parlay::make_slice(P), less);
				cur_leaf->update_mem();
				return;
			}
			else{
				auto new_points = parlay::merge(cur_leaf->records, parlay::make_slice(P), less);
				x = build(std::span<Point>(new_points.begin(), new_points.end()), b);
				return;
			}
		}
		auto splitter = split_by_bit(P, 0, P.size(), b);
		auto cur_inte = static_cast<InteNode*>(x.get());
		auto insert_left = [&](){ if (0 < splitter){ batch_insert_sorted_node(cur_inte->l_son, P.subspan(0, splitter), b - 1); }; };
		auto insert_right = [&](){ if (splitter < P.size()){ batch_insert_sorted_node(cur_inte->r_son, P.subspan(splitter), b - 1); }; };
		par_do_if(P.size() >= granularity_cutoff,
			insert_left,
			insert_right); 
		delete_merge_nodes(cur_inte->l_son, cur_inte->r_son, cur_inte);
	}

	void Tree::batch_insert_sorted(std::span<Point> P){
		if (P.empty()) return;
		if (root == nullptr) build(P);
		else batch_insert_sorted_node(root, P, 64);
	}

	shared_ptr<BaseNode> Tree::node_copy(shared_ptr<BaseNode> &x){
		if (!x) return nullptr;
		if (x->is_leaf()){
			// cout << "copying leaf node" << endl;
			// parlay::internal::timer t("leaf time");

			auto cur_leaf = static_cast<LeafNode*>(x.get());
			shared_ptr<LeafNode> new_leaf(new LeafNode(cur_leaf->records));

			// utils::fetch_and_add(&zd_leaf_copy_time, t.next_time());
			// zd_leaf_copy_time += t.next_time();
			return new_leaf;
		}
		else{
			// cout << "copying inte node" << endl;
			// parlay::internal::timer t("inte time");

			auto cur_inte = static_cast<InteNode*>(x.get());
			shared_ptr<InteNode> new_inte(new InteNode());
			new_inte->num_pts = cur_inte->num_pts;
			new_inte->l_son = cur_inte->l_son;
			new_inte->r_son = cur_inte->r_son;

			// utils::fetch_and_add(&zd_inte_copy_time, t.next_time());
			// zd_inte_copy_time += t.next_time();
			return new_inte;
		}
	}

	template<typename Func>
	shared_ptr<BaseNode> Tree::node_copy(shared_ptr<BaseNode> &x, Func &f){
		if (!x) return nullptr;
		if (x->is_leaf()){
			// cout << "copying leaf node" << endl;
			auto cur_leaf = static_cast<LeafNode*>(x.get());

			shared_ptr<LeafNode> new_leaf(new LeafNode(cur_leaf->records, f));
			return new_leaf;
		}
		else{
			// cout << "copying inte node" << endl;
			auto cur_inte = static_cast<InteNode*>(x.get());
			shared_ptr<InteNode> new_inte(new InteNode());
			new_inte->num_pts = cur_inte->num_pts;
			new_inte->l_son = cur_inte->l_son;
			new_inte->r_son = cur_inte->r_son;
			return new_inte;
		}
	}

	void Tree::batch_delete_sorted_node(shared_ptr<BaseNode> &x, std::span<Point> P, size_t b){
		if (x == nullptr){
			return;
		}
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			auto tmp_records = parlay::make_slice(P);
			auto [added, removed] = merge_pts(cur_leaf->records, tmp_records);
			cur_leaf->records = removed;
			// cur_leaf->records = parlay::filter(cur_leaf->records, not_in);
			return;
		}
		auto splitter = split_by_bit(P, 0, P.size(), b);
		auto cur_inte = static_cast<InteNode*>(x.get());
		auto delete_left = [&](){ if (0 < splitter){ batch_delete_sorted_node(cur_inte->l_son, P.subspan(0, splitter), b - 1); }; };
		auto delete_right = [&](){ if (splitter < P.size()){ batch_delete_sorted_node(cur_inte->r_son, P.subspan(splitter), b - 1); }; };
		par_do_if(P.size() >= granularity_cutoff,
			delete_left,
			delete_right); 

		delete_merge_nodes(cur_inte->l_son, cur_inte->r_son, cur_inte);
	}

	void Tree::batch_delete_sorted(std::span<Point> P){
		if (P.empty() || root == nullptr) return;
		else batch_delete_sorted_node(root, P, 64);
	}

	int Tree::tree_hash(shared_ptr<BaseNode> &x){
		node_cnt++;
		if (x->is_leaf()){
			leaf_cnt++;
			int ret = 0;
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			// cout << "leaf size = " << cur_leaf->records.size() << endl;
			record_cnt += cur_leaf->records.size();
			for (auto &pt: cur_leaf->records){
				int x_val = static_cast<int>(pt.x);
				int y_val = static_cast<int>(pt.y);
				ret ^= parlay::hash64(x_val ^ y_val);
			}
			return parlay::hash64(ret);
		}
		int L = 0, R = 0;
		auto cur_inte = static_cast<InteNode*>(x.get());
		if (cur_inte->l_son != nullptr) L = tree_hash(cur_inte->l_son);
		if (cur_inte->r_son != nullptr) R = tree_hash(cur_inte->r_son);
		return parlay::hash64(L ^ R);
	}

	const unsigned long long MOD = 1e9 + 7;
	int Tree::tree_hash_non_associative(shared_ptr<BaseNode> &x){
		node_cnt++;
		if (x->is_leaf()){
			leaf_cnt++;
			int ret = 0;
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			// cout << "leaf size = " << cur_leaf->records.size() << endl;
			record_cnt += cur_leaf->records.size();
			auto cur = 1;
			for (auto &pt: cur_leaf->records){
				int x_val = static_cast<int>(pt.x);
				int y_val = static_cast<int>(pt.y);
				ret = (ret + (cur * x_val % MOD) + (cur * y_val % MOD)) % MOD;
				cur++;
			}
			// cout << ret << endl;
			return parlay::hash64(ret);
		}
		int L = 0, R = 0;
		auto cur_inte = static_cast<InteNode*>(x.get());
		if (cur_inte->l_son != nullptr) L = tree_hash_non_associative(cur_inte->l_son);
		if (cur_inte->r_son != nullptr) R = tree_hash_non_associative(cur_inte->r_son);
		R = (2 * R) % MOD;
		return parlay::hash64(L ^ R);
	}

	template<typename LEAF_SZ_MAP>
	void collect_leaf_size(shared_ptr<BaseNode> &x, LEAF_SZ_MAP &mmp){
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			mmp[cur_leaf->get_num_points()]++;
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		if (cur_inte->l_son != nullptr) collect_leaf_size(cur_inte->l_son, mmp);
		if (cur_inte->r_son != nullptr) collect_leaf_size(cur_inte->r_son, mmp);
	}

	template<typename NODE_MAP>
	void node_traverse(shared_ptr<BaseNode> &x, NODE_MAP &mmp){
		mmp[x] = 1;	//	node appears
		if (x->is_leaf()){
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		if (cur_inte->l_son != nullptr) node_traverse(cur_inte->l_son, mmp);
		if (cur_inte->r_son != nullptr) node_traverse(cur_inte->r_son, mmp);
	}

	//	traverse sequentially
	template<typename NodeMap>
	void Tree::get_tree_statistics(shared_ptr<BaseNode> &x, tree_stat &stat, NodeMap &mmp){
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			if (!mmp[x]) {
				stat.num_leaf_nodes++;
				stat.mem_leaf_nodes += sizeof(LeafNode) + cur_leaf->records.capacity() * sizeof(Point);
				mmp[x] = 1;
			}
			return;
		}
		if (!mmp[x]){
			stat.num_inte_nodes++;
			stat.mem_inte_nodes += sizeof(InteNode);
			mmp[x] = 1;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		if (cur_inte->l_son) get_tree_statistics(cur_inte->l_son, stat, mmp);
		if (cur_inte->r_son) get_tree_statistics(cur_inte->r_son, stat, mmp);
	}

	auto Tree::get_tree_statistics(){
		tree_stat stat;
		unordered_map<shared_ptr<BaseNode>, bool> node_map;
		for (auto &rt: multi_version_roots){
			if (!rt){
				cout << "Released version" << endl;
				continue;
			}
			get_tree_statistics(rt, stat, node_map);
		}
		return stat;
	}
	
	auto Tree::leaf_size_status(){
		unordered_map<int, int> leaf_map = {};
		for (auto &rt: multi_version_roots){
			if (rt == nullptr){
				cout << "Released version" << endl;
				continue;
			}
			collect_leaf_size(rt, leaf_map);
		}
		return leaf_map;
	}

	auto Tree::leaf_size_status(shared_ptr<BaseNode>& rt){
		unordered_map<int, int> leaf_map = {};
		if (rt == nullptr) return leaf_map;
		collect_leaf_size(rt, leaf_map);
		return leaf_map;
	}

	auto Tree::print_leaf(shared_ptr<BaseNode> &x){
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			cur_leaf->print_records();
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		if (cur_inte->l_son != nullptr) print_leaf(cur_inte->l_son);
		if (cur_inte->r_son != nullptr) print_leaf(cur_inte->r_son);
	}

	

	size_t Tree::num_of_nodes(){
		cout << "[INFO] zdtree # of versions: " << multi_version_roots.size() << endl;
		unordered_map<shared_ptr<BaseNode>, bool> node_map;
		for (auto &rt: multi_version_roots){
			if (rt == nullptr){
				cout << "Released version" << endl;
				continue;
			}
			node_traverse(rt, node_map);
			// cout << "current version tree nodes: " << node_map.size() << endl;
		}
		return node_map.size();
	}

}

// --- From mvzd/src/mvq/query.hpp ---
namespace mvq {
	size_t Tree::range_count_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter){
		// print_mbr(cur_mbr);
		// print_mbr(query_mbr);
		int flag = mbr_mbr_relation(cur_mbr, query_mbr);
		if (flag < 0) return 0;
		if (flag > 0) {
			// cout << x->get_num_points() << endl;
			return x->get_num_points();
		}
		if (x->is_leaf()){	// we have to scan the leaf to report the number of points;
			size_t ret = 0;
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (auto &p: cur_leaf->records){
				if (point_in_mbr(p, query_mbr)){
					ret += 1;
				}
			}
			// cout << ret << endl;
			return ret;
		}
		else{
			auto cur_inte = static_cast<InteNode*>(x.get());
			auto [L_box, R_box, rx_prefix, ry_prefix] = geobase::compute_cur_box(cur_mbr, x_prefix, y_prefix, b, x_splitter);
			size_t ret_L = 0, ret_R = 0;
			if (cur_inte->l_son != nullptr) {
				ret_L = range_count_node(cur_inte->l_son, query_mbr, L_box, x_prefix, y_prefix, b - 1, !x_splitter);
			}
			if (cur_inte->r_son != nullptr) {
				ret_R = range_count_node(cur_inte->r_son, query_mbr, R_box, rx_prefix, ry_prefix, b - 1, !x_splitter);
			}
			return ret_L + ret_R;
		}
		return -1; // unexpected error happens if the code runs to here.
	}

	size_t Tree::range_count(Bounding_Box &query_mbr, Bounding_Box &cur_mbr){
		// size_t ret = range_count_node(root, query_mbr, cur_mbr, true);		
		size_t ret = range_count_node(root, query_mbr, cur_mbr, 0.0, 0.0, 64, true);
		return ret;
	}

	size_t Tree::range_count(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr){
		// size_t ret = range_count_node(root, query_mbr, cur_mbr, true);		
		size_t ret = range_count_node(x, query_mbr, cur_mbr, 0.0, 0.0, 64, true);
		return ret;
	}

	template <class Out> 
	void Tree::range_report_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter, size_t &cnt, Out &out){
		mvq::query_nodes_touched.fetch_add(1, std::memory_order_relaxed);
		if (!x){
			return;
		}
		// print_mbr(cur_mbr);
		// print_mbr(query_mbr);
		auto flag = mbr_mbr_relation(cur_mbr, query_mbr);
		if (flag < 0) return;

		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (auto &p: cur_leaf->records){
				if (point_in_mbr(p, query_mbr)){
					out[cnt++] = p;
				}
			}
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		auto [L_box, R_box, rx_prefix, ry_prefix] = compute_cur_box(cur_mbr, x_prefix, y_prefix, b, x_splitter);
		range_report_node(cur_inte->l_son, query_mbr, L_box, x_prefix, y_prefix, b - 1, !x_splitter, cnt, out);
		range_report_node(cur_inte->r_son, query_mbr, R_box, rx_prefix, ry_prefix, b - 1, !x_splitter, cnt, out);


	}

	template <class Out>
	void Tree::range_report(Bounding_Box &query_mbr, Bounding_Box &cur_mbr, size_t &cnt, Out &out){
		range_report_node(root, query_mbr, cur_mbr, 0.0, 0.0, 64, true, cnt, out);
	}

	template <class Out>
	void Tree::range_report(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, size_t &cnt, Out &out){
		range_report_node(x, query_mbr, cur_mbr, 0.0, 0.0, 64, true, cnt, out);
	}

	#ifdef USE_MBR
	template <class T>
	void Tree::knn_report_node(shared_ptr<BaseNode> &x, size_t &k, Point* query_point, T &nn_res){
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (auto &p: cur_leaf->records){
				auto cur_sqrdis = point_point_sqrdis(p, query_point);
				if (nn_res.size() < k){
					nn_res.push({p, cur_sqrdis});
				}
				else {
					nn_pair current_pair = {p, cur_sqrdis};
					nn_pair_cmp cmp;
					if (cmp(current_pair, nn_res.top())) {
						nn_res.pop();
						nn_res.push(current_pair);
					}
				}
			}
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		auto l_son_sqrdis = FT_INF_MAX, r_son_sqrdis = FT_INF_MAX;
		if (cur_inte->l_son != nullptr) l_son_sqrdis = point_mbr_sqrdis(query_point, cur_inte->l_son->mbr);
		if (cur_inte->r_son != nullptr) r_son_sqrdis = point_mbr_sqrdis(query_point, cur_inte->r_son->mbr);
		if (l_son_sqrdis <= r_son_sqrdis){ // first go left
			if (nn_res.size() < k || l_son_sqrdis <= nn_res.top().second){
				knn_report_node(cur_inte->l_son, k, query_point, nn_res);
			}
			if (nn_res.size() < k || r_son_sqrdis <= nn_res.top().second){
				knn_report_node(cur_inte->r_son, k, query_point, nn_res);
			}
		}
		else{	// first go right
			if (nn_res.size() < k || r_son_sqrdis <= nn_res.top().second){
				knn_report_node(cur_inte->r_son, k, query_point, nn_res);
			}
			if (nn_res.size() < k || l_son_sqrdis <= nn_res.top().second){
				knn_report_node(cur_inte->l_son, k, query_point, nn_res);
			}
		}
		return;
	}

	auto Tree::knn_report(size_t &k, Point* query_point){
		priority_queue<nn_pair, vector<nn_pair>, nn_pair_cmp> nn_res;
		knn_report_node(root, k, query_point, nn_res);
		return nn_res;
	}
	#endif

	auto Tree::knn_report(size_t &k, Point query_point, Bounding_Box &cur_mbr){
		priority_queue<nn_pair, vector<nn_pair>, nn_pair_cmp> nn_res;
		knn_report_node(root, k, query_point, cur_mbr, 0.0, 0.0, 64, true, nn_res);
		return nn_res;
	}

	template <class T>
	void Tree::knn_report_node(shared_ptr<BaseNode> &x, size_t &k, Point query_point, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter, T &nn_res){
		mvq::query_nodes_touched.fetch_add(1, std::memory_order_relaxed);
		if (!x) return;
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (auto &p: cur_leaf->records){
				auto cur_sqrdis = point_point_sqrdis(p, query_point);
				if (nn_res.size() < k){
					nn_res.push({p, cur_sqrdis});
				}
				else {
					nn_pair current_pair = {p, cur_sqrdis};
					nn_pair_cmp cmp;
					if (cmp(current_pair, nn_res.top())) {
						nn_res.pop();
						nn_res.push(current_pair);
					}
				}
			}
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		auto l_son_sqrdis = FT_INF_MAX, r_son_sqrdis = FT_INF_MAX;
		auto [L_box, R_box, rx_prefix, ry_prefix] = geobase::compute_cur_box(cur_mbr, x_prefix, y_prefix, b, x_splitter);
		
		if (cur_inte->l_son != nullptr) {
			l_son_sqrdis = point_mbr_sqrdis(query_point, L_box);
		}
		if (cur_inte->r_son != nullptr) {
			r_son_sqrdis = point_mbr_sqrdis(query_point, R_box);
		}
		
		if (l_son_sqrdis <= r_son_sqrdis){ // first go left
			if (nn_res.size() < k || l_son_sqrdis <= nn_res.top().second){
				knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix, y_prefix, b - 1, !x_splitter, nn_res);
			}
			if (nn_res.size() < k || r_son_sqrdis <= nn_res.top().second){
				knn_report_node(cur_inte->r_son, k, query_point, R_box, rx_prefix, ry_prefix, b - 1, !x_splitter, nn_res);
			}
		}
		else{	// first go right
			if (nn_res.size() < k || r_son_sqrdis <= nn_res.top().second){
				knn_report_node(cur_inte->r_son, k, query_point, R_box, rx_prefix, ry_prefix, b - 1, !x_splitter, nn_res);
			}
			if (nn_res.size() < k || l_son_sqrdis <= nn_res.top().second){
				knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix, y_prefix, b - 1, !x_splitter, nn_res);
			}
		}
		return;
	}

}

// --- From mvzd/src/mvq/versioning.hpp ---
namespace mvq {
	auto Tree::commit(shared_ptr<BaseNode> &old_version, sequence<Point> &P_insert, sequence<Point> &P_delete){
		// parlay::internal::timer t("zdtree breakdown", true);
		auto P_set = get_sorted_points(P_delete);
		auto tmp_ver = multi_version_batch_delete_sorted(P_set, old_version);
		P_set = get_sorted_points(P_insert);
		auto new_ver = multi_version_batch_insert_sorted(P_set, tmp_ver);
		gc_nocheck(tmp_ver);

		// cout << "base version size: " << collect_records(old_version).size() << endl;
		// auto new_ver = multi_version_batch_delete_sorted(P_set, old_version);
		// t.next("delete time");
		// P_set = get_sorted_points(P_insert);
		// new_ver = multi_version_batch_insert_sorted(P_set, new_ver);
		// t.next("insert time");

		// t.next("garbage-collect time");
		// cout << "base version size: " << collect_records(old_version).size() << endl;

		return new_ver;
	}

	template<typename DIFF>
	auto Tree::leaf_leaf_diff(shared_ptr<BaseNode> &leaf1, shared_ptr<BaseNode> &leaf2, Bounding_Box &query_mbr, DIFF &ret_diff){
		auto cur_leaf1 = static_cast<LeafNode*>(leaf1.get());
		auto cur_leaf2 = static_cast<LeafNode*>(leaf2.get());
		auto f = [&](auto pt){
			return point_in_mbr(pt, query_mbr);
		};
		return geobase::merge_pts(cur_leaf1->records, cur_leaf2->records, f, ret_diff);
	}

	auto Tree::leaf_leaf_diff(shared_ptr<BaseNode> &leaf1, shared_ptr<BaseNode> &leaf2){
		// sequence<Point> add = {}, remove = {};

		auto cur_leaf1 = static_cast<LeafNode*>(leaf1.get());
		auto cur_leaf2 = static_cast<LeafNode*>(leaf2.get());

		return merge_pts(cur_leaf1->records, cur_leaf2->records);
	}

	//	note, this returns {add, remove} w.r.t. leaf, pay attention to the order (especially when swapped) 
	//	leaf must not be nullptr
	auto Tree::leaf_inte_diff(std::span<Point> P, shared_ptr<BaseNode> &x, size_t b){
		if (!x){	//	inte could be nullptr if only one branch exists
			sequence<Point> add = {}, remove = parlay::to_sequence(P);
			return make_tuple(add, remove);
		}
		if (x->is_leaf()){	//	Base case when we meet two leaf nodes
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			return merge_pts(P, cur_leaf->records);
		}
		//	divide into two subcases
		auto splitter = split_by_bit(P, 0, P.size(), b);
		auto cur_inte = static_cast<InteNode*>(x.get());
		sequence<Point> L_add, L_remove, R_add, R_remove;

		auto diff_left = [&](){ tie(L_add, L_remove) = leaf_inte_diff(P.subspan(0, splitter), cur_inte->l_son, b - 1); };
		auto diff_right = [&](){ tie(R_add, R_remove) = leaf_inte_diff(P.subspan(splitter), cur_inte->r_son, b - 1); };

		diff_left();
		diff_right();
		L_add.append(R_add);
		L_remove.append(R_remove);

		return make_tuple(L_add, L_remove);
	}


	template<typename DIFF>
	auto Tree::leaf_inte_diff(std::span<Point> P, shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, 
								FT x_prefix, FT y_prefix, size_t b, bool x_splitter, DIFF &ret_diff, bool reverse){
		auto flag = mbr_mbr_relation(cur_mbr, query_mbr);
		if (flag < 0) {
			return;	//	no intersection
		}

		auto f = [&](auto x){	//	filter a sequence of points in the query mbr
			return point_in_mbr(x, query_mbr);
		};

		// intersected, need further recursive processing
		if (!x){	//	inte could be nullptr if only one branch exists
			for (auto& pt : P){
				if (f(pt)){
					ret_diff.remove_point(pt, reverse);
					// parlay::assign_uninitialized(ret_diff.remove[ret_diff.remove_cnt++], pt);
				}
			}
			return;
		}

		if (x->is_leaf()){	//	Base case when we meet two leaf nodes
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			return merge_pts(P, cur_leaf->records, f, ret_diff, reverse);
		}
		//	divide into two subcases
		auto splitter = split_by_bit(P, 0, P.size(), b);
		auto cur_inte = static_cast<InteNode*>(x.get());

		auto [L_box, R_box, rx_prefix, ry_prefix] = compute_cur_box(cur_mbr, x_prefix, y_prefix, b, x_splitter);

		sequence<Point> L_add, L_remove, R_add, R_remove;
		leaf_inte_diff(P.subspan(0, splitter), cur_inte->l_son, query_mbr, L_box, x_prefix, y_prefix, b - 1, !x_splitter, ret_diff, reverse);
		leaf_inte_diff(P.subspan(splitter), cur_inte->r_son, query_mbr, R_box, rx_prefix, ry_prefix,  b - 1, !x_splitter, ret_diff, reverse);
		// auto diff_left = [&](){ tie(L_add, L_remove) = leaf_inte_diff(P, l, splitter, cur_inte->l_son, query_mbr, L_box, x_prefix, y_prefix, b - 1, !x_splitter, ret_diff); };
		// auto diff_right = [&](){ tie(R_add, R_remove) = leaf_inte_diff(P, splitter, r, cur_inte->r_son, query_mbr, R_box, rx_prefix, ry_prefix,  b - 1, !x_splitter, ret_diff); };
		// diff_left();
		// diff_right();
		// L_add.append(R_add);
		// L_remove.append(R_remove);
		return;
	}

	auto Tree::diff(shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2, size_t b){
		if (node1 == node2){	//	case0, no difference, all empty or point to the same (shared) node
			// utils::fetch_and_add(&case0, 1);
			sequence<Point> add = {}, remove = {};
			return make_tuple(add, remove);
		}
		if (!node1){	//	case1, node1 is empty, we need to add all points in node2
			// utils::fetch_and_add(&case1, 1);
			sequence<Point> add = {}, remove = {};

			add = collect_records(node2);			

			return make_tuple(add, remove);
		}
		if (!node2){	//	case1, node2 is empty, we need to remove all points in node1
			// utils::fetch_and_add(&case1, 1);
			sequence<Point> add = {}, remove = {};
			remove = collect_records(node1); 
			return make_tuple(add, remove);
		}
		// both are not empty, but they correspond to the same region
		if (node1->is_leaf() && node2->is_leaf()){	// case2, two leafs
			// utils::fetch_and_add(&case2, 1);
			auto [add, remove] = leaf_leaf_diff(node1, node2);
			return make_tuple(add, remove);
		}
		if (node1->is_leaf() && !node2->is_leaf()){	//	case3, a leaf and an inte
			// utils::fetch_and_add(&case3, 1);
			auto P = static_cast<LeafNode*>(node1.get())->records;
			auto [add, remove] = leaf_inte_diff(P, node2, b);
			return make_tuple(add, remove);
		}
		if (!node1->is_leaf() && node2->is_leaf()){	//	case3, an inte and a leaf
			// utils::fetch_and_add(&case3, 1);
			auto P = static_cast<LeafNode*>(node2.get())->records;
			auto [remove, add] = leaf_inte_diff(P, node1, b);	//	inversed due to the second one is leaf
			return make_tuple(add, remove);
		}
		//	case4, two inte
		// utils::fetch_and_add(&case4, 1);
		sequence<Point> L_add, L_remove, R_add, R_remove;
		auto cur_inte1 = static_cast<InteNode*>(node1.get());
		auto cur_inte2 = static_cast<InteNode*>(node2.get());
		auto diff_left = [&](){tie(L_add, L_remove) = diff(cur_inte1->l_son, cur_inte2->l_son, b - 1); };	// both go left
		auto diff_right = [&](){tie(R_add, R_remove) = diff(cur_inte1->r_son, cur_inte2->r_son, b - 1);	};	// both go right

		par_do_if(node1->get_num_points() + node2->get_num_points() >= granularity_cutoff,
			diff_left,
			diff_right);

		L_add.append(R_add);
		L_remove.append(R_remove);
		return make_tuple(L_add, L_remove);
	}

	template<typename DIFF>
	auto Tree::spatial_diff(shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, 
							FT x_prefix, FT y_prefix, size_t b, bool x_splitter, DIFF &ret_diff){
		auto flag = mbr_mbr_relation(cur_mbr, query_mbr);
		if (flag < 0){	// no intersection with the query region, skip
			return;
		}

		//	case overlapped, need to handle multiple cases
		/* Case 0: no difference */
		if (node1 == node2) return;
		/* Case 1: node1 is empty, add all points from node2 */
		if (!node1){
			range_report_node(node2, query_mbr, cur_mbr, x_prefix, y_prefix, b, x_splitter, ret_diff.add_cnt, ret_diff.add);
			return;
		}
		/* Case 1: node2 is empty, remove all points from node1 */
		if (!node2){
			range_report_node(node1, query_mbr, cur_mbr, x_prefix, y_prefix, b, x_splitter, ret_diff.remove_cnt, ret_diff.remove);
			return;
		}
		/* Case 2: two leaf nodes, need to scan both */
		if (node1->is_leaf() && node2->is_leaf()){
			leaf_leaf_diff(node1, node2, query_mbr, ret_diff);
			return;	
		}
		/* Case 3: one of the nodes is leaf */
		if (node1->is_leaf() && !node2->is_leaf()){
			auto P = static_cast<LeafNode*>(node1.get())->records;
			leaf_inte_diff(P, node2, query_mbr, cur_mbr, x_prefix, y_prefix, b, x_splitter, ret_diff);
			return;
		}
		/* Case 3: the other possibility */
		if (!node1->is_leaf() && node2->is_leaf()){
			auto P = static_cast<LeafNode*>(node2.get())->records;
			leaf_inte_diff(P, node1, query_mbr, cur_mbr, x_prefix, y_prefix, b, x_splitter, ret_diff, true);
			return;
		}
		/* Case 4: two interior nodes, not fully overlapped */
		auto cur_inte1 = static_cast<InteNode*>(node1.get());
		auto cur_inte2 = static_cast<InteNode*>(node2.get());

		auto [L_box, R_box, rx_prefix, ry_prefix] = compute_cur_box(cur_mbr, x_prefix, y_prefix, b, x_splitter);
		spatial_diff(cur_inte1->l_son, cur_inte2->l_son, query_mbr, L_box, x_prefix, y_prefix, b - 1, !x_splitter, ret_diff);
		spatial_diff(cur_inte1->r_son, cur_inte2->r_son, query_mbr, R_box, rx_prefix, ry_prefix, b - 1, !x_splitter, ret_diff);

		return;
	}

	auto Tree::two_version_diff(shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2){
		// case0 = 0, case1 = 0, case2 = 0, case3 = 0, case4 = 0;
		// case0_time = 0, case1_time = 0, case2_time = 0, case3_time = 0, case4_time = 0;
		auto [add, remove] = diff(node1, node2, 64);
		return make_tuple(add, remove);
	}
	
	template<typename DIFF>
	auto Tree::spatial_two_version_diff(shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, DIFF &ret_diff){
		return spatial_diff(node1, node2, query_mbr, cur_mbr, 0.0, 0.0, 64, true, ret_diff);
	}

	//	current return conflict only
	auto Tree::merge(shared_ptr<BaseNode> &base, shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2){
		auto [add1, remove1] = two_version_diff(base, node1);
		auto [add2, remove2] = two_version_diff(base, node2); 

		auto [insert1, delete1, update1] = filter_diff_results(add1, remove1);
		auto [insert2, delete2, update2] = filter_diff_results(add2, remove2);

		auto [no_conflict_insert, conflict_insert] = merge_by_id_with_conflict(insert1, insert2);	// insert conflict
		auto [no_conflict_update, conflict_update] = merge_by_id_with_conflict(update1, update2);	// update conflict
		//	delate = delete + update
		auto [no_conflict_delete, conflict_delete] = merge_by_id_with_conflict(delete1, update2);	//	delete and update conflict
		auto [no_conflict_delete1, conflict_delete1] = merge_by_id_with_conflict(delete2, update1);	//	delete and update conflict, another case

		no_conflict_delete.append(no_conflict_delete1);
		conflict_delete.append(conflict_delete1);

		auto new_ver = commit(base, no_conflict_insert, no_conflict_delete);

		return make_tuple(new_ver, conflict_insert, conflict_update, conflict_delete);
	}

	// spatial join w.r.t. function F
	template<typename JOIN_RES>
	void Tree::spatial_join(shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2, FT &point_dis, JOIN_RES &join_res,
							Bounding_Box &mbr1, FT x_prefix1, FT y_prefix1, size_t b1, bool x_splitter1,
							Bounding_Box &mbr2, FT x_prefix2, FT y_prefix2, size_t b2, bool x_splitter2){
		
		if (!node1 || !node2) return;

		if (!mbr_mbr_within_dis(mbr1, mbr2, point_dis)) return;
		// Base Case
		if (node1->is_leaf() && node2->is_leaf()){
			auto leaf1 = static_cast<LeafNode*>(node1.get());
			auto leaf2 = static_cast<LeafNode*>(node2.get());
			for (auto &pt1: leaf1->records){
				for (auto &pt2: leaf2->records){
					if (pt1.id != pt2.id && dcmp(point_point_sqrdis(pt1, pt2) - point_dis * point_dis) <= 0){
						join_res.emplace_back(pt1, pt2);
					}
				}
			}
			return;
		}

		if (!node1->is_leaf() && !node2->is_leaf()){
			auto cur_inte1 = static_cast<InteNode*>(node1.get());
			auto cur_inte2 = static_cast<InteNode*>(node2.get());
	
			auto [L_box1, R_box1, rx_prefix1, ry_prefix1] = compute_cur_box(mbr1, x_prefix1, y_prefix1, b1, x_splitter1);
			auto [L_box2, R_box2, rx_prefix2, ry_prefix2] = compute_cur_box(mbr2, x_prefix2, y_prefix2, b2, x_splitter2);
			/* Case LL */
			spatial_join(cur_inte1->l_son, cur_inte2->l_son, point_dis, join_res, 
						L_box1, x_prefix1, y_prefix1, b1 - 1, !x_splitter1, 
						L_box2, x_prefix2, y_prefix2, b2 - 1, !x_splitter2);
			/* Case LR */
			spatial_join(cur_inte1->l_son, cur_inte2->r_son, point_dis, join_res, 
				L_box1, x_prefix1, y_prefix1, b1 - 1, !x_splitter1, 
				R_box2, rx_prefix2, ry_prefix2, b2 - 1, !x_splitter2);
			/* Case RL*/
			spatial_join(cur_inte1->r_son, cur_inte2->l_son, point_dis, join_res, 
				R_box1, rx_prefix1, ry_prefix1, b1 - 1, !x_splitter1, 
				L_box2, x_prefix2, y_prefix2, b2 - 1, !x_splitter2);
			/* Case RR*/
			spatial_join(cur_inte1->r_son, cur_inte2->r_son, point_dis, join_res, 
				R_box1, rx_prefix1, ry_prefix1, b1 - 1, !x_splitter1, 
				R_box2, rx_prefix2, ry_prefix2, b2 - 1, !x_splitter2);
		}
		else if (!node1->is_leaf()){
			auto cur_inte1 = static_cast<InteNode*>(node1.get());
			auto [L_box1, R_box1, rx_prefix1, ry_prefix1] = compute_cur_box(mbr1, x_prefix1, y_prefix1, b1, x_splitter1);
			/* Case L-cur_Leaf */
			spatial_join(cur_inte1->l_son, node2, point_dis, join_res, 
				L_box1, x_prefix1, y_prefix1, b1 - 1, !x_splitter1, 
				mbr2, x_prefix2, y_prefix2, b2, x_splitter2);
			/* Case R-cur_leaf */
			spatial_join(cur_inte1->r_son, node2, point_dis, join_res, 
				R_box1, rx_prefix1, ry_prefix1, b1 - 1, !x_splitter1, 
				mbr2, x_prefix2, y_prefix2, b2, x_splitter2);
		}
		else if (!node2->is_leaf()){
			auto cur_inte2 = static_cast<InteNode*>(node2.get());
			auto [L_box2, R_box2, rx_prefix2, ry_prefix2] = compute_cur_box(mbr2, x_prefix2, y_prefix2, b2, x_splitter2);
			/* Case cur_Leaf-L */
			spatial_join(node1, cur_inte2->l_son, point_dis, join_res, 
				mbr1, x_prefix1, y_prefix1, b1, x_splitter1,	
				L_box2, x_prefix2, y_prefix2, b2 - 1, !x_splitter2);
			/* Case cur_leaf-R */
			spatial_join(node1, cur_inte2->r_son, point_dis, join_res, 
				mbr1, x_prefix1, y_prefix1, b1, x_splitter1,	
				R_box2, rx_prefix2, ry_prefix2, b2 - 1, !x_splitter2);			
		}
		return;
	}

	template<typename JOIN_RES>
	void Tree::two_version_spatial_join(shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2, FT &point_dis, JOIN_RES &join_res,
									Bounding_Box mbr1, Bounding_Box mbr2){
		spatial_join(node1, node2, point_dis, join_res,
						mbr1, 0.0, 0.0, 64, true,
						mbr2, 0.0, 0.0, 64, true);
	}
	void Tree::merge_nodes(shared_ptr<BaseNode> &L, shared_ptr<BaseNode> &R, shared_ptr<InteNode> &cur_node){
		// deal with MBR, covered points of parent
		auto L_num_pts = L == nullptr ? 0 : L->get_num_points();
		auto R_num_pts = R == nullptr ? 0 : R->get_num_points();

		cur_node->num_pts = L_num_pts + R_num_pts;
		cur_node->l_son = move(L);
		cur_node->r_son = move(R);
		// if (L == nullptr) {
		// 	cur_node->num_pts = R->get_num_points();
		// }
		// else if (R == nullptr){
		// 	cur_node->num_pts = L->get_num_points();
		// } 
		// else {
		// 	cur_node->num_pts = L->get_num_points() + 
		// 						R->get_num_points();
		// }
		// cur_node->l_son = move(L);
		// cur_node->r_son = move(R);
	}

	void Tree::multi_version_batch_insert_sorted_node(shared_ptr<BaseNode> &x, std::span<Point> P, size_t b){
		if (x == nullptr){
			x = build(P, b);
			return;
		}
		auto less = [&](auto lhs, auto rhs){
			return lhs.morton_id < rhs.morton_id ||
					(lhs.morton_id == rhs.morton_id && lhs.id < rhs.id);
		};
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			if (!b || cur_leaf->records.size() + P.size() <= leaf_size){	// current leaf is not full
				cur_leaf->records = parlay::merge(cur_leaf->records, parlay::make_slice(P), less);
				cur_leaf->update_mem();
				return;
			}
			else{
				auto new_points = parlay::merge(cur_leaf->records, parlay::make_slice(P), less);
				x = build(std::span<Point>(new_points.begin(), new_points.end()), b);
				return;
			}
		}
		auto splitter = split_by_bit(P, 0, P.size(), b);
		
		auto cur_inte = static_cast<InteNode*>(x.get());
		auto new_l_son(cur_inte->l_son);
		auto new_r_son(cur_inte->r_son);

		auto insert_left = [&](){ 
			if (0 < splitter){ 
				new_l_son = node_copy(cur_inte->l_son);
				multi_version_batch_insert_sorted_node(new_l_son, P.subspan(0, splitter), b - 1);
			};
		};
		auto insert_right = [&](){ 
			if (splitter < P.size()){ 
				new_r_son = node_copy(cur_inte->r_son);
				multi_version_batch_insert_sorted_node(new_r_son, P.subspan(splitter), b - 1); 
			};
		};
		par_do_if(P.size() >= 256,
			insert_left,
			insert_right); 
			
		delete_merge_nodes(new_l_son, new_r_son, cur_inte);
		// delete_merge_nodes(new_l_son, new_r_son, cur_inte);
	}

	shared_ptr<BaseNode> Tree::multi_version_batch_insert_sorted(std::span<Point> P, shared_ptr<BaseNode> &old_root){
		auto new_root = node_copy(old_root);
		multi_version_batch_insert_sorted_node(new_root, P, 64);
		// multi_version_roots.emplace_back(new_root);
		return new_root;
	}

	void Tree::multi_version_batch_delete_sorted_node(shared_ptr<BaseNode> &x, std::span<Point> P, size_t b){
		if (x == nullptr){
			return;
		}
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			auto tmp_records = parlay::to_sequence(P);
			cur_leaf->records = get_delete_p(cur_leaf->records, tmp_records, 0, tmp_records.size());
			cur_leaf->update_mem();
			if (!cur_leaf->records.size()) x.reset();
			return;
		}
		// visited_inte++;

		auto splitter = split_by_bit(P, 0, P.size(), b);
		auto cur_inte = static_cast<InteNode*>(x.get());
		auto new_l_son(cur_inte->l_son);
		auto new_r_son(cur_inte->r_son);

		auto delete_left = [&](){ 
			if (0 < splitter){ 
				new_l_son = node_copy(cur_inte->l_son);
				multi_version_batch_delete_sorted_node(new_l_son, P.subspan(0, splitter), b - 1); 
			}; 
		};
		auto delete_right = [&](){
			if (splitter < P.size()){ 
				new_r_son = node_copy(cur_inte->r_son);
				multi_version_batch_delete_sorted_node(new_r_son, P.subspan(splitter), b - 1); 
			}; 
		};

		par_do_if(P.size() >= 256,
			delete_left,
			delete_right); 

		auto less = [&](auto lhs, auto rhs){
			return lhs.morton_id < rhs.morton_id ||
					(lhs.morton_id == rhs.morton_id && lhs.id < rhs.id);
		};

		cur_inte->l_son = new_l_son;
		cur_inte->r_son = new_r_son;
		if (!new_l_son && !new_r_son) x.reset();
		else{
			if (!new_l_son) {
				if (new_r_son->get_num_points() <= leaf_size) x = move(new_r_son);
			}
			else {
				if (!new_r_son) {
					if (new_l_son->get_num_points() <= leaf_size) x = move(new_l_son);
				}
				else{
			 		if (new_l_son->get_num_points() + new_r_son->get_num_points() <= leaf_size){
						auto L = static_cast<LeafNode*>(new_l_son.get());
						auto R = static_cast<LeafNode*>(new_r_son.get());
						auto cur_records = parlay::merge(L->records, R->records, less);
						x = create_leaf(std::span<Point>(cur_records.begin(), cur_records.end()), 0);
						x = create_leaf(std::span<Point>(cur_records.begin(), cur_records.end()), 0);
					}
					else{
						delete_merge_nodes(new_l_son, new_r_son, cur_inte);
					}
				}
			}
			auto L_num_pts = cur_inte->l_son == nullptr ? 0 : cur_inte->l_son->get_num_points();
			auto R_num_pts = cur_inte->r_son == nullptr ? 0 : cur_inte->r_son->get_num_points();
			cur_inte->num_pts = L_num_pts + R_num_pts;
		}
	}

	shared_ptr<BaseNode> Tree::multi_version_batch_delete_sorted(std::span<Point> P, shared_ptr<BaseNode> &old_root){
		auto new_root = node_copy(old_root);
		// cout << "copy finished" << endl;
		multi_version_batch_delete_sorted_node(new_root, P, 64);
		// multi_version_roots.emplace_back(new_root);
		return new_root;
	}

}
