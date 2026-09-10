#pragma once 

#include <bits/stdc++.h>

#include <cpam/cpam.h>
// #include <pam/pam.h>
#include <parlay/primitives.h>
#include "geo/point.hpp"
#include "geo/operations.hpp"
#include "geo/io.hpp"
#include <parlay/internal/get_time.h>
#include "pam/utils.h"

#define BR_MBR
#define SEQ


namespace CPAMBB{
	using namespace std;
	using namespace geobase;
	using parlay::sequence;
	using parlay::par_do;
	using parlay::par_do_if;

	using key_type = pair<unsigned long long, unsigned long long>;	// morton_id, id
	using val_type = Point;
	using aug_type = pair<Bounding_Box, size_t>;

	//	CPAM entry
	struct entry {
		using key_t = key_type;
		using val_t = val_type;
		using aug_t = aug_type;

		static inline bool comp(key_t a, key_t b) { return a < b; }
		static aug_t get_empty() { return make_pair(Bounding_Box{Point(1e60, 1e60), Point(-1, -1)}, 0); }
		static aug_t from_entry(key_t k, val_t v) { return make_pair(Bounding_Box(v, v), 1); }
		static aug_t combine(aug_t a, aug_t b) { return make_pair(merge_mbr(a.first, b.first), a.second + b.second); }
	};

	// using zmap = cpam::aug_map<entry, 32>;
	using zmap = cpam::aug_map<entry, 16>;
	using par = std::tuple<entry::key_t, entry::val_t>;

	template<class T, class MBR>
	auto filter_range(T &tree, MBR query_mbr, bool use_hilbert = false){
		auto BL = query_mbr.first;
		auto UR = query_mbr.second;

		auto fpt = [&](par cur){ 
			if (BL.x <= get<1>(cur).x && get<1>(cur).x <= UR.x &&
				BL.y <= get<1>(cur).y && get<1>(cur).y <= UR.y){
				return true;
			}
			return false;
		};

		auto fbb = [&](auto cur){ 
			return !mbr_exclude_mbr(query_mbr, cur.first);
		};

		// return zmap::aug_filter(tree, fbb);
		// return  zmap::filter(zmap::aug_filter(tree, fbb), fpt);
		return zmap::aug_filter2(tree, fpt, fbb);
	}

	template<typename M, typename DIFF>
	auto map_spatial_diff(M &lhs, M &rhs, Bounding_Box &query, DIFF &ret_diff){
		auto filtered_lhs = filter_range(lhs, query);
		auto filtered_rhs = filter_range(rhs, query);
		ret_diff.add = zmap::values(zmap::map_difference(filtered_rhs, filtered_lhs));
		ret_diff.add_cnt = ret_diff.add.size();
		ret_diff.remove = zmap::values(zmap::map_difference(filtered_lhs, filtered_rhs));
		ret_diff.remove_cnt = ret_diff.remove.size();
		return;
	}

	template<typename M>
	auto map_diff(M &lhs, M &rhs){
		auto add = zmap::values(zmap::map_difference(rhs, lhs));
		auto remove = zmap::values(zmap::map_difference(lhs, rhs));
		return make_tuple(add, remove);
	}

	template<class PT>
	auto map_init(PT &P, bool use_hilbert = false){
		size_t n = P.size();
		
		parlay::parallel_for(0, n, [&](int i){
			P[i].morton_id = use_hilbert ? P[i].overlap_bits() : P[i].interleave_bits();
		});
		// auto P_set = get_sorted_points(P);
		parlay::sequence<par> entries(n);
		parlay::parallel_for(0, n, [&](int i){
			entries[i] = {{P[i].morton_id, P[i].id}, P[i]};
			// entries[i] = {P[i]->id, P[i]};
		});
		zmap m1(entries);
		// auto vals = zmap::values(m1);
		// for (auto i = 1; i < vals.size(); i++){
		// 	if (vals[i].morton_id < vals[i - 1].morton_id){
		// 		cout << "[ERROR]: not sorted" << endl;
		// 	}
		// }
		return m1;
	}	

	template<typename PT, typename M>
	auto map_insert(PT &P, M &mmp, bool use_hilbert = false){
		size_t n = P.size();

		parlay::internal::timer t("debug", false);
		parlay::parallel_for(0, n, [&](int i){
			P[i].morton_id = use_hilbert ? P[i].overlap_bits() : P[i].interleave_bits();
		});
		auto insert_pts = parlay::sequence<par>::uninitialized(n);
		parlay::parallel_for(0, n, [&](int i){
			insert_pts[i] = {{P[i].morton_id, P[i].id}, P[i]};
			// insert_pts[i] = {P[i]->id, P[i]};
		});
		t.next("init time");
		auto m2 = zmap::multi_insert(mmp, insert_pts);
		t.next("insert time");

		return m2;
	}

	template<typename PT, typename M>
	auto map_delete(PT &P, M &mmp, bool use_hilbert = false){
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

	template<class T, class MBR>
	auto range_count(T &zCPAM, MBR &query_mbr, bool use_hilbert = false){
		auto f = [&](auto cur){ 
			return mbr_mbr_relation(cur, query_mbr);
		};

		auto f2 = [&](auto cur){ 
			return point_in_mbr(cur, query_mbr);
		};

		// auto res = zmap::range_count_filter(zCPAM, f, f2);
		auto res = zmap::range_count_filter2(zCPAM, f, f2);
		return res;
	}

	template<class T, class MBR>
	auto range_report(T &tree, MBR &query_mbr, parlay::sequence<Point> &out, bool use_hilbert = false){
		// auto ret = zmap::values(filter_range(tree, query_mbr, use_hilbert));
		auto f = [&](auto &cur){ 
			return mbr_mbr_relation(cur, query_mbr);
		};

		int64_t ret = 0;
		// for (auto i = 0; i < 20; i++){
			// ret = 0;
			// zmap::range_report_filter(tree, f, ret, out);
			zmap::range_report_filter2(tree, f, ret, out);
		// }
		return ret;
	}

	template<typename M, typename DIFF>
	auto plain_map_spatial_diff(M &lhs, M &rhs, Bounding_Box &query, DIFF &ret_diff, parlay::sequence<Point> &l_pts, parlay::sequence<Point> &r_pts){
		// auto l_pts = zmap::values(filter_range(lhs, query));
		// auto r_pts = zmap::values(filter_range(rhs, query));

		auto ret1 = range_report(lhs, query, l_pts);
		auto ret2 = range_report(rhs, query, r_pts);
		l_pts.resize(ret1);
		r_pts.resize(ret2);
		// print_Pset_info(l_pts, "lpts");
		// for (auto &pt: l_pts){
		// 	cout << pt.morton_id << endl;
		// }
		// print_Pset_info(r_pts, "rpts");
		// for (auto &pt: r_pts){
		// 	cout << pt.morton_id << endl;
		// }
		// cout << "debug: " << l_pts.size() << ", " << r_pts.size() << endl;
		// auto [add, remove] = merge_pts(l_pts, r_pts);
		// cout << "l size = " << l_pts.size() << endl;
		// cout << "r size = " << r_pts.size() << endl;
		
		merge_pts(l_pts, r_pts, ret_diff);
		
		// geobase::print_Pset_info(add, "add");
		// geobase::print_Pset_info(remove, "remove");
		return;
	}

	template<typename T>
	auto knn(T &tree, geobase::Point &query_point, size_t &k){
		auto f = [&](auto cur_pt){ return point_point_sqrdis(cur_pt, query_point); };

		auto f2 = [&](auto cur_mbr){ return point_mbr_sqrdis(query_point, cur_mbr); };

		priority_queue<nn_pair, vector<nn_pair>, nn_pair_cmp> nn_res;
		zmap::knn_filter(tree, f, f2, k, nn_res);
		return nn_res;
	}

	template<typename PT, typename M>
	auto map_commit(M &mmp, PT &P_insert, PT &P_delete){
		// parlay::internal::timer t("CPAM-BB breakdown", true);
		auto new_ver = map_delete(P_delete, mmp);	//	new	version
		// t.next("delete time");
		new_ver = map_insert(P_insert, new_ver); 
		// t.next("insert time");
		return new_ver;
	}

	template<typename M>
	auto map_merge(M &base, M &v1, M &v2){
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

		auto new_ver = map_commit(base, no_conflict_insert, no_conflict_delete);
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

	auto print_stats(){
		zmap::GC::print_stats();
	}

}

namespace BRTree{

	using namespace std;
	using namespace geobase;
	using parlay::sequence;
	using parlay::par_do;
	using parlay::par_do_if;

	struct BaseNode{
		#ifdef BR_MBR
		Bounding_Box mbr;
		#endif
		BaseNode(){}

		virtual~BaseNode() = default;
		virtual bool is_leaf(){ return false; }
		virtual size_t get_num_points(){ return 0; }
	};

	struct InteNode: BaseNode{
		shared_ptr<BaseNode> l_son, r_son;
		size_t num_pts;

		InteNode(): l_son(nullptr), r_son(nullptr), num_pts(0){}
		// InteNode(InteNode &x): l_son(x->l_son), r_son(x->r_son), num_pts(x->num_pts){}
		
		virtual bool is_leaf(){ return false; }
		virtual size_t get_num_points(){ return num_pts; }
	};

	struct LeafNode: BaseNode{
		sequence<Point> records = {};
		// parlay::slice<Point*, Point*> records {nullptr, nullptr};

		template<typename Records>
		LeafNode(Records &r){
			records.resize(r.size());
			for (size_t i = 0; i < r.size(); i++){
				records[i] = r[i];
			}
		}

		// LeafNode(LeafNode &x): records(x->records){}

		virtual bool is_leaf(){ return true; }
		virtual size_t get_num_points(){ return records.size(); }

		#ifdef BR_MBR
		void print_mbr(){ 
			cout << "[(" << mbr.first.x << ", " << mbr.first.y << "), (" << mbr.second.x << ", " << mbr.second.y << ")" << "]" << endl;
		}
		#endif

		void print_records(){
			cout << records.size() << endl;
			for (size_t i = 0; i < records.size(); i++){
				cout << "(" << records[i].x << ", " << records[i].y << ")" << endl; 
			}
		}
	};

	class Tree{
	private:
		size_t visited_leaf = 0;
		size_t visited_inte = 0;
		size_t granularity_cutoff = 1000;

	public:
		size_t height = 0;
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
		// ~Tree();

		shared_ptr<BaseNode> node_copy(shared_ptr<BaseNode> &x);

		// entrance of building zdtree & tree construction
		shared_ptr<BaseNode> build(sequence<Point> &P, size_t l, size_t r, size_t b);
		void build(sequence<Point> &P);
		void clear();

		void in_order_traverse(shared_ptr<BaseNode> &x);
		void in_order_traverse(shared_ptr<BaseNode> &x, unsigned x_prefix, unsigned y_prefix, size_t b, bool x_splitter);
		int  tree_hash(shared_ptr<BaseNode> &x);
		void tree_compress(shared_ptr<BaseNode> &x);

		void tree_node_count(shared_ptr<BaseNode> &x, set<uintptr_t> &s);	//	store the nodes of a root in a set s
		size_t node_count();	//	count the total number of nodes among all versions

		void merge_nodes(shared_ptr<BaseNode> &lhs, shared_ptr<BaseNode> &rhs, shared_ptr<InteNode> &cur); 	// merge two sons to current node.
		void delete_merge_nodes(shared_ptr<BaseNode> &L, shared_ptr<BaseNode> &R, InteNode* cur_node);

		shared_ptr<InteNode> create_internal(shared_ptr<BaseNode> &L, shared_ptr<BaseNode> &R);	//	create an internal node, do not store pointers to original records.
		shared_ptr<LeafNode> create_leaf(sequence<Point> &P, size_t l, size_t r, size_t b); // create a leaf, store all (pointers of) records.

		// size_t num_covered_points(shared_ptr<BaseNode> x);

		// intervals are left close right open: i.e., [L, R)
		// 2D range report
		#ifdef BR_MBR
		template <class Out> 
		void range_report_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, size_t &cnt, Out &out);
		template <class Out>
		void range_report(Bounding_Box &query_mbr, size_t &cnt, Out &out);
		#endif

		template <class Out> 
		void range_report_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter, size_t &cnt, Out &out);
		template <class Out>
		void range_report(Bounding_Box &query_mbr, Bounding_Box &cur_mbr, size_t &cnt, Out &out);

		// 2D range count
		#ifdef BR_MBR
		size_t range_count_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr);
		size_t range_count(Bounding_Box &query_mbr);
		size_t range_count_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, bool x_splitter);
		#endif

		size_t range_count_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter);
		size_t range_count(Bounding_Box &query_mbr, Bounding_Box &cur_mbr);

		// 2D k nearest neighbor report
		#ifdef BR_MBR
		template<class T> 
		void knn_report_node(shared_ptr<BaseNode> &x, size_t &k, Point query_point, T &nn_res);
		auto knn_report(size_t &k, Point query_point);
		#endif

		template<class T> 
		void knn_report_node(shared_ptr<BaseNode> &x, size_t &k, Point query_point, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter, T &nn_res);
		auto knn_report(size_t &k, Point query_point, Bounding_Box &cur_mbr);
	};
		
	Tree::Tree(size_t _leaf_sz){
		leaf_size = _leaf_sz;
	}

	void Tree::clear(){
		root.reset();
	}

	void Tree::delete_merge_nodes(shared_ptr<BaseNode> &L, shared_ptr<BaseNode> &R, InteNode* cur_node){
		// deal with MBR, covered points of parent
		if (L == nullptr) {
			#ifdef BR_MBR
			cur_node->mbr = R->mbr;
			#endif
			cur_node->num_pts = R->get_num_points();
		}
		else if (R == nullptr){
			#ifdef BR_MBR
			cur_node->mbr = L->mbr;
			#endif
			cur_node->num_pts = L->get_num_points();
		} 
		else {
			#ifdef BR_MBR
			cur_node->mbr = merge_mbr(L->mbr, R->mbr);
			#endif
			cur_node->num_pts = L->get_num_points() + 
								R->get_num_points();
		}
		cur_node->l_son = move(L);
		cur_node->r_son = move(R);
		// cur_node->l_son = L;
		// cur_node->r_son = R;
	}

	void Tree::merge_nodes(shared_ptr<BaseNode> &L, shared_ptr<BaseNode> &R, shared_ptr<InteNode> &cur_node){
		// deal with MBR, covered points of parent
		if (L == nullptr) {
			#ifdef USE_MBR
			cur_node->mbr = R->mbr;
			#endif
			cur_node->num_pts = R->get_num_points();
		}
		else if (R == nullptr){
			#ifdef BR_MBR
			cur_node->mbr = L->mbr;
			#endif
			cur_node->num_pts = L->get_num_points();
		} 
		else {
			#ifdef BR_MBR
			cur_node->mbr = merge_mbr(L->mbr, R->mbr);
			#endif
			cur_node->num_pts = L->get_num_points() + 
								R->get_num_points();
		}
		cur_node->l_son = move(L);
		cur_node->r_son = move(R);
	}

	shared_ptr<InteNode> Tree::create_internal(shared_ptr<BaseNode> &L, shared_ptr<BaseNode> &R){
		shared_ptr<InteNode> cur_node(new InteNode());
		
		// augmented changes happen here
		merge_nodes(L, R, cur_node);

		return cur_node;		
	}

	shared_ptr<LeafNode> Tree::create_leaf(sequence<Point> &P, size_t l, size_t r, size_t b){
		auto cur_records = parlay::make_slice(&P[l], &P[r]);
		shared_ptr<LeafNode> cur_node(new LeafNode(cur_records));

		#ifdef BR_MBR
		cur_node->mbr = get_mbr(cur_records);
		#endif

		return cur_node;
	}

	shared_ptr<BaseNode> Tree::build(sequence<Point> &P, size_t l, size_t r, size_t b){
		if (r - l <= leaf_size){
			height = max(height, 64 - b);
			return create_leaf(P, l, r, b);
		}
		// auto splitter = split_by_bit(P, l, r, b);
		auto splitter = (r + l) >> 1;   //  build a BST
		shared_ptr<BaseNode> L = nullptr;
		shared_ptr<BaseNode> R = nullptr;
		#ifdef SEQ
		if (l < splitter) L = build(P, l, splitter, b - 1);
		if (splitter < r) R = build(P, splitter, r, b - 1);
		#else
		auto build_left = [&](){ if (l < splitter) L = build(P, l, splitter, b - 1); }; 
		auto build_right = [&](){ if (splitter < r) R = build(P, splitter, r, b - 1); }; 
		par_do_if(r - l >= granularity_cutoff,
			build_left, 
			build_right);
		#endif
		return create_internal(L, R);
	}

	void Tree::build(sequence<Point> &P){
		if (!P.size()) return;
		root = build(P, 0, P.size(), 64);
	}

	void Tree::tree_compress(shared_ptr<BaseNode> &x){
		if (x->is_leaf()){
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		if (cur_inte->l_son != nullptr) tree_compress(cur_inte->l_son);
		if (cur_inte->r_son != nullptr) tree_compress(cur_inte->r_son);

		if (cur_inte->l_son == nullptr && cur_inte->r_son == nullptr) x.reset();
		else if (cur_inte->l_son == nullptr) x = move(cur_inte->r_son);
		else if (cur_inte->r_son == nullptr) x = move(cur_inte->l_son);
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

	//	store all pointers of a root in a set
	void Tree::tree_node_count(shared_ptr<BaseNode> &x, set<uintptr_t> &s){
		if (!x) return;
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			auto node_addr = reinterpret_cast<uintptr_t>(cur_leaf); 
			if (s.find(node_addr) != s.end()) return;
			s.insert(node_addr);
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		auto node_addr = reinterpret_cast<uintptr_t>(cur_inte);
		if (s.find(node_addr) != s.end()) return;
		s.insert(node_addr);
		tree_node_count(cur_inte->l_son, s);
		tree_node_count(cur_inte->r_son, s);
	}

	size_t Tree::node_count(){
		set<uintptr_t > s;
		// cout << "set size = " << multi_version_roots.size() << endl;
		for (auto &root: multi_version_roots){
			tree_node_count(root, s);
		}
		return s.size();
	}

	void Tree::in_order_traverse(shared_ptr<BaseNode> &x){
		#ifdef BR_MBR
		print_mbr(x->mbr);
		#endif
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			if (cur_leaf->get_num_points() > leaf_size){
			// 	print_mbr(x->mbr);
			// 	cout << cur_leaf->get_num_points() << endl;
			}
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		if (cur_inte->l_son != nullptr) in_order_traverse(cur_inte->l_son);
		if (cur_inte->r_son != nullptr) in_order_traverse(cur_inte->r_son);
	}
	
	void Tree::in_order_traverse(shared_ptr<BaseNode> &x, unsigned x_prefix, unsigned y_prefix, size_t b, bool x_splitter){
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			cur_leaf->print_records();
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		if (x_splitter){
			if (cur_inte->l_son != nullptr) {
				in_order_traverse(cur_inte->l_son, x_prefix, y_prefix, b, !x_splitter);
			}
			if (cur_inte->r_son != nullptr){
				in_order_traverse(cur_inte->r_son, x_prefix + (1 << (b - 1)), y_prefix, b, !x_splitter);
			}
		}
		else{
			if (cur_inte->l_son != nullptr){
				in_order_traverse(cur_inte->l_son, x_prefix, y_prefix, b - 1, !x_splitter);
			}
			if (cur_inte->r_son != nullptr){
 				in_order_traverse(cur_inte->r_son, x_prefix, y_prefix + (1 << (b - 1)), b - 1, !x_splitter);
			}
		}
	}
	
// range count using mbr
#ifdef BR_MBR
	size_t Tree::range_count_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr){
		// visited_inte++;
		int flag = mbr_mbr_relation(x->mbr, query_mbr);
		if (flag < 0) return 0;
		if (flag > 0) return x->get_num_points();
		if (x->is_leaf()){	// we have to scan the leaf to report the number of points;
			// visited_leaf++;
			size_t ret = 0;
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (auto &p: cur_leaf->records){
				if (point_in_mbr(p, query_mbr)){
					ret += 1;
				}
			}
			return ret;
		}
		else{
			auto cur_inte = static_cast<InteNode*>(x.get());
			size_t ret_L = 0, ret_R = 0;
			if (cur_inte->l_son != nullptr) ret_L = range_count_node(cur_inte->l_son, query_mbr);
			if (cur_inte->r_son != nullptr) ret_R = range_count_node(cur_inte->r_son, query_mbr);
			return ret_L + ret_R;
		}
		return -1; // unexpected error happens if the code runs to here.
	}

	size_t Tree::range_count(Bounding_Box &query_mbr){
		// visited_inte = 0, visited_leaf = 0;
		size_t ret = range_count_node(root, query_mbr);		
		// visited_inte -= visited_leaf;
		// cout << "visited nodes (interior/leaf): " << visited_inte << ", " << visited_leaf << endl;
		return ret;
	}

	size_t Tree::range_count_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, bool x_splitter){
		int flag = mbr_mbr_relation(x->mbr, query_mbr);
		if (flag < 0) return 0;
		if (flag > 0) return x->get_num_points();
		if (x->is_leaf()){	// we have to scan the leaf to report the number of points;
			size_t ret = 0;
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (auto &p: cur_leaf->records){
				if (point_in_mbr(p, query_mbr)){
					ret += 1;
				}
			}
			return ret;
		}
		else{
			auto cur_inte = static_cast<InteNode*>(x.get());
			size_t ret_L = 0, ret_R = 0;
			if (cur_inte->l_son != nullptr) {
				auto L_box = cur_mbr;
				if (x_splitter) L_box.second.x = cur_inte->l_son->mbr.second.x;
				else L_box.second.y = cur_inte->l_son->mbr.second.y;
				ret_L = range_count_node(cur_inte->l_son, query_mbr, L_box, !x_splitter);
			}
			if (cur_inte->r_son != nullptr) {
				auto R_box = cur_mbr;
				if (x_splitter) R_box.first.x = cur_inte->r_son->mbr.first.x;
				else R_box.first.y = cur_inte->r_son->mbr.first.y;
				ret_R = range_count_node(cur_inte->r_son, query_mbr, R_box, !x_splitter);
			}
			return ret_L + ret_R;
		}
		return -1; // unexpected error happens if the code runs to here.
	}
#endif

	size_t Tree::range_count_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter){
		int flag = mbr_mbr_relation(cur_mbr, query_mbr);
		if (flag < 0) return 0;
		if (flag > 0) return x->get_num_points();
		if (x->is_leaf()){	// we have to scan the leaf to report the number of points;
			size_t ret = 0;
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (auto &p: cur_leaf->records){
				if (point_in_mbr(p, query_mbr)){
					ret += 1;
				}
			}
			return ret;
		}
		else{
			auto cur_inte = static_cast<InteNode*>(x.get());
			FT splitter = 1.0 * (1u << (b - 1));
			size_t ret_L = 0, ret_R = 0;
			if (cur_inte->l_son != nullptr) {
				auto L_box = cur_mbr;
				if (x_splitter){
					L_box.second.x = min(x_prefix + splitter - FT_EPS, L_box.second.x);
					ret_L = range_count_node(cur_inte->l_son, query_mbr, L_box, x_prefix, y_prefix, b, !x_splitter);
				}
				else{
					L_box.second.y = min(y_prefix + splitter - FT_EPS, L_box.second.y);
					ret_L = range_count_node(cur_inte->l_son, query_mbr, L_box, x_prefix, y_prefix, b - 1, !x_splitter);
				}
			}
			if (cur_inte->r_son != nullptr) {
				auto R_box = cur_mbr;
				if (x_splitter){
					R_box.first.x = max(x_prefix + splitter, R_box.first.x);
					ret_R = range_count_node(cur_inte->r_son, query_mbr, R_box, x_prefix + splitter, y_prefix, b, !x_splitter);
				}
				else{
					R_box.first.y = max(y_prefix + splitter, R_box.first.y);
					ret_R = range_count_node(cur_inte->r_son, query_mbr, R_box, x_prefix, y_prefix + splitter, b - 1, !x_splitter);
				}
			}
			return ret_L + ret_R;
		}
		return -1; // unexpected error happens if the code runs to here.
	}

	size_t Tree::range_count(Bounding_Box &query_mbr, Bounding_Box &cur_mbr){
		// size_t ret = range_count_node(root, query_mbr, cur_mbr, true);		
		size_t ret = range_count_node(root, query_mbr, cur_mbr, 0.0, 0.0, 32, true);		
		return ret;
	}

	#ifdef BR_MBR
	template <class Out> 
	void Tree::range_report_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, size_t &cnt, Out &out){
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
		if (cur_inte->l_son != nullptr){
			if (mbr_mbr_relation(cur_inte->l_son->mbr, query_mbr) >= 0) range_report_node(cur_inte->l_son, query_mbr, cnt, out);
		}
		if (cur_inte->r_son != nullptr){
 			if (mbr_mbr_relation(cur_inte->r_son->mbr, query_mbr) >= 0) range_report_node(cur_inte->r_son, query_mbr, cnt, out);
		}
	}

	template <class Out>
	void Tree::range_report(Bounding_Box &query_mbr, size_t &cnt, Out &out){
		range_report_node(root, query_mbr, cnt, out);
	}
	#endif

	template <class Out> 
	void Tree::range_report_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter, size_t &cnt, Out &out){
		// print_mbr(cur_mbr);
		// print_mbr(query_mbr);
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
		FT splitter = 1.0 * (1u << (b - 1));
		if (cur_inte->l_son != nullptr){
			auto L_box = cur_mbr;
			if (x_splitter){
				L_box.second.x = min(x_prefix + splitter - FT_EPS, L_box.second.x);
				if (mbr_mbr_relation(L_box, query_mbr) >= 0) range_report_node(cur_inte->l_son, query_mbr, L_box, x_prefix, y_prefix, b, !x_splitter, cnt, out);
			}
			else{
				L_box.second.y = min(y_prefix + splitter - FT_EPS, L_box.second.y);
				if (mbr_mbr_relation(L_box, query_mbr) >= 0) range_report_node(cur_inte->l_son, query_mbr, L_box, x_prefix, y_prefix, b - 1, !x_splitter, cnt, out);
			}
		}
		if (cur_inte->r_son != nullptr){
			auto R_box = cur_mbr;
			if (x_splitter){
				R_box.first.x = max(x_prefix + splitter, R_box.first.x);
 				if (mbr_mbr_relation(R_box, query_mbr) >= 0) range_report_node(cur_inte->r_son, query_mbr, R_box, x_prefix + splitter, y_prefix, b, !x_splitter, cnt, out);
			}
			else{
				R_box.first.y = max(y_prefix + splitter, R_box.first.y);
 				if (mbr_mbr_relation(R_box, query_mbr) >= 0) range_report_node(cur_inte->r_son, query_mbr, R_box, x_prefix, y_prefix + splitter, b - 1, !x_splitter, cnt, out);
			}
		}
	}

	template <class Out>
	void Tree::range_report(Bounding_Box &query_mbr, Bounding_Box &cur_mbr, size_t &cnt, Out &out){
		range_report_node(root, query_mbr, cur_mbr, 0.0, 0.0, 32, true, cnt, out);
	}

	#ifdef BR_MBR
	template <class T>
	void Tree::knn_report_node(shared_ptr<BaseNode> &x, size_t &k, Point query_point, T &nn_res){
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (auto &p: cur_leaf->records){
				auto cur_sqrdis = point_point_sqrdis(p, query_point);
				if (nn_res.size() < k){
					nn_res.push({p, cur_sqrdis});
				}
				else if (cur_sqrdis < nn_res.top().second){
					nn_res.pop();
					nn_res.push({p, cur_sqrdis});
				}
			}
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		auto l_son_sqrdis = FT_INF_MAX, r_son_sqrdis = FT_INF_MAX;
		if (cur_inte->l_son != nullptr) l_son_sqrdis = point_mbr_sqrdis(query_point, cur_inte->l_son->mbr);
		if (cur_inte->r_son != nullptr) r_son_sqrdis = point_mbr_sqrdis(query_point, cur_inte->r_son->mbr);
		if (l_son_sqrdis <= r_son_sqrdis){ // first go left
			if (nn_res.size() < k || l_son_sqrdis < nn_res.top().second){
				knn_report_node(cur_inte->l_son, k, query_point, nn_res);
			}
			if (nn_res.size() < k || r_son_sqrdis < nn_res.top().second){
				knn_report_node(cur_inte->r_son, k, query_point, nn_res);
			}
		}
		else{	// first go right
			if (nn_res.size() < k || r_son_sqrdis < nn_res.top().second){
				knn_report_node(cur_inte->r_son, k, query_point, nn_res);
			}
			if (nn_res.size() < k || l_son_sqrdis < nn_res.top().second){
				knn_report_node(cur_inte->l_son, k, query_point, nn_res);
			}
		}
		return;
	}

	auto Tree::knn_report(size_t &k, Point query_point){
		priority_queue<nn_pair, vector<nn_pair>, nn_pair_cmp> nn_res;
		knn_report_node(root, k, query_point, nn_res);
		return nn_res;
	}
	#endif

	auto Tree::knn_report(size_t &k, Point query_point, Bounding_Box &cur_mbr){
		priority_queue<nn_pair, vector<nn_pair>, nn_pair_cmp> nn_res;
		knn_report_node(root, k, query_point, cur_mbr, 0.0, 0.0, 32, true, nn_res);
		return nn_res;
	}

	template <class T>
	void Tree::knn_report_node(shared_ptr<BaseNode> &x, size_t &k, Point query_point, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter, T &nn_res){
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (auto &p: cur_leaf->records){
				auto cur_sqrdis = point_point_sqrdis(p, query_point);
				if (nn_res.size() < k){
					nn_res.push({p, cur_sqrdis});
				}
				else if (cur_sqrdis < nn_res.top().second){
					nn_res.pop();
					nn_res.push({p, cur_sqrdis});
				}
			}
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		auto l_son_sqrdis = FT_INF_MAX, r_son_sqrdis = FT_INF_MAX;
		FT splitter = 1.0 * (1u << (b - 1));
		auto L_box = cur_mbr, R_box = cur_mbr;
		if (cur_inte->l_son != nullptr) {
			if (x_splitter){
				L_box.second.x = min(x_prefix + splitter - FT_EPS, L_box.second.x);
			}
			else{
				L_box.second.y = min(y_prefix + splitter - FT_EPS, L_box.second.y);
			}
			l_son_sqrdis = point_mbr_sqrdis(query_point, L_box);
		}
		if (cur_inte->r_son != nullptr) {
			if (x_splitter){
				R_box.first.x = max(x_prefix + splitter, R_box.first.x);
			}
			else{
				R_box.first.y = max(y_prefix + splitter, R_box.first.y);
			}
			r_son_sqrdis = point_mbr_sqrdis(query_point, R_box);
		}
		// auto go_left = [&](){
		// 	if (nn_res.size() < k || l_son_sqrdis < nn_res.top().second){
		// 		if (x_splitter) knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix, y_prefix, b, !x_splitter, nn_res);
		// 		else knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix, y_prefix, b - 1, !x_splitter, nn_res);
		// 	}
		// };
		// auto go_right = [&](){
		// 	if (nn_res.size() < k || r_son_sqrdis < nn_res.top().second){
		// 		if (x_splitter) knn_report_node(cur_inte->r_son, k, query_point, R_box, x_prefix + splitter, y_prefix, b, !x_splitter, nn_res);
		// 		else knn_report_node(cur_inte->r_son, k, query_point, R_box, x_prefix, y_prefix + splitter, b - 1, !x_splitter, nn_res);
		// 	}
		// };
		if (l_son_sqrdis <= r_son_sqrdis){ // first go left
			// go_left();
			// go_right();
			if (nn_res.size() < k || l_son_sqrdis < nn_res.top().second){
				if (x_splitter) knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix, y_prefix, b, !x_splitter, nn_res);
				else knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix, y_prefix, b - 1, !x_splitter, nn_res);
			}
			if (nn_res.size() < k || r_son_sqrdis < nn_res.top().second){
				if (x_splitter) knn_report_node(cur_inte->r_son, k, query_point, R_box, x_prefix + splitter, y_prefix, b, !x_splitter, nn_res);
				else knn_report_node(cur_inte->r_son, k, query_point, R_box, x_prefix, y_prefix + splitter, b - 1, !x_splitter, nn_res);
			}
		}
		else{	// first go right
			// go_right();
			// go_left();
			if (nn_res.size() < k || r_son_sqrdis < nn_res.top().second){
				if (x_splitter) knn_report_node(cur_inte->r_son, k, query_point, R_box, x_prefix + splitter, y_prefix, b, !x_splitter, nn_res);
				else knn_report_node(cur_inte->r_son, k, query_point, R_box, x_prefix, y_prefix + splitter, b - 1, !x_splitter, nn_res);
				
			}
			if (nn_res.size() < k || l_son_sqrdis < nn_res.top().second){
				if (x_splitter) knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix, y_prefix, b, !x_splitter, nn_res);
				else knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix, y_prefix, b - 1, !x_splitter, nn_res);
			}
		}
		return;
	}

}