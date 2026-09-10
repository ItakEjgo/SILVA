#include <parlay/primitives.h>
#include <parlay/sequence.h>
#pragma once
#include "point.hpp"
namespace geobase {
    template <typename PT>
    auto filter_diff_results(PT &add, PT &remove)
    {
        auto id_cmp = [&](auto lhs, auto rhs)
        { return lhs.id < rhs.id; };
        auto sorted_add = parlay::sort(add, id_cmp);
        auto sorted_remove = parlay::sort(remove, id_cmp);

        PT insert_points = {}, delete_points = {}, update_points = {};
        size_t i = 0, j = 0;
        while (i < sorted_add.size() && j < sorted_remove.size())
        {
            if (sorted_add[i].id < sorted_remove[j].id)
            { //	we should insert sorted_add[i]
                insert_points.emplace_back(sorted_add[i++]);
            }
            else if (sorted_add[i].id == sorted_remove[j].id)
            { //	exist in both add and remove, indicate it is an update point
                update_points.emplace_back(sorted_add[i++]);
                j++;
            }
            else
            { //	we should remove sorted_remove[j]
                delete_points.emplace_back(sorted_remove[j++]);
            }
        }
        while (i < sorted_add.size())
            insert_points.emplace_back(sorted_add[i++]);
        while (j < sorted_remove.size())
            delete_points.emplace_back(sorted_remove[j++]);
        return make_tuple(insert_points, delete_points, update_points);
    }

    template <class MBR>
    bool is_same_mbr(MBR &lhs, MBR &rhs)
    {
        // cout << "[(" << mbr.first.x << ", " << mbr.first.y << "), (" << mbr.second.x << ", " << mbr.second.y << ")" << "]" << endl;
        return lhs.first.x == rhs.first.x && lhs.first.y == rhs.first.y && lhs.second.x == rhs.second.x && lhs.second.y == rhs.second.y;
    }

    template <class T>
    size_t split_by_bit(T &P, size_t l, size_t r, size_t b){
        unsigned int splitter = (1u << ((b - 1) / 2));
        auto less = [&](auto pt){
            return b % 2 ? (static_cast<unsigned int>(pt.y) & splitter) == 0 : (static_cast<unsigned int>(pt.x) & splitter) == 0;
        };
        size_t start = l, end = r, mid = start;
        while (end - start > 16){
            mid = (start + end) / 2;
            if (!less(P[mid]))
                end = mid;
            else
                start = mid + 1;
        }
        mid = end;
        for (auto i = start; i < end; i++){
            if (!less(P[i])){
                mid = i;
                return i;
            }
        }
        return mid;
    }

    // merge two bounding boxes, please make sure *both bounding boxes are valid*.
    template <class MBR>
    MBR merge_mbr(MBR &a, MBR &b)
    {
        return {Point(min(a.first.x, b.first.x), min(a.first.y, b.first.y)), Point(max(a.second.x, b.second.x), max(a.second.y, b.second.y))};
    }

    // check whether a given point inside a mbr (boundary included).
    template <class MBR>
    bool point_in_mbr(Point &p, MBR &mbr)
    {
        return mbr.first.x <= p.x && p.x <= mbr.second.x &&
               mbr.first.y <= p.y && p.y <= mbr.second.y;
    }

    // check whether a given mbr inside a mbr, e,g., check whether small_mbr in large_mbr (boundary included).
    template <class MBR>
    bool mbr_in_mbr(MBR &small_mbr, MBR &large_mbr)
    {
        return small_mbr.first.x >= large_mbr.first.x && small_mbr.second.x <= large_mbr.second.x &&
               small_mbr.first.y >= large_mbr.first.y && small_mbr.second.y <= large_mbr.second.y;
    }

    // check whether two given mbrs exclude each other (do not check boundaries).
    template <class MBR>
    bool mbr_exclude_mbr(MBR &small_mbr, MBR &large_mbr)
    {
        return (small_mbr.first.x > large_mbr.second.x || small_mbr.second.x < large_mbr.first.x) ||
               (small_mbr.first.y > large_mbr.second.y || small_mbr.second.y < large_mbr.first.y);
    }

    // relations between two mbrs: -1: excluded; 0: intersected; 1: the smaller one is contained by the larger one;
    template <class MBR>
    int mbr_mbr_relation(MBR &small_mbr, MBR &large_mbr)
    {
        auto minc_x = max(small_mbr.first.x, large_mbr.first.x);
        auto minc_y = max(small_mbr.first.y, large_mbr.first.y);
        auto maxc_x = min(small_mbr.second.x, large_mbr.second.x);
        auto maxc_y = min(small_mbr.second.y, large_mbr.second.y);
        if (minc_x <= maxc_x && minc_y <= maxc_y)
        {
            if (minc_x == small_mbr.first.x && maxc_x == small_mbr.second.x &&
                minc_y == small_mbr.first.y && maxc_y == small_mbr.second.y)
                return 1;
            return 0;
        }
        return -1;
    }

    template <class Records>
    auto get_mbr(Records &P)
    {
        FT x_min = FT_INF_MAX, x_max = FT_INF_MIN, y_min = FT_INF_MAX, y_max = FT_INF_MIN;
        for (auto &p : P)
        {
            x_min = min(x_min, p.x);
            x_max = max(x_max, p.x);
            y_min = min(y_min, p.y);
            y_max = max(y_max, p.y);
        }
        return Bounding_Box({Point(x_min, y_min), Point(x_max, y_max)});
    }

    auto mbr_mbr_within_dis(Bounding_Box &mbr1, Bounding_Box &mbr2, FT &point_dis){
        return !(dcmp(mbr1.second.x + point_dis - mbr2.first.x) < 0 ||
                dcmp(mbr1.second.y + point_dis - mbr2.first.y) < 0 ||
                dcmp(mbr2.second.x + point_dis - mbr1.first.x) < 0 ||
                dcmp(mbr2.second.y + point_dis - mbr1.first.y) < 0);
    }

    // return the sqr distance between a point and a mbr
    template <class MBR>
    auto point_mbr_sqrdis(Point &p, MBR &mbr)
    {
        FT dx = max(max(mbr.first.x - p.x, (FT)0.0), max(p.x - mbr.second.x, (FT)0.0));
        FT dy = max(max(mbr.first.y - p.y, (FT)0.0), max(p.y - mbr.second.y, (FT)0.0));
        return dx * dx + dy * dy;
    }

    // return the sqr distance between two points
    auto point_point_sqrdis(Point &lhs, Point &rhs)
    {
        return (lhs.x - rhs.x) * (lhs.x - rhs.x) + (lhs.y - rhs.y) * (lhs.y - rhs.y);
    }

    // check whether two given MBRs are intersected.
    template <class MBR>
    bool mbr_intersects_mbr(MBR &a, MBR &b)
    {
        if (a.first.x > b.second.x || a.second.x < b.first.x)
            return false;
        if (a.first.y > b.second.y || a.second.y < b.first.y)
            return false;
        return true;
    }

typedef pair<Point, FT> nn_pair;

    struct nn_pair_cmp
    {
        bool operator()(const nn_pair &lhs, const nn_pair &rhs) const
        {
            return lhs.second < rhs.second ||
                   (lhs.second == rhs.second && lhs.first.id > rhs.first.id);
        }
    };

    template <class Pset>
    FT knn_bf(size_t &k, Point &q, Pset &P)
    {
        vector<FT> q_sqrdis = {};
        for (size_t i = 0; i < P.size(); i++)
        {
            auto cur_sqrdis = point_point_sqrdis(q, P[i]);
            q_sqrdis.emplace_back(cur_sqrdis);
        }
        sort(q_sqrdis.begin(), q_sqrdis.end());
        return q_sqrdis[k - 1];
    }

    template <class Pset>
    void morton_sort(Pset &P)
    {
        sort(P.begin(), P.end(), [&](auto lhs, auto rhs)
             {
		    auto msd = 0;
		    if (geobase::less_msb(static_cast<unsigned int>(lhs.x) ^ static_cast<unsigned int>(rhs.x), static_cast<unsigned int>(lhs.y) ^ static_cast<unsigned int>(rhs.y))) 
			    msd = 1;
		    return !msd ? lhs.x < rhs.x : lhs.y < rhs.y; });
    }

    template <class P>
    bool morton_less(const P &lhs, const P &rhs)
    {
        auto msd = 0;
        if (geobase::less_msb(static_cast<unsigned int>(lhs->x) ^ static_cast<unsigned int>(rhs->x), static_cast<unsigned int>(lhs->y) ^ static_cast<unsigned int>(rhs->y)))
            msd = 1;
        return !msd ? lhs->x < rhs->x : lhs->y < rhs->y;
    }

    template <typename Pset>
    auto morton_merge(Pset &lhs, Pset &rhs)
    {
        size_t i = 0, j = 0, k = 0, n = lhs.size(), m = rhs.size();
        Pset ret;
        ret.resize(n + m);
        while (i < n && j < m)
        {
            if (morton_less(lhs[i], rhs[j]))
            {
                ret[k++] = lhs[i++];
            }
            else
            {
                ret[k++] = rhs[j++];
            }
        }
        while (i < n)
            ret[k++] = lhs[i++];
        while (j < m)
            ret[k++] = rhs[j++];
        return ret;
    }

    template <typename Pset>
    auto shuffle_point(Pset &P, size_t substr_size = 0){
        auto n = P.size();
        parlay::sequence<Point> rand_p(n);
        // random_device rd;
        auto rand_seed = 233666;
        mt19937 g(rand_seed);
        parlay::parallel_for(0, n, [&](size_t i)
                             { rand_p[i] = P[i]; });
        shuffle(rand_p.begin(), rand_p.end(), g);
        if (substr_size > 0){
            rand_p = rand_p.substr(0, substr_size);
        }
        return rand_p;
    }

    /* Insert ratio is from 0 to 10. Insert point will add bias to id to make them unique. */
    template <typename Pset>
    auto split_insert_delete(Pset &P, size_t &insert_ratio, size_t id_bias)
    {
        size_t insert_num = P.size() / 10 * insert_ratio;
        size_t delete_num = P.size() - insert_num;
        auto P_insert = P.substr(0, insert_num);
        parlay::parallel_for(0, insert_num, [&](size_t i)
                             { P_insert[i].id += id_bias; });
        auto P_delete = P.substr(P.size() - delete_num, delete_num);
        return make_tuple(P_insert, P_delete);
    }

    template <typename Pset>
    auto print_Pset_info(Pset &P, string name = "debug", size_t end_idx = 0){
        cout << name << ": " << endl;
        if (end_idx == 0) end_idx = P.size();
        for (size_t i = 0; i < end_idx; i++){
            cout << P[i] << endl;
        }
    }

    template <typename Pset>
    auto collect_newver_point(Pset &P, Pset &P_insert, Pset &P_delete){
        unordered_map<size_t, bool> removed = {};
        for (auto &pt : P_delete){
            removed[pt.id] = true;
        }
        auto ret = parlay::sequence<Point>::uninitialized(P.size() + P_insert.size());
        auto cnt = 0;
        for (auto &pt : P){
            if (!removed[pt.id])
                parlay::assign_uninitialized(ret[cnt++], pt);
        }
        for (auto &pt : P_insert){
            // if (!removed[pt.id])
                parlay::assign_uninitialized(ret[cnt++], pt);
        }
        ret.resize(cnt);
        return ret;
    }

    auto compute_cur_box(Bounding_Box &cur_mbr, FT &x_prefix, FT &y_prefix, size_t b, bool x_splitter){
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
        return make_tuple(L_box, R_box, rx_prefix, ry_prefix);
    }

    template <typename Func>
    auto get_delete_p(parlay::sequence<Point> &lhs, parlay::sequence<Point> &P, size_t l, size_t r, Func &f)
    {
        auto ret = parlay::sequence<Point>::uninitialized(lhs.size());
        auto pt_cmp = [&](auto &pt1, auto &pt2){
            if (pt1.morton_id == pt2.morton_id){
                if (pt1.id == pt2.id){
                    return 0;
                }
                else{
                    if (pt1.id < pt2.id)
                        return 1;
                    else
                        return 2;
                }
            }
            else {
                if (pt1.morton_id < pt2.morton_id)
                    return 1;
                else
                    return 2;
            }
        };
        size_t i = 0, j = l, cnt = 0;
        while (i < lhs.size() && j < r){
            auto flag = pt_cmp(lhs[i], P[j]);
            if (!flag){ // point should be deleted
                i++, j++;
            }
            else if (flag == 1){ // first smaller, in A not in B
                if (!f(lhs[i])){
                    i++;
                    continue;
                }
                parlay::assign_uninitialized(ret[cnt++], lhs[i++]);
            }
            else{ //  second smaller, in B not in A
                j++;
            }
        }
        while (i < lhs.size()){
            if (!f(lhs[i])){
                i++;
                continue;
            }
            parlay::assign_uninitialized(ret[cnt++], lhs[i++]);
        }
        ret.resize(cnt);
        return ret;
    }

    auto get_delete_p(parlay::sequence<Point> &lhs, parlay::sequence<Point> &P, size_t l, size_t r)
    {
        auto f_true = [](const auto &) { return true; };
        return get_delete_p(lhs, P, l, r, f_true);
    }

    // delete rhs from lhs, merge by morton_id and point id, apply F
    template <typename T, typename P, typename F, typename DIFF>
    auto merge_pts(P &a, T &b, F &f, DIFF &ret_diff, bool reverse = false){
        size_t i = 0, j = 0;
        //  0: same point, 1: first is smaller, 2: second is smaller
        auto pt_cmp = [&](const auto &pt1, const auto &pt2){
            if (pt1.morton_id == pt2.morton_id){
                if (pt1.id == pt2.id){
                    return 0;
                }
                else{
                    return pt1.id < pt2.id ? 1 : 2;
                }
            }
            else{
                return pt1.morton_id < pt2.morton_id ? 1 : 2;
            }
        };

        while (i < a.size() && j < b.size()){
            auto flag = pt_cmp(a[i], b[j]);

            if (!flag){ // same point
                i++, j++;
            }
            else if (flag == 1){ // first smaller, in A not in B
                if (f(a[i])) {
                    ret_diff.remove_point(a[i], reverse);
                }

                i++;
            }
            else{ //  second smaller, in B not in A
                if (f(b[j])) {
                    ret_diff.add_point(b[j], reverse);
                }
                j++;
            }
        }
        while (i < a.size()){
            if (f(a[i])) {
                ret_diff.remove_point(a[i], reverse);
            }
            i++;
        }
        while (j < b.size()){
            if (f(b[j])){
                ret_diff.add_point(b[j], reverse);
            } 
            j++;
        }
        
        return;
    }

    // delete rhs from lhs, merge by morton_id and point id
    template <typename T, typename P, typename DIFF>
    auto merge_pts(P &a, T &b, DIFF &ret_diff, bool reverse = false){
        auto f_true = [](const auto &) { return true; };
        return merge_pts(a, b, f_true, ret_diff, reverse);
    }

    // delete rhs from lhs, merge by morton_id and point id
    template <typename T, typename P>
    auto merge_pts(P &a, T &b){
        diff_type ret_diff(b.size(), a.size());
        auto f_true = [](const auto &) { return true; };
        merge_pts(a, b, f_true, ret_diff, false);
        ret_diff.compact();
        return make_tuple(ret_diff.add, ret_diff.remove);
    }
    

    // generate n random points within the largest bounding box
    template <typename Pset>
    auto random_sample(size_t n, Pset P)
    {
        unsigned seed = 233666;
        default_random_engine e(seed);
        vector<size_t> sampled_ids(P.size());
        parlay::parallel_for(0, P.size(), [&](size_t i)
                             { sampled_ids[i] = i; });
        shuffle(sampled_ids.begin(), sampled_ids.end(), e);
        auto P2 = P.substr(0, n);
        for (size_t i = 0; i < n; i++)
        {
            P2[i] = P[sampled_ids[i]];
        }
        return P2;
    }

    //  return a set of sorted points
    template <typename PT>
    auto get_sorted_points(PT &P, bool use_hilbert = false)
    {
        auto n = P.size();
        
        if (n < 256) {
            for(size_t i = 0; i < n; i++) {
                P[i].morton_id = use_hilbert ? P[i].overlap_bits() : P[i].interleave_bits();
            }
        } else {
            parlay::parallel_for(0, n, [&](int i)
                                 { P[i].morton_id = use_hilbert ? P[i].overlap_bits() : P[i].interleave_bits(); });
        }

        auto morton_cmp = [&](auto lhs, auto rhs) {
            return lhs.morton_id < rhs.morton_id ||
                   (lhs.morton_id == rhs.morton_id && lhs.id < rhs.id); 
        };

        if (n < 256) {
            PT P_copy = P;
            std::sort(P_copy.begin(), P_copy.end(), morton_cmp);
            return P_copy;
        } else {
            return parlay::sort(P, morton_cmp);
        }
        // if (use_hilbert){
        // }
        // else{
        //     return parlay::sort(P,  [&](auto lhs, auto rhs){
        // 	    auto msd = 0;
        // 	    if (geobase::less_msb(static_cast<unsigned int>(lhs.x) ^ static_cast<unsigned int>(rhs.x), static_cast<unsigned int>(lhs.y) ^ static_cast<unsigned int>(rhs.y)))
        // 		    msd = 1;
        // 	    return !msd ? lhs.x < rhs.x : lhs.y < rhs.y;
        //     });
        // }
    }

    template <typename PT>
    auto get_sorted_points_hilbert(PT &P, bool use_hilbert = true)
    {
        auto n = P.size();
        
        if (n < 256) {
            for(size_t i = 0; i < n; i++) {
                P[i].morton_id = use_hilbert ? P[i].overlap_bits() : P[i].interleave_bits();
            }
        } else {
            parlay::parallel_for(0, n, [&](int i)
                                 { P[i].morton_id = use_hilbert ? P[i].overlap_bits() : P[i].interleave_bits(); });
        }

        auto morton_cmp = [&](auto lhs, auto rhs) {
            return lhs.morton_id < rhs.morton_id ||
                   (lhs.morton_id == rhs.morton_id && lhs.id < rhs.id); 
        };

        if (n < 256) {
            PT P_copy = P;
            std::sort(P_copy.begin(), P_copy.end(), morton_cmp);
            return P_copy;
        } else {
            return parlay::sort(P, morton_cmp);
        }
    }

    //	return a set of sorted address, do not modify the input;
    template <typename PT>
    auto get_sorted_address(PT &P)
    {
        parlay::sequence<geobase::Point *> P_set(P.size());
        parlay::parallel_for(0, P.size(), [&](int i)
                             { P_set[i] = &P[i]; });
        P_set = parlay::sort(P_set, [&](auto lhs, auto rhs)
                             {
			auto msd = 0;
			if (geobase::less_msb(static_cast<unsigned int>(lhs->x) ^ static_cast<unsigned int>(rhs->x), static_cast<unsigned int>(lhs->y) ^ static_cast<unsigned int>(rhs->y))) 
				msd = 1;
			return !msd ? lhs->x < rhs->x : lhs->y < rhs->y; });
        return P_set;
    }

    template <typename PT>
    auto get_unsorted_address(PT &P)
    {
        parlay::sequence<Point *> P_set(P.size());
        parlay::parallel_for(0, P.size(), [&](int i)
                             { P_set[i] = &P[i]; });
        return P_set;
    }

    template <typename PT>
    auto record_check(PT &P)
    {
        map<int, int> mmp = {};
        auto cnt = 0;
        auto num = 0;
        for (size_t i = 0; i < P.size(); i++)
        {
            mmp[P[i]->id]++;
        }
        for (auto &key_val : mmp)
        {
            if (key_val.second > 1)
            {
                num++;
                cnt += key_val.second;
                // cout << key_val.first << ", " << key_val.second << endl;
            }
        }
        return make_tuple(num, cnt);
    }

    // return a set contains all ids in a sequence of point pointer
    template <typename PT>
    auto get_point_id(PT &P)
    {
        set<int> ret = {};
        for (auto pt : P)
        {
            ret.insert(pt->id);
        }
        return ret;
    }

    template <typename PT>
    auto count_duplicate_zvalues(PT &P)
    {
        unordered_map<unsigned long long, int> zvalue_map = {};
        for (auto &pt : P)
        {
            zvalue_map[pt.morton_id]++;
        }
        auto cnt = 0, total_cnt = 0;
        for (const auto &pair : zvalue_map)
        {
            if (pair.second > 1)
            {
                // cout << "Number " << pair.first << " appears " << pair.second << " times." << endl;
                cnt++;
                total_cnt += pair.second;
            }
        }
        cout << "Total # of duplicate keys: " << cnt << ", " << total_cnt << endl;
        return cnt;
    }

    template <typename T>
    auto diff_by_id(T &a, T &b){
        auto id_cmp = [&](auto &lhs, auto &rhs) { return lhs.id < rhs.id; };
        
        auto sorted_a = parlay::sort(a, id_cmp);
        auto sorted_b = parlay::sort(b, id_cmp);

        auto added = T::uninitialized(b.size());
        auto removed = T::uninitialized(a.size());

        size_t cnt_add = 0, cnt_remove = 0;
        size_t i = 0, j = 0;

        while (i < a.size() && j < b.size()){
            // cout << i << ", " << j << endl;
            if (sorted_a[i].id == sorted_b[j].id){ // same point, check coordinate
                if (sorted_a[i] == sorted_b[j]){ //  same coordinate, do nothing
                    i++, j++;
                }
                else{
                    parlay::assign_uninitialized(removed[cnt_remove++], sorted_a[i++]);
                    parlay::assign_uninitialized(added[cnt_add++], sorted_b[j++]);
                }
            }
            else if (sorted_a[i].id < sorted_b[j].id){ // first smaller, in A not in B
                parlay::assign_uninitialized(removed[cnt_remove++], sorted_a[i++]);
            }
            else{ //  second smaller, in B not in A
                parlay::assign_uninitialized(added[cnt_add++], sorted_b[j++]);
            }
        }
        while (i < sorted_a.size()){
            parlay::assign_uninitialized(removed[cnt_remove++], sorted_a[i++]);
        }
        while (j < sorted_b.size()){
            parlay::assign_uninitialized(added[cnt_add++], sorted_b[j++]);
        }
        removed.resize(cnt_remove);
        added.resize(cnt_add);

        return make_tuple(added, removed);
    }

    //  merge to sequence by id, return two sequences. The first one has no conflict, the second one has conflict
    template <typename PT>
    auto merge_by_id_with_conflict(PT &a, PT &b)
    {
        size_t i = 0, j = 0;
        auto id_cmp = [&](auto lhs, auto rhs)
        { return lhs.id < rhs.id; };
        auto sorted_a = parlay::sort(a, id_cmp);
        auto sorted_b = parlay::sort(b, id_cmp);
        PT ret_noconflict = {}, ret_conflict = {};
        while (i < sorted_a.size() && j < sorted_b.size())
        {
            if (sorted_a[i].id == sorted_b[j].id)
            { //  two points have the same id
                if (!(sorted_a[i] == sorted_b[j]))
                { //  coordinate does not match, conflict
                    ret_conflict.emplace_back(sorted_a[i++]);
                    ret_conflict.emplace_back(sorted_b[j++]);
                }
                else
                { //  no conflict
                    ret_noconflict.emplace_back(sorted_a[i++]);
                    j++;
                }
            }
            else
            { //  no conflict, store the one with smaller id
                if (sorted_a[i].id < sorted_b[j].id)
                {
                    ret_noconflict.emplace_back(sorted_a[i++]);
                }
                else
                {
                    ret_noconflict.emplace_back(sorted_b[j++]);
                }
            }
        }
        while (i < sorted_a.size())
            ret_noconflict.emplace_back(sorted_a[i++]);
        while (j < sorted_b.size())
            ret_noconflict.emplace_back(sorted_b[j++]);
        return make_tuple(ret_noconflict, ret_conflict);
    }
}
