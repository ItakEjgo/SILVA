#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <silva/geo/point.hpp>
#include <silva/geo/operations.hpp>

struct Operation {
    char op_type;
    size_t node_id;
    double lon;
    double lat;
};

struct Event {
    int type; // 0 = START, 1 = END
    std::string timestamp;
    size_t changeset_id;
    std::vector<Operation> ops;
    int year; // keep track of the year this event belongs to
};

class DatasetLoader {
public:
    std::vector<geobase::Point> base_data;
    std::vector<Event> all_events;
    geobase::Bounding_Box largest_mbr;

    void load(const std::string& dataset_dir, const std::string& commits_dir, int start_year, int end_year) {
        std::cout << "[Init] Pre-scanning for Dynamic MBR & Loading memory workload..." << std::endl;
        double min_x = 1e18, min_y = 1e18, max_x = -1e18, max_y = -1e18;
        
        auto update_mbr = [&](double lon, double lat) {
            if (lon < min_x) min_x = lon;
            if (lat < min_y) min_y = lat;
            if (lon > max_x) max_x = lon;
            if (lat > max_y) max_y = lat;
        };

        // 1. Load Base Snapshot
        std::ifstream fin_base(dataset_dir + "/00_build/base_snapshot_" + std::to_string(start_year - 1) + ".csv");
        if (!fin_base.is_open()) {
            std::cerr << "Failed to open base snapshot!" << std::endl;
            return;
        }
        std::string line;
        std::getline(fin_base, line);
        while (std::getline(fin_base, line)) {
            std::stringstream ss(line);
            std::string nid_str, ver_str, lon_str, lat_str;
            std::getline(ss, nid_str, ','); std::getline(ss, ver_str, ',');
            std::getline(ss, lon_str, ','); std::getline(ss, lat_str, ',');
            double lon = std::stod(lon_str) * 1000000.0;
            double lat = std::stod(lat_str) * 1000000.0;
            update_mbr(lon, lat);
            base_data.push_back(geobase::Point(std::stoull(nid_str), lon, lat));
        }
        fin_base.close();

        // 2. Load Commits into Memory
        for (int year = start_year; year <= end_year; year++) {
            std::ifstream fin(dataset_dir + "/" + commits_dir + "/commits_" + std::to_string(year) + ".csv");
            if (!fin.is_open()) continue;
            std::getline(fin, line);
            
            Event current_event;
            current_event.changeset_id = 0;
            current_event.year = year;
            bool first = true;
            
            while(std::getline(fin, line)) {
                std::stringstream ss(line);
                std::string cs_str, start_str, end_str, op_str, nid_str, ver_str, lon_str, lat_str;
                std::getline(ss, cs_str, ','); std::getline(ss, start_str, ','); std::getline(ss, end_str, ',');
                std::getline(ss, op_str, ','); std::getline(ss, nid_str, ','); std::getline(ss, ver_str, ',');
                std::getline(ss, lon_str, ','); std::getline(ss, lat_str, ',');
                
                size_t cs = std::stoull(cs_str);
                Operation op;
                op.op_type = op_str[0];
                op.node_id = std::stoull(nid_str);
                if (op.op_type != 'D') {
                    op.lon = std::stod(lon_str) * 1000000.0;
                    op.lat = std::stod(lat_str) * 1000000.0;
                    update_mbr(op.lon, op.lat);
                } else {
                    op.lon = 0; op.lat = 0;
                }
                
                if (first) {
                    current_event.changeset_id = cs;
                    first = false;
                }
                
                if (current_event.changeset_id != cs) {
                    all_events.push_back(current_event);
                    current_event.ops.clear();
                    current_event.changeset_id = cs;
                }
                current_event.ops.push_back(op);
            }
            if (!current_event.ops.empty()) {
                all_events.push_back(current_event);
            }
        }

        // 1% Padding for MBR
        double dx = (max_x - min_x) * 0.01;
        double dy = (max_y - min_y) * 0.01;
        largest_mbr = geobase::Bounding_Box(
            geobase::Point(min_x - dx, min_y - dy), 
            geobase::Point(max_x + dx, max_y + dy)
        );
        std::cout << "[Init] Load complete. MBR: X[" << min_x-dx << ", " << max_x+dx << "] Y[" << min_y-dy << ", " << max_y+dy << "]" << std::endl;
        std::cout << "[Init] Base data size: " << base_data.size() << " points." << std::endl;
        std::cout << "[Init] Total changesets loaded: " << all_events.size() << std::endl;
    }
};
