#ifndef PROFILER_H
#define PROFILER_H

#include <chrono>
#include <string>
#include <unordered_map>
#include <iostream>
#include <iomanip>

#ifdef PROFILING_ENABLED
#define PROFILE_FUNCTION() ProfileTimer timer(__FUNCTION__)
#define PROFILE_SCOPE(name) ProfileTimer timer(name)
#else
#define PROFILE_FUNCTION()
#define PROFILE_SCOPE(name)
#endif

class Profiler {
private:
    std::unordered_map<std::string, long long> total_times;
    std::unordered_map<std::string, int> call_counts;
    
    Profiler() = default;
    
public:
    static Profiler& getInstance() {
        static Profiler instance;
        return instance;
    }
    
    void addTime(const std::string& name, long long microseconds) {
        total_times[name] += microseconds;
        call_counts[name]++;
    }
    
    void printReport() {
        std::cout << "\n=== PROFILING REPORT ===" << std::endl;
        std::cout << "Function Name                    | Total Time (ms) | Calls | Avg Time (μs)" << std::endl;
        std::cout << "--------------------------------|-----------------|-------|---------------" << std::endl;
        
        for (const auto& pair : total_times) {
            const std::string& name = pair.first;
            long long total_us = pair.second;
            int calls = call_counts[name];
            double avg_us = static_cast<double>(total_us) / calls;
            
            std::cout << std::left << std::setw(32) << name 
                      << "| " << std::setw(15) << (total_us / 1000.0)
                      << "| " << std::setw(5) << calls
                      << "| " << std::setw(13) << avg_us << std::endl;
        }
        std::cout << "=========================" << std::endl;
    }
    
    void reset() {
        total_times.clear();
        call_counts.clear();
    }
};

class ProfileTimer {
private:
    std::string name;
    std::chrono::high_resolution_clock::time_point start_time;
    
public:
    ProfileTimer(const std::string& function_name) : name(function_name) {
        start_time = std::chrono::high_resolution_clock::now();
    }
    
    ~ProfileTimer() {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        Profiler::getInstance().addTime(name, duration.count());
    }
};

#endif
