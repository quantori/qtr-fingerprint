#pragma once

#include <unordered_map>
#include <string>
#include <stdexcept>

class SearchEngineConfig {
public:
    SearchEngineConfig() = default;
    
    inline void set(const std::string& key, const std::string& value) {
        _config[key] = value;
    }
    
    inline bool has(const std::string& key) const {
        return _config.find(key) != _config.end();
    }
    
    inline std::string getString(const std::string& key, const std::string& defaultValue = "") const {
        auto it = _config.find(key);
        if (it != _config.end()) {
            return it->second;
        }
        return defaultValue;
    }
    
    inline int getInt(const std::string& key, int defaultValue = 0) const {
        auto it = _config.find(key);
        if (it != _config.end()) {
            try {
                return std::stoi(it->second);
            } catch (const std::exception& e) {
                return defaultValue;
            }
        }
        return defaultValue;
    }
    
    inline double getDouble(const std::string& key, double defaultValue = 0.0) const {
        auto it = _config.find(key);
        if (it != _config.end()) {
            try {
                return std::stod(it->second);
            } catch (const std::exception& e) {
                return defaultValue;
            }
        }
        return defaultValue;
    }
    
    inline bool getBool(const std::string& key, bool defaultValue = false) const {
        auto it = _config.find(key);
        if (it != _config.end()) {
            std::string value = it->second;
            if (value == "true" || value == "1" || value == "yes" || value == "on") {
                return true;
            } else if (value == "false" || value == "0" || value == "no" || value == "off") {
                return false;
            }
        }
        return defaultValue;
    }
    
private:
    std::unordered_map<std::string, std::string> _config;
}; 