#pragma once

#include <string>
#include <vector>
#include <memory>
#include <fstream>

class CSVTable {
public:
    virtual ~CSVTable() = default;
    
    virtual size_t getColumnCount() const = 0;
    
    virtual size_t getRowCount() const = 0;
    
    virtual const std::vector<std::string>& getColumnNames() const = 0;
    
    virtual std::string getValue(size_t row, size_t col) const = 0;
};

class CSVWriter {
public:
    CSVWriter();
    explicit CSVWriter(char delimiter);
    ~CSVWriter();
    
    bool writeToFile(const std::string& filename, const CSVTable& table, bool includeHeader = true);
    
    bool writeToStream(std::ostream& stream, const CSVTable& table, bool includeHeader = true);
    
    void setDelimiter(char delimiter);
    
private:
    char delimiter_;
    
    std::string escapeField(const std::string& field) const;
};
