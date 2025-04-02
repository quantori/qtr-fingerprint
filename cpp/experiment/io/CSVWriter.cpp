#include "CSVWriter.h"
#include <sstream>

CSVWriter::CSVWriter() : delimiter_(',') {
}

CSVWriter::CSVWriter(char delimiter) : delimiter_(delimiter) {
}

CSVWriter::~CSVWriter() {
}

bool CSVWriter::writeToFile(const std::string& filename, const CSVTable& table, bool includeHeader) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    return writeToStream(file, table, includeHeader);
}

bool CSVWriter::writeToStream(std::ostream& stream, const CSVTable& table, bool includeHeader) {
    if (!stream) {
        return false;
    }
    
    // Write header if requested
    if (includeHeader) {
        const std::vector<std::string>& columnNames = table.getColumnNames();
        for (size_t i = 0; i < columnNames.size(); ++i) {
            if (i > 0) {
                stream << delimiter_;
            }
            stream << escapeField(columnNames[i]);
        }
        stream << "\n";
    }
    
    // Write data rows
    for (size_t row = 0; row < table.getRowCount(); ++row) {
        for (size_t col = 0; col < table.getColumnCount(); ++col) {
            if (col > 0) {
                stream << delimiter_;
            }
            stream << escapeField(table.getValue(row, col));
        }
        stream << "\n";
    }
    
    return stream.good();
}

void CSVWriter::setDelimiter(char delimiter) {
    delimiter_ = delimiter;
}

std::string CSVWriter::escapeField(const std::string& field) const {
    // If the field contains delimiter, newline, or quotes, we need to escape it
    bool needsQuotes = field.find(delimiter_) != std::string::npos ||
                       field.find('\n') != std::string::npos ||
                       field.find('"') != std::string::npos;
    
    if (!needsQuotes) {
        return field;
    }
    
    // Escape by wrapping in quotes and doubling any quotes inside
    std::ostringstream escaped;
    escaped << '"';
    
    for (char c : field) {
        if (c == '"') {
            escaped << "\"\""; // Double quotes to escape them
        } else {
            escaped << c;
        }
    }
    
    escaped << '"';
    return escaped.str();
}
