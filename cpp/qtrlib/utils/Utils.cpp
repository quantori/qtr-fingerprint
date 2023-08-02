#include "Utils.h"

#include <stdexcept>

using namespace std;

namespace qtr {

    bool endsWith(const string &a, const string &b) {
        return a.size() >= b.size() &&
               (0 == a.compare(a.size() - b.size(), b.size(), b));
    }

    int chexToInt(char letter) {
        return letter >= '0' && letter <= '9' ? letter - '0' : letter - 'a' + 10;
    }

    vector <filesystem::path> findFiles(const filesystem::path &pathToDir, string extension) {
        vector<filesystem::path> files;
        if (!extension.empty() && extension[0] != '.')
            extension = "." + extension;
        for (const auto &entry: filesystem::recursive_directory_iterator(pathToDir)) {
            if (entry.path().stem() != ".DS_Store" && entry.path().extension() == extension) {
                files.push_back(entry.path());
            }
        }
        return files;
    }

    void askAboutContinue(const string &question) {
        cout << question << ". Continue? (y/N): ";
        string userAnswer;
        cin >> userAnswer;
        if (userAnswer != "Y" && userAnswer != "y") {
            cout << "abort\n";
            exit(-1);
        }
    }

    void
    initLogging(char **argv, google::LogSeverity severity, const char *base_filename, bool alsoLogToStderr) {
        google::InitGoogleLogging(argv[0]);
        google::SetLogDestination(severity, base_filename);
        FLAGS_alsologtostderr = alsoLogToStderr;
    }

    void logErrorAndExit(const string &message) {
        LOG(ERROR) << message;
        exit(-1);
    }

    void copyFileAndCheck(const std::filesystem::path &from, const std::filesystem::path &to) {
        if (!filesystem::copy_file(from, to))
            logErrorAndExit("Cannot copy file from " + from.string() + " to " + to.string());
    }
} // namespace qtr