

#include <logger.h>

Logger::Logger(std::string filename) {
    lFilename = filename;
    lLogFile = std::ofstream(filename);
}

Logger::~Logger() {
    Close();
}

void Logger::WriteLine(std::string line) {
    lLogFile << line << '\n';
}

void Logger::Flush() {
    lLogFile << std::flush;
}

void Logger::Close() {
    lLogFile.close();
}
