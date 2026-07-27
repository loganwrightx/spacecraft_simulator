#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>
#include <iostream>
#include <iomanip>

class Logger {
public:
    Logger(std::string filename);
    ~Logger();

    void WriteLine(std::string line);

    void Flush();

    void Close();

private:
    std::string lFilename;
    std::ofstream lLogFile;
};

#endif // LOGGER_H
