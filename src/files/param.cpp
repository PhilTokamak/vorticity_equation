#include "param.h"
#include "utility.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include <functional>
#include <stdexcept>

// Default value for global parameters
int N = 128;
double L = 2 * M_PI;
double dt = 0.01;
int STEPS = 1000;
int SAVE_EVERY = 100;
double nu = 5.0e-7;
Scheme scheme = ARAKAWA_2;

double dx = L / N;


void load_parameters(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error: cannot open " << filename << std::endl;
        return;
    }

    // Create a mapping from parameter name to the assignment function
    std::unordered_map<std::string, std::function<void(const std::string&)>> setters;
    setters["N"] = [](const std::string& value){ N = std::stoi(value); };
    setters["L"] = [](const std::string& value){ L = std::stod(value); };
    setters["dt"] = [](const std::string& value){ dt = std::stod(value); };
    setters["STEPS"] = [](const std::string& value){ STEPS = std::stoi(value); };
    setters["SAVE_EVERY"] = [](const std::string& value){ SAVE_EVERY = std::stoi(value); };
    setters["nu"] = [](const std::string& value){ nu = std::stod(value); };
    setters["scheme"] = [](const std::string& value) {
        if (value == "ARAKAWA_2") {
            scheme = ARAKAWA_2;
        } else if (value == "ARAKAWA_4") {
            scheme = ARAKAWA_4;
        } else if (value == "CENTERED_2") {
            scheme = CENTERED_2;
        } else {
            throw std::invalid_argument("Unknown scheme: " + value +
                                        " (in file " + __FILE__ +
                                        ", line " + std::to_string(__LINE__) + ")");
        };
    };

    std::string line;
    while (std::getline(infile, line)) {
        // Neglect empty lines or lines beginning with # or //
        if (line.empty() ||
            line[0] == '#' ||
            (line[0] == '/' && line[1] == '/')) continue;

        std::istringstream iss(line);
        std::string key;
        if (std::getline(iss, key, '=')) {
            std::string value;
            if (std::getline(iss, value)) {
                // Remove leading and trailing white spaces
                key = trim_spaces(key);
                value = trim_spaces(value);

                // Find the mapping, and call the coresponding setup function
                auto it = setters.find(key);
                if (it != setters.end()) {
                    it->second(value);
                } else {
                    throw std::invalid_argument("Unkown parameter: " + key +
                                                " (in file " + __FILE__ +
                                                ", line " + std::to_string(__LINE__) + ")");
                }
            }
        }
    }

    // Update the value of dx
    dx = L / N;
}