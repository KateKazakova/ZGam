#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <string>

#include "ConfigReader.h"


ConfigReader::ConfigReader(std::string filename)	{
	std::ifstream config_file(filename);

    if (config_file.is_open())  {
        std::string line;

        while(getline(config_file, line)) {
            line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
            if(line.empty() || line[0] == '#')  {
                continue;
            }

            auto delimiterPos = line.find(":");
            auto name = line.substr(0, delimiterPos);
            auto value = line.substr(delimiterPos + 1);

            this -> config[name] = value;
        }
    }

    else    {
        std::cerr << "Couldn't open config file for reading.\n";
    }
};

bool ConfigReader::getBool(std::string varname)	{
	std::string var = this -> config[varname];

	if (var == "true")	return true;
	else return false;
};

int ConfigReader::getInt(std::string varname) {
    std::string var = this -> config[varname];

    return stoi(var);
};

double ConfigReader::getDouble(std::string varname) {
    std::string var = this -> config[varname];
  
    return stod(var);
};

int ConfigReader::getHex(std::string varname) {
    std::string var = this -> config[varname];

    return stoul(var, NULL, 16);
};

std::string ConfigReader::getString(std::string varname) {
    std::string var = this -> config[varname];

    return var;
};
