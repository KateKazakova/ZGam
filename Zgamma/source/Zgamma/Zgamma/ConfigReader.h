#include <string>
#include <map>


class ConfigReader  {
    public:
        ConfigReader(std::string);

        bool getBool(std::string);
        int  getInt(std::string);
        int  getHex(std::string);
        double  getDouble(std::string);
        std::string getString(std::string);
    private:
        std::map <std::string, std::string> config;
};
