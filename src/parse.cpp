#include "parse.hpp"

std::vector <long double> line_to_double_array(std::string str) {
    std::istringstream input(str);
    std::vector <long double> coords;
    while (!input.eof()) {
        std::string substring;
        input >> substring;
        coords.push_back(std::stold(substring));
    }
    return coords;
}

std::vector <std::vector <long double>> parse_file(std::string file_name) {
    std::vector <std::vector <long double>> coords;

    std::ifstream input(file_name);

    std::string line;
    while (!input.eof()) {
        std::getline(input, line);
        std::getline(input, line);
        for(size_t i = 0; i < line.length(); i++) {
            if(line[i] == '=')line[i] = ' ';
            if(line[i] != '-' && !isdigit(line[i]) && line[i] != '.' && line[i] != 'E' && line[i] != '+') line[i] = ' ';
        }
        coords.push_back(line_to_double_array(line));
    }
    input.close();
    return coords;
}
/*int main() {
    std::vector <std::vector <long double>> result;
    result = parse_file("2023bu.txt");
    for(int i = 0; i < result.size(); i++) {
        for(int j = 0; j < result[i].size(); j++) {
            std::cout<<result[i][j] << " ";
        }
        std::cout<<std::endl;
    }
    return 0;
}*/