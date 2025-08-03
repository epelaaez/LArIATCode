#include <iostream>
#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;

int main() {
    std::string directory_path = "/pnfs/lariat/persistent/users/epelaez/reco_files/"; 
    std::string output_file = "files.list";
    
    std::ofstream list_file(output_file);
    if (!list_file) {
        std::cerr << "Error opening output file." << std::endl;
        return 1;
    }
    
    for (const auto& entry : fs::directory_iterator(directory_path)) {
        if (entry.is_regular_file()) {
            list_file << directory_path + entry.path().filename().string() << std::endl;
        }
    }
    
    list_file.close();
    std::cout << "File list saved to " << output_file << std::endl;
    
    std::string nn_directory_path = "/pnfs/lariat/persistent/users/epelaez/reco_nn_files/";
    std::string output_nn_file = "nn_files.list";

    std::ofstream nn_list_file(output_nn_file);
    if (!nn_list_file) {
        std::cerr << "Error opening NN output file." << std::endl;
        return 1;
    }

    for (const auto& entry : fs::directory_iterator(nn_directory_path)) {
        if (entry.is_regular_file()) {
            nn_list_file << nn_directory_path + entry.path().filename().string() << std::endl;
        }
    }

    nn_list_file.close();
    std::cout << "NN file list saved to " << output_nn_file << std::endl;

    return 0;
}
