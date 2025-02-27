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
    
    return 0;
}
