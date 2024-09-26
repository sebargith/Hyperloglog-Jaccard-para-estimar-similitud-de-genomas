#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <cstdint>
#include "hyperloglog.h"

// Funcion para extraer minimizadores de una secuencia genómica
std::unordered_set<std::string> extractMinimizers(const std::string& genome, int k, int w) {
    std::unordered_set<std::string> minimizers;
    for (size_t i = 0; i <= genome.size() - k; ++i) {
        std::string kmer = genome.substr(i, k);
        std::string minimizer = kmer.substr(0, w);
        for (size_t j = 1; j <= k - w; ++j) {
            std::string window = kmer.substr(j, w);
            if (window < minimizer) {
                minimizer = window;
            }
        }
        minimizers.insert(minimizer);
    }
    return minimizers;
}

// Function para computar la similitud de Jaccard entre dos conjuntos de minimizadores
double computeJaccardSimilarity(const std::unordered_set<std::string>& set1, const std::unordered_set<std::string>& set2) {
    std::unordered_set<std::string> intersection;
    std::unordered_set<std::string> unionSet = set1;

    for (const auto& minimizer : set2) {
        if (set1.find(minimizer) != set1.end()) {
            intersection.insert(minimizer);
        }
        unionSet.insert(minimizer);
    }

    return static_cast<double>(intersection.size()) / unionSet.size();
}

// Funcion para leer los genomas del archivo especificado
std::vector<std::string> readGenomesFromFile(const std::string& filename, int numGenomes) {
    std::vector<std::string> genomes;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return genomes;
    }

    std::string line, genome;
    while (std::getline(file, line) && genomes.size() < numGenomes) {
        if (line.empty() || line[0] == '>') {
            if (!genome.empty()) {
                genomes.push_back(genome);
                genome.clear();
            }
        } else {
            genome += line;
        }
    }
    if (!genome.empty() && genomes.size() < numGenomes) {
        genomes.push_back(genome);
    }

    file.close();
    return genomes;
}

int main() {
    std::string filename = "GCF_001969825.1_ASM196982v1_genomic.fna";
    int numGenomes = 5;  // Procesar al menos 5 genomas (o mas?)
    int k = 20;  // Valor de k para los k-mers según enunciado de la profe (ponerlo en comando?)
    int w = 10;  // Tamaño de la ventana de los minimizer

    // Leer los genomas del archivo
    std::vector<std::string> genomes = readGenomesFromFile(filename, numGenomes);

    if (genomes.size() < 2) {
        std::cerr << "No hay suficientes genomas para comparar." << std::endl;
        return 1;
    }

    // Calcular los minimizadores para cada genoma
    std::vector<std::unordered_set<std::string>> minimizersList;
    for (const auto& genome : genomes) {
        minimizersList.push_back(extractMinimizers(genome, k, w));
    }

    // Computar Jaccard entre cada par de genomas
    for (size_t i = 0; i < minimizersList.size(); ++i) {
        for (size_t j = i + 1; j < minimizersList.size(); ++j) {
            double similarity = computeJaccardSimilarity(minimizersList[i], minimizersList[j]);
            std::cout << "Jaccard similarity between genome " << i << " and genome " << j << ": " << similarity << std::endl;
        }
    }

    return 0;
}
