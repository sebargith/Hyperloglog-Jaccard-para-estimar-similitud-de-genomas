#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>

using namespace std;

// Function to read genomes from a file
vector<string> readGenomesFromFile(const string& filename) {
    vector<string> genomes;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "No se pudo abrir el archivo: " << filename << endl;
        return genomes;
    }

    string line, genome;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '>') {
            if (!genome.empty()) {
                genomes.push_back(genome);
                genome.clear();
            }
        } else {
            genome += line;
        }
    }
    if (!genome.empty()) {
        genomes.push_back(genome);
    }

    file.close();
    return genomes;
}

// Function to extract minimizers from a genome sequence
unordered_set<string> extractMinimizers(const string& genome, int k, int w) {
    unordered_set<string> minimizers;
    if (genome.size() < k) {
        return minimizers; // Return empty set if genome is too short
    }
    for (size_t i = 0; i <= genome.size() - w; ++i) {
        string window = genome.substr(i, w);
        string minimizer = window.substr(0, k);
        for (size_t j = 1; j <= window.size() - k; ++j) {
            string kmer = window.substr(j, k);
            if (kmer < minimizer) {
                minimizer = kmer;
            }
        }
        minimizers.insert(minimizer);
    }
    return minimizers;
}

// Function to compute Jaccard similarity between two sets of minimizers
double computeJaccardSimilarity(const unordered_set<string>& set1, const unordered_set<string>& set2) {
    unordered_set<string> intersection;
    unordered_set<string> unionSet = set1;

    for (const auto& minimizer : set2) {
        if (set1.find(minimizer) != set1.end()) {
            intersection.insert(minimizer);
        }
        unionSet.insert(minimizer);
    }

    return static_cast<double>(intersection.size()) / unionSet.size();
}

int main() {
    string filename = "GCF_001969825.1_ASM196982v1_genomic.fna";
    int k = 20;  // Valor de k para los k-mers
    int w = 30;  // Tamaño de la ventana de los minimizer

    // Leer los genomas del archivo
    vector<string> genomes = readGenomesFromFile(filename);

    if (genomes.size() < 2) {
        cerr << "No hay suficientes genomas para comparar." << endl;
        return 1;
    }

    // Calcular los minimizadores para cada genoma
    vector<unordered_set<string>> minimizersList;
    for (const auto& genome : genomes) {
        minimizersList.push_back(extractMinimizers(genome, k, w));
    }

    // Computar Jaccard entre cada par de genomas
    for (size_t i = 0; i < minimizersList.size(); ++i) {
        for (size_t j = i + 1; j < minimizersList.size(); ++j) {
            double similarity = computeJaccardSimilarity(minimizersList[i], minimizersList[j]);
            cout << "Jaccard similarity between genome " << i << " and genome " << j << ": " << similarity << endl;
        }
    }

    return 0;
}