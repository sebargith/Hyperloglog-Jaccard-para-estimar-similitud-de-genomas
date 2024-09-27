#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <cmath>

using namespace std;

// Lee genomas del archivo
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

// Funcion para extraer minimizadores de una secuencia de genoma
unordered_set<string> extractMinimizers(const string& genome, int k, int w) {
    unordered_set<string> minimizers;
    if (genome.size() < k) {
        return minimizers; // Devuelve set vacio si el genoma es muy corto
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

// Funcion para calcular Jaccard entre dos conjuntos de minimizadores
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

// Funcion para calcular el verdadero Jaccard entre dos conjuntos de k-mers
double computeTrueJaccardSimilarity(const unordered_set<string>& set1, const unordered_set<string>& set2) {
    unordered_set<string> intersection;
    unordered_set<string> unionSet = set1;

    for (const auto& kmer : set2) {
        if (set1.find(kmer) != set1.end()) {
            intersection.insert(kmer);
        }
        unionSet.insert(kmer);
    }

    return static_cast<double>(intersection.size()) / unionSet.size();
}

// Funcion para calcular el error relativo medio
double computeMeanRelativeError(const vector<double>& trueValues, const vector<double>& estimatedValues) {
    double totalRelativeError = 0.0;
    size_t validPairs = 0;
    for (size_t i = 0; i < trueValues.size(); ++i) {
        if (trueValues[i] != 0) {
            double relativeError = fabs(trueValues[i] - estimatedValues[i]) / trueValues[i];
            totalRelativeError += relativeError;
            validPairs++;
        }
    }
    return validPairs > 0 ? totalRelativeError / validPairs : 0.0;
}

// Funcion para computar el error absoluto medio
double computeMeanAbsoluteError(const vector<double>& trueValues, const vector<double>& estimatedValues) {
    double totalAbsoluteError = 0.0;
    for (size_t i = 0; i < trueValues.size(); ++i) {
        double absoluteError = fabs(trueValues[i] - estimatedValues[i]);
        totalAbsoluteError += absoluteError;
    }
    return totalAbsoluteError / trueValues.size();
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

    // Calcular los k-mers para cada genoma
    vector<unordered_set<string>> kmersList;
    for (const auto& genome : genomes) {
        kmersList.push_back(extractMinimizers(genome, k, k)); // Usar k como ventana para k-mers
    }

    // Computar Jaccard y errores relativos
    vector<double> trueJaccardValues;
    vector<double> estimatedJaccardValues;
    double totalRelativeError = 0.0;
    double totalAbsoluteError = 0.0;
    size_t pairCount = 0;

    for (size_t i = 0; i < minimizersList.size(); ++i) {
        for (size_t j = i + 1; j < minimizersList.size(); ++j) {
            double trueSimilarity = computeTrueJaccardSimilarity(kmersList[i], kmersList[j]);
            double estimatedSimilarity = computeJaccardSimilarity(minimizersList[i], minimizersList[j]);
            trueJaccardValues.push_back(trueSimilarity);
            estimatedJaccardValues.push_back(estimatedSimilarity);
            if (trueSimilarity != 0) {
                double relativeError = fabs(trueSimilarity - estimatedSimilarity) / trueSimilarity;
                totalRelativeError += relativeError;
                pairCount++;
            }
            double absoluteError = fabs(trueSimilarity - estimatedSimilarity);
            totalAbsoluteError += absoluteError;
            cout << "Estimacion de similitud de Jaccard entre genoma " << i << " y genoma " << j << ": " << estimatedSimilarity << endl;
        }
    }

    // Calcular el error relativo medio
    double meanRelativeError = pairCount > 0 ? totalRelativeError / pairCount : 0.0;
    cout << "Error relativo medio: " << meanRelativeError << endl;

    // Calcular el error absoluto medio
    double meanAbsoluteError = totalAbsoluteError / trueJaccardValues.size();
    cout << "Error absoluto medio: " << meanAbsoluteError << endl;

    return 0;
}