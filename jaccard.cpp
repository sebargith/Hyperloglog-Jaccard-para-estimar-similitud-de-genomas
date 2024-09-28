#include <iostream>
#include <fstream>
#include <unordered_set>
#include <vector>
#include <cmath>  
#include "hyperloglog.h"

// Función para generar k-mers de una secuencia
std::unordered_set<std::string> generateKMers(const std::string& sequence, int k) {
    std::unordered_set<std::string> kmers;
    for (size_t i = 0; i <= sequence.size() - k; ++i) {
        kmers.insert(sequence.substr(i, k));
    }
    return kmers;
}

// Función para calcular la similitud de Jaccard real
double realJaccard(const std::unordered_set<std::string>& kmersA, const std::unordered_set<std::string>& kmersB) {
    int intersectionSize = 0;
    for (const auto& kmer : kmersA) {
        if (kmersB.find(kmer) != kmersB.end()) {
            intersectionSize++;
        }
    }
    int unionSize = kmersA.size() + kmersB.size() - intersectionSize;
    return static_cast<double>(intersectionSize) / unionSize;
}

// Función para calcular la similitud de Jaccard estimada usando HyperLogLog
double jaccardSimilarity(HyperLogLog& hllA, HyperLogLog& hllB) {
    double estimateA = hllA.estimate();
    double estimateB = hllB.estimate();

    // Fusionar ambos HyperLogLog para estimar la unión
    HyperLogLog hllUnion = hllA;
    hllUnion.merge(hllB);
    double estimateUnion = hllUnion.estimate();

    // Calcular la similitud de Jaccard asegurando que no sea negativa
    double jaccard = std::max(0.0, (estimateA + estimateB - estimateUnion) / estimateUnion);

    return jaccard;
}

// Funcion para leer el archivo de genomas
std::vector<std::string> readGenomesFromFile(const std::string& filename, int numGenomes) {
    std::ifstream infile(filename);
    std::string line;
    std::vector<std::string> genomes;
    std::string currentGenome;
    int genomeCount = 0;

    while (std::getline(infile, line)) {
        if (line.empty()) continue;  // Saltar líneas vacías

        // Si encontramos un nuevo genoma (línea que empieza con '>')
        if (line[0] == '>') {
            if (!currentGenome.empty()) {
                genomes.push_back(currentGenome);
                currentGenome.clear();
                genomeCount++;
            }

            if (genomeCount >= numGenomes) {
                break;  // Solo procesamos el número de genomas solicitado
            }
        } else {
            // Concatenar las líneas que contienen la secuencia de ADN
            currentGenome += line;
        }
    }

    // Agregar el último genoma si es necesario
    if (!currentGenome.empty() && genomeCount < numGenomes) {
        genomes.push_back(currentGenome);
    }

    return genomes;
}

// Función para calcular ERM y EAM para una comparación específica
void CalculodeErrores(double realJaccard, double estimatedJaccard) {
    double erm = 0.0;
    double eam = 0.0;

    // Evitar dividir por 0 en ERM
    if (realJaccard != 0) {
        erm = fabs(estimatedJaccard - realJaccard) / realJaccard;
    }

    // Calcular EAM
    eam = fabs(estimatedJaccard - realJaccard);

    std::cout << "Error Relativo Medio (ERM): " << (erm == 0.0 ? "No se puede calcular" : std::to_string(erm)) << std::endl;
    std::cout << "Error Absoluto Medio (EAM): " << eam << std::endl;
}

int main() {
    std::string filename = "GCF_001969825.1_ASM196982v1_genomic.fna";
    int numGenomes = 5;  // Procesar al menos 5 genomas
    int k = 20;  // Valor de k para los k-mers

    // Leer los genomas del archivo
    std::vector<std::string> genomes = readGenomesFromFile(filename, numGenomes);

    if (genomes.size() < 2) {
        std::cerr << "No hay suficientes genomas para comparar." << std::endl;
        return 1;
    }

    // Comparar los genomas de manera par a par
    for (size_t i = 0; i < genomes.size(); ++i) {
        for (size_t j = i + 1; j < genomes.size(); ++j) {
            std::cout << "Comparando genoma " << i + 1 << " con genoma " << j + 1 << std::endl;

            // Generamos los k-mers para ambos genomas
            auto kmersA = generateKMers(genomes[i], k);
            auto kmersB = generateKMers(genomes[j], k);

            // Calcular Jaccard real
            double realJ = realJaccard(kmersA, kmersB);
            std::cout << "Similitud de Jaccard real entre genoma " << i + 1 << " y genoma " << j + 1 << ": " << realJ << std::endl;

            // Creamos instancias de HyperLogLog para cada genoma
            HyperLogLog hllA, hllB;
            for (const auto& kmer : kmersA) {
                hllA.add(kmer);
            }
            for (const auto& kmer : kmersB) {
                hllB.add(kmer);
            }

            // Calcular Jaccard estimado
            double estimatedJ = jaccardSimilarity(hllA, hllB);
            std::cout << "Similitud de Jaccard estimada entre genoma " << i + 1 << " y genoma " << j + 1 << ": " << estimatedJ << std::endl;

            // Calcular y mostrar errores
            CalculodeErrores(realJ, estimatedJ);
        }
    }

    return 0;
}
