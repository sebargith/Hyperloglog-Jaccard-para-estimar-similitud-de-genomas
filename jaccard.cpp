#include <iostream>
#include <fstream>
#include <unordered_set>
#include <vector>
#include "hyperloglog.h"

// Función para generar k-mers de una secuencia
std::unordered_set<std::string> generateKMers(const std::string& sequence, int k) {
    std::unordered_set<std::string> kmers;
    for (size_t i = 0; i <= sequence.size() - k; ++i) {
        kmers.insert(sequence.substr(i, k));
    }
    return kmers;
}

// Función para calcular la similitud de Jaccard
double jaccardSimilarity(HyperLogLog& hllA, HyperLogLog& hllB) {
    double estimateA = hllA.estimate();
    double estimateB = hllB.estimate();

    // Fusionar ambos HyperLogLog para estimar la unión
    HyperLogLog hllUnion = hllA;  // Copiamos hllA en hllUnion
    hllUnion.merge(hllB);         // Fusionamos hllB en hllUnion
    double estimateUnion = hllUnion.estimate();

    // Calcular la similitud de Jaccard
    return (estimateA + estimateB - estimateUnion) / estimateUnion;
}

// Funcion para leer el archivo
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
                break;  // Solo procesamos el número de genomas solicitado (ver esto, por si es necesario)
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

int main() {
    std::string filename = "GCF_001969825.1_ASM196982v1_genomic.fna";
    int numGenomes = 5;  // Procesar al menos 5 genomas (o mas?)
    int k = 20;  // Valor de k para los k-mers según enunciado de la profe (ponerlo en comando?)

    // Leer los genomas del archivo
    std::vector<std::string> genomes = readGenomesFromFile(filename, numGenomes);

    if (genomes.size() < 2) {
        std::cerr << "No hay suficientes genomas para comparar." << std::endl;
        return 1;
    }

    // Vamos a comparar los genomas de manera par a par
    for (size_t i = 0; i < genomes.size(); ++i) {
        for (size_t j = i + 1; j < genomes.size(); ++j) {
            std::cout << "Comparando genoma " << i + 1 << " con genoma " << j + 1 << std::endl;

            // Generamos los k-mers para ambos genomas
            auto kmersA = generateKMers(genomes[i], k);
            auto kmersB = generateKMers(genomes[j], k);

            // Creamos instancias de HyperLogLog para cada genoma
            HyperLogLog hllA, hllB;
            for (const auto& kmer : kmersA) {
                hllA.add(kmer);
            }
            for (const auto& kmer : kmersB) {
                hllB.add(kmer);
            }

            // Calculamos la similitud de Jaccard entre ambos
            double jaccard = jaccardSimilarity(hllA, hllB);

            std::cout << "Similitud de Jaccard entre genoma " << i + 1 << " y genoma " << j + 1 << ": " << jaccard << std::endl;
        }
    }

    return 0;
}
