#ifndef HYPERLOGLOG_H
#define HYPERLOGLOG_H

#include <vector>
#include <string>

class HyperLogLog {
private:
    std::vector<int> registers;
    static const int p = 14; // 2^14 buckets
    static const int m = 1 << p;

public:
    HyperLogLog();  // Constructor

    // Añadir un elemento al HyperLogLog
    void add(const std::string &data);

    // Estimar la cardinalidad
    double estimate();

    // Función hash utilizando SpookyHash
    static uint32_t hash(const std::string &data);

    // Contar ceros a la izquierda
    static int countLeadingZeros(uint32_t hashValue);

    // Método para fusionar dos HyperLogLog
    void merge(const HyperLogLog &other);
};

#endif
