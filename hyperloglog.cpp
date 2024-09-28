#include <iostream>
#include <vector>
#include <cmath>
#include <cstdint>
#include <climits>
#include <algorithm> 
#include "hyperloglog.h"
#include "Spooky.h" 

// Constructor inicializando los registros con ceros
HyperLogLog::HyperLogLog() : registers(m, 0) {}

// Usamos __buildtin_clz__ para contar los ceros a la izquierda
int HyperLogLog::countLeadingZeros(uint32_t hashValue) {
    return __builtin_clz(hashValue);
}

// Función hash utilizando SpookyHash
uint32_t HyperLogLog::hash(const std::string &data) {
    // Usar SpookyHash32 para calcular el hash
    return SpookyHash::Hash32(data.c_str(), data.size(), 0);
}

// Añadir un elemento al HyperLogLog
void HyperLogLog::add(const std::string &data) {
    uint32_t hashValue = hash(data);

    // Extraer el índice del registro (los primeros p bits del hash)
    int registerIndex = hashValue >> (32 - p);

    // Contar los ceros en el resto de los bits sin desplazamiento
    int leadingZeros = countLeadingZeros(hashValue << p);  // Usa los bits restantes sin perder información

    // Actualizamos el registro con el máximo número de ceros encontrados
    registers[registerIndex] = std::max(registers[registerIndex], leadingZeros);
}


// Estimar la cardinalidad usando HyperLogLog
double HyperLogLog::estimate() {
    // Estimador para la cardinalidad
    double alphaMM = 0.7213 / (1 + 1.079 / m) * m * m;
    double harmonicSum = 0.0;

    for (int reg : registers) {
        harmonicSum += 1.0 / (1 << reg);
    }

    double rawEstimate = alphaMM / harmonicSum;

    // Aplicar correcciones
    if (rawEstimate <= (5.0 / 2.0) * m) {
        int zeroCount = std::count(registers.begin(), registers.end(), 0);
        if (zeroCount > 0) {
            return m * std::log(static_cast<double>(m) / zeroCount);
        }
    } else {
        uint64_t max_value = uint64_t(1) << 32;
        if (rawEstimate > max_value / 30.0) {
            return -double(max_value) * std::log(1.0 - rawEstimate / max_value);
        }
    }

    return rawEstimate;
}

// Función para fusionar dos HyperLogLog
void HyperLogLog::merge(const HyperLogLog &other) {
    for (size_t i = 0; i < registers.size(); ++i) {
        registers[i] = std::max(registers[i], other.registers[i]);
    }
}
