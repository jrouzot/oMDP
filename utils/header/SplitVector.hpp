#ifndef SPLIT_VECTOR
#define SPLIT_VECTOR

#include <vector>
#include <algorithm>

#include "PrintVector.hpp"

template <typename T>
std::vector<std::vector<T>> splitVector(const std::vector<T>& source, const std::vector<int>& indexes) {
    std::vector<std::vector<T>> result;

    if(!indexes.size())
        return {source};

    result.reserve(indexes.size());
    auto anchor_front = source.begin();
    for (auto b : indexes) {
        auto anchor_end = std::next(source.begin(), b + 1);
        result.emplace_back(anchor_front, anchor_end);
        anchor_front = anchor_end;
    }
    return result;
}

#endif