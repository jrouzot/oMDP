#ifndef PRINT_VECTOR
#define PRINT_VECTOR

#include <vector>
#include "../header/SimulationSingleWindow.hpp"

/*
 *  Print 1D vector
 */
template <typename T>
void print_vec(vector<T> vec)
{
    for(auto item : vec)
    {
        std::cout << item << " "; 
    }
    std::cout << std::endl;
}

/*
 *  Print 2D vector
 */
template <typename T>
void print_2Dvec(vector<vector<T>> vec)
{
    for(auto item : vec)
    {
        for(auto i : item)
        {
            std::cout << i << " "; 
        }
        std::cout << std::endl;
    }
}

/*
 *  Print event vector
 */
template <typename T>
void print_events(vector<vector<T>> vec)
{
    for(auto item : vec)
    {
        for(auto i : item)
        {
            std::cout << "(" << i.time << ", " << i.value << ") "; 
        }
        std::cout << std::endl;
    }
}

#endif