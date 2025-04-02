#ifndef STATE
#define STATE

/*
 *  State for a particular buffer
 */
template <typename DataT>
struct State
{
    int pointer;        // Pointer to a particular event
    DataT usage;        // Memory usage of a particular buffer
    DataT transferrate; // Transfert rate of a particular buffer
};

#endif