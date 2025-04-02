#ifndef ROUND_AND_CAST
#define ROUND_AND_CAST

/*
 *  Round floating values to next int
 */
template <typename U, typename T>
T round_and_cast_up(const U x)
{
    return static_cast<T>(ceil(x));
}

/*
 *  Round floating values to last int
 */
template <typename U, typename T>
T round_and_cast_down(const U x)
{
    return static_cast<T>(floor(x));
}

#endif