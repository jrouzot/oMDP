#ifndef EVENT
#define EVENT

/*
 * An event occur at a time t for a particular buffer b
 */
template <typename TimeT>
class Event
{

public:
    TimeT time;
    int buffer;

    Event(TimeT t, int b);

    bool operator<(const Event<TimeT> &e) const;
    bool operator==(const Event<TimeT> &e) const;
    bool operator<=(const Event<TimeT> &e) const;
    bool operator>(const Event<TimeT> &e) const;
    bool operator>=(const Event<TimeT> &e) const;
};

template <typename TimeT>
Event<TimeT>::Event(TimeT t, int b) : time(t), buffer(b) {}

template <typename TimeT>
bool Event<TimeT>::operator<(const Event<TimeT> &e) const
{
    return this->time < e.time;
}

template <typename TimeT>
bool Event<TimeT>::operator==(const Event<TimeT> &e) const
{
    return this->time == e.time;
}

template <typename TimeT>
bool Event<TimeT>::operator<=(const Event<TimeT> &e) const
{
    return this->time <= e.time;
}

template <typename TimeT>
bool Event<TimeT>::operator>(const Event<TimeT> &e) const
{
    return this->time > e.time;
}

template <typename TimeT>
bool Event<TimeT>::operator>=(const Event<TimeT> &e) const
{
    return this->time >= e.time;
}

#endif