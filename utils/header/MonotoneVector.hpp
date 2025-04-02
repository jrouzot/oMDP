
#ifndef __MONOTONEVECTOR_HPP
#define __MONOTONEVECTOR_HPP

#include <vector>

/**********************************************
 * MONOTONEVECTOR
 **********************************************/
/// Sparse set representation


template <typename T> class MonotoneVector {
protected:
  std::vector<T> owner_;
  std::size_t    end_{0};
public:
  MonotoneVector() 
  {
     // do some range checks
  }

  size_t size() const {
    return end_;
  }
	
	T& operator[](const int i) {
		return owner_[i];
	}

	void resize(const size_t s) {
		end_ = s;
		if(owner_.size() < end_)
			owner_.resize(end_);
	}
	
	typename std::vector<T>::iterator begin() {
		return owner_.begin();
	}
	
	typename std::vector<T>::iterator end() {
		return owner_.begin() + end_;
	}
	
	typename std::vector<T>::reverse_iterator rbegin() {
          return owner_.rend() - end_;
        }

        typename std::vector<T>::reverse_iterator rend() {
          return owner_.rend();
        }

        // void push_back(T& e) {
        // 	if(end_ < owner_.size())
        // 		owner_[_end++] = e;
        // 	else
        // 		owner_.push_back(e);
        // }
}; 

#endif
