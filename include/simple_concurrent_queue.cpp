/* This is a simple thread-safe queue implemented using the C++ standard library.
 * It's a simpler alternative to the tbb::concurrent_queue
 * Taken directly from: https://juanchopanzacpp.wordpress.com/2013/02/26/concurrent-queue-c11/
*/

//
// Copyright (c) 2013 Juan Palacios juan.palacios.puyana@gmail.com
// Subject to the BSD 2-Clause License
// - see < http://opensource.org/licenses/BSD-2-Clause>
//

#ifndef CONCURRENT_QUEUE_
#define CONCURRENT_QUEUE_

#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>

template <typename T>
class SimpleConcurrentQueue {
 public:

  T pop(){
    std::unique_lock<std::mutex> mlock(mutex_);
    while (queue_.empty()) {
      cond_.wait(mlock);
    }
    auto val = queue_.front();
    queue_.pop();
    return val;
  }

  void pop(T& item){
    std::unique_lock<std::mutex> mlock(mutex_);
    while (queue_.empty())
    {
      cond_.wait(mlock);
    }
    item = queue_.front();
    queue_.pop();
  }

  void push(const T& item){
    std::unique_lock<std::mutex> mlock(mutex_);
    queue_.push(item);
    mlock.unlock();
    cond_.notify_one();
  }

  bool empty(){
    std::unique_lock<std::mutex> mlock(mutex_);
    bool res = queue_.empty();
    mlock.unlock();
    return res;
  }

  int64_t size(){
    std::unique_lock<std::mutex> mlock(mutex_);
    return queue_.size();
  }

  SimpleConcurrentQueue()=default;
  SimpleConcurrentQueue(const SimpleConcurrentQueue&) = delete;            // disable copying
  SimpleConcurrentQueue& operator=(const SimpleConcurrentQueue&) = delete; // disable assignment
  
 private:
  std::queue<T> queue_;
  std::mutex mutex_;
  std::condition_variable cond_;
};

#endif
