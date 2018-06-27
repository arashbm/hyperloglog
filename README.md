# HyperLogLog
HyperLogLog++ Implementation with C++14

## Getting Started
Clone the library:
```bash
$ git clone https://github.com/arashbm/HyperLogLog.git
```

Compile the library:
```
$ cd HyperLogLog
$ make
```


```
// in "example.cpp"

#include <iostream>
#include <hyperloglog.hpp>

int main() {
  auto h  = hll::HyperLogLog<18, 25>();

  for (unsigned long int i = 1; i <= 10'000'000UL; i++)
    h.insert(i);

  std::cout << "Estimate is " << h.estimate() << std::endl;
  return 0;
}
```

Assuming you cloned this library in `/path/to/HyperLogLog` and have already
compiled it, you can compile `example.cpp` with:

```bash
$ g++ -std=c++14 -I/path/to/HyperLogLog  -c -o example.o example.cpp
$ g++ example.o /path/to/HyperLogLog/MurmurHash3.o -o example

```
