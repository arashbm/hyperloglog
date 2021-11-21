# hyperloglog [![Actions Status](https://github.com/arashbm/hyperloglog/workflows/Tests/badge.svg)](https://github.com/arashbm/hyperloglog/actions)
hyperloglog++ implementation with C++14

## Getting Started

Clone the library:
```bash
$ git clone https://github.com/arashbm/hyperloglog.git
```

Compile and run the tests:
```
$ cd hyperloglog
$ make
```


```
// in "example.cpp"

#include <iostream>
#include <hyperloglog.hpp>

int main() {
  auto h  = hll::hyperloglog<18, 25>();

  for (unsigned long int i = 1; i <= 10'000'000UL; i++)
    h.insert(i);

  std::cout << "Estimate is " << h.estimate() << std::endl;
  return 0;
}
```

Assuming you cloned this library in `/path/to/hyperloglog` you can compile
and run `example.cpp` with:

```bash
$ g++ -std=c++14 \
    -I/path/to/hyperloglog/include \
    -isystem /path/to/hyperloglog/dep/MurmurHash3/include \
    -c -o example.o example.cpp
$ g++ example.o -o example
$ ./example
Estimate is 1.00296e+07

```
