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
$ make check
```


```
// in "example.cpp"

#include <iostream>
#include <hll/hyperloglog.hpp>

int main() {
  hll::hyperloglog<18, 25> h;

  for (std::size_t i = 1; i <= 10'000'000ul; i++)
    h.insert(i);

  std::cout << "Estimate is " << h.estimate() << std::endl;
  return 0;
}
```

Assuming you cloned this library in `/path/to/hyperloglog` you can compile
and run `example.cpp` with:

```bash
$ g++ -std=c++14 \
    -I /path/to/hyperloglog/include \
    -isystem /path/to/hyperloglog/dep/MurmurHash3/include \
    -o example.o example.cpp
$ g++ example.o -o example
$ time ./example
Estimate is 1.00296e+07

real    0m0.669s
user    0m0.661s
sys     0m0.000s
```

On a relativly recent comodity CPU `hll::hyperloglog` enables you to insert and
calculate cardinality of ten million items in less than a second.

See more examples of `hll::hyperloglog` in the tests located at
`src/hll_tests.cpp`.
