# `hll::hyperloglog`: HyperLogLog++ with C++14 [![Actions Status](https://github.com/arashbm/hyperloglog/workflows/Tests/badge.svg)](https://github.com/arashbm/hyperloglog/actions)

[HyperLogLog][hll] is a probabilistic data structure that can help you estimate
cardinality of very large multisets with a pre-determined accuracy using
constant space.

Typical cardinality calculation methods, say using Pythons `collections.Counter`
or using C++'s `std::unordered_map`, need O(n) space to count n unique elements.
HyperLogLog on the other hand requires a constant O(1) amount of space for
estimating cardinality of multisets with billions of items within a
pre-determined range of accuracy. For example, you can use HyperLogLog to
estimate the number of unique IP addresses that connect to your web server or
the number of unique words in a book to within a percent of the actual value
all with a few kilobytes of memory.

The tradeoff here is between the amount of constant space you allocate to the
HyperLogLog data structure and the accuracy as determined by the relative error.
The more *registers* you use, the more accurate your estimations are going to
be. Each register is represented here by a `uint8_t`, although a maximum of 6
bits of each register is ever used. A HyperLogLog data structure with
2<sup>m</sup> registers has a relative error or 1.04/√m. This means that a
HyperLogLog data structure with m=12 uses 2<sup>12</sup> registers (~4kB) and
has a relative error of 1.6% for large multisets.

[hll]: https://en.wikipedia.org/wiki/HyperLogLog

## Installation
This package relies on cmake. Make sure a moderatly recent version (3.14 or
newer) is already installed. You should also have a C++ compiler with C++14
support.


Here is how to build and run the tests:
```bash
$ git clone https://github.com/arashbm/hyperloglog.git
$ cd hyperloglog
$ mkdir build  # make directory to build in
$ cd build
$ cmake ..
$ cmake --build . --target hyperloglog_tests  # build the tests
$ ./hyperloglog_tests
```


At this point you should see "All tests passed" if all the steps are successful.
You can continue by install the library on your system:
```bash
cmake --build . --target install
```
This final step might require admin previlages.

## Example
```cpp
// in "example.cpp"

#include <iostream>
#include <hll/hyperloglog.hpp>

int main() {
  hll::hyperloglog<std::size_t, 18, 25> h;

  for (std::size_t i = 1; i <= 10'000'000ul; i++)
    h.insert(i);

  std::cout << "Estimate is " << h.estimate() << std::endl;
  return 0;
}
```

Assuming you have installed the library as instructed above, you can compile
and run `example.cpp` with:

```bash
$ g++ -std=c++14 -lhyperloglog -o example example.cpp
$ time ./example
Estimate is 1.00296e+07

real    0m0.669s
user    0m0.661s
sys     0m0.000s
```

On a relativly recent comodity CPU `hll::hyperloglog` enables you to insert and
calculate cardinality of ten million items in less than a second.

See more examples of `hll::hyperloglog` in the tests located at
`src/tests.cpp`.
