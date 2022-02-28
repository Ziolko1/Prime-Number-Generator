#include <iostream>     //IO
#include <cmath>        //std::sqrt
#include <chrono>       //time
#include <vector>       //container

/* We are skipping evens. Here are examples and necessary calculations
Index:    1   2   3   4   5   6   7   8   9   10  11  12  13  14  15
Value:    3   5   7   9   11  13  15  17  19  21  23  25  27  29  31
Calculating index from value            -> value = 2*index +1
Calcularing value from index            -> index = value/2 (integer division)
Calculating index of squared value      -> 2*index*index + 2*index -> 2 * index * (index+1)
Calculating jump to the next value      -> 2*index+1

Jumping over more primes (3, 5) does not gain any time on my PC,
This is most likely caused by additional overhead from more complex math

Optimizations:
Sieve with only primes.
Start sieveing from prime^2, lower multiples are taken care of by lower primes
This also means we can end sieving with sqrt of the limit.
TO DO:
Make windowed sieve multithreaded!
*/

constexpr uint64_t UPPER_BOUND {2'147'483'647}; //2^31-1
constexpr uint32_t L1_CACHE_SIZE {32'768};      //2^15

class Sieve
{
    std::vector<bool> m_sieve;

public:
    void print (uint64_t start, uint64_t end)
    {
        if (start <=2)
        {
            std::cout << "2 ";
            start = 2;
        }
        for (uint64_t i{start/2}; i<(end+1)/2; ++i)
            if (m_sieve.at(i))
                std::cout << i*2 +1 << ' ';
    }
    uint32_t CountPrimes (uint64_t start, uint64_t end)
    {
        uint32_t answer{0};
        if (start <=2 && 2 <= end)
        {
            ++answer;
            start = 2;
        }
        for (uint64_t i{start/2}; i<(end+1)/2; ++i)
        {
            if (m_sieve.at(i))
                ++answer;
        }
        return answer;
    }
    Sieve (uint64_t limit)
    {
        //Skip evens to save half of the memory
        uint64_t memory {limit/2+1};
        uint64_t checkingRange {static_cast<uint64_t> (std::sqrt(static_cast<double>(limit)))/2+1};
        std::vector<bool> answer (limit/2+1, 1);
        for (uint64_t i{1}; i<checkingRange; ++i)
        {
            if (answer[i])
                for (uint64_t j{(i+1)*i*2}; j<memory; j+=2*i+1)
                    answer[j] = false;
        }
        m_sieve = std::move(answer);
    }

    Sieve (const uint64_t limit, const uint32_t cache_size)
    {
        uint32_t sqrt {static_cast<uint32_t> (std::sqrt(static_cast<double> (limit))) +1};
        //Skip evens
        uint64_t memory {limit/2+1};
        uint32_t segment_size {8*cache_size}, index1 {1}, copyIndex {1};

        //Auxiliary sieve holds all primes that will be used to sieve the actual sieve
        //Maybe generating it "recursively" could save some time
        Sieve auxiliarySieve {Sieve (sqrt)};
        std::vector <char> sieve (cache_size);      //Bool vectors are about 40% slower
        std::vector <bool> answer (memory, 0);      //Collects all primes
        std::vector <uint32_t> primes;              //Holds values that we need to add
        std::vector <uint32_t> next;                //Holds starting values for the next frame

        for (uint64_t low = 0; low <= memory; low += cache_size)
        {
            uint64_t high {std::min(low + cache_size -1, memory)};
            std::fill(sieve.begin(), sieve.end(), 1);
            for (; (index1+1)*index1*2 <= high; ++index1)
            {
                if (auxiliarySieve.m_sieve[index1])
                {
                    primes.push_back(index1);
                    next.push_back((index1+1)*index1*2 - low);
                }
            }
            for (std::size_t i{0}; i<primes.size(); ++i)
            {
                uint32_t j{next[i]};
                for (uint32_t increment {primes[i]*2+1}; j<cache_size ; j+=increment)
                {
                    sieve[j] = false;
                }
                next[i] = j - cache_size;
            }

            for (; copyIndex <= high; ++copyIndex)
                if (sieve[copyIndex-low])
                    answer[copyIndex] = true;
        }
        m_sieve = std::move(answer);
    }
    static void tests ()
    {
        std::vector<uint32_t> testValues      {96,97,100,101,144,1000,26341,46341};
        std::vector<uint32_t> expectedResults {24,25, 25, 26, 34, 168, 2894, 4792};
        std::cout << "Standard sieve tests:\n";
        for (int i{0}; i<testValues.size(); ++i)
        {
            std::cout << "Prime count from 0 up to " << testValues.at(i) << " (expected  " << expectedResults.at(i)
                      << ')' << Sieve(testValues.at(i)).CountPrimes(0, testValues.at(i)) << '\n';
        }

        std::cout << "Segmented sieve tests:\n";
        for (int i{0}; i<testValues.size(); ++i)
        {
            std::cout << "Prime count from 0 up to " << testValues.at(i) << " (expected  " << expectedResults.at(i)
                      << ')' << Sieve(testValues.at(i), L1_CACHE_SIZE).CountPrimes(0, testValues.at(i))<< '\n';
        }
    }
};

int main()
{
    std::ios_base::sync_with_stdio(false);      //For faster IO
    Sieve::tests();
    for (int i{0}; i<5; ++i)
    {
        auto start = std::chrono::steady_clock::now();
        Sieve (UPPER_BOUND, L1_CACHE_SIZE);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds {end-start};
        std::cout << "The " << UPPER_BOUND << " element test " << i << " done in " << elapsed_seconds.count() << " s\n";
    }
    return 0;
}
