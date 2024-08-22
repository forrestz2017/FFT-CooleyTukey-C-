/*
 * Project: Fast Fourier Transform (FFT) Implementation
 * Author: Forrest Zhang
 * Date: 8/8/2024
 * Description:
 * This program implements the Cooley-Tukey Fast Fourier Transform (FFT) algorithm,
 * which efficiently computes the Discrete Fourier Transform (DFT) of a sequence.
 * The algorithm operates with a time complexity of O(n log n) and is widely used in
 * signal processing, image analysis, and other fields where frequency analysis is needed.
 */

#include <iostream>
#include <complex>
#include <vector>
#include <cmath>

 // Define a constant for PI (used in the FFT calculations)
const double PI = 3.141592653589793;

// Define complex number type using C++'s standard library
typedef std::complex<double> Complex;
typedef std::vector<Complex> CArray;

/*
 * description:
 * Reverses the order of bits in an integer.
 *
 * Parameters:
 *   - x: The integer whose bits need to be reversed.
 *   - n: The number of bits to consider in the reversal (log2 of array size).
 *
 * Returns:
 *   - An integer with the bits reversed in comparison to the input integer.
 *
 * Explanation:
 *   Bit reversal is essential in the Cooley-Tukey FFT algorithm to reorder
 *   the input array before performing the FFT. This reordering ensures that
 *   the recursive combination of smaller DFTs is done correctly.
 */
unsigned int reverse_bits(unsigned int x, unsigned int n) {
    unsigned int result = 0;
    for (unsigned int i = 0; i < n; ++i) {
        if (x & (1 << i))  // Check if the i-th bit is set
            result |= 1 << (n - 1 - i);  // Set the corresponding bit in result
    }
    return result;
}

/*
 * description:
 * Computes the Fast Fourier Transform (FFT) of the input array.
 *
 * Parameters:
 *   - x: A reference to the array of complex numbers representing the input signal.
 *
 *   The FFT function processes the input array in three main steps:
 *   1. Reorder the input array elements using bit-reversal to prepare for the FFT.
 *   2. Iteratively combine pairs of elements to form the DFT for increasingly larger subproblems.
 *   3. Perform the final combination to obtain the full DFT of the input signal.
 */
void fft(CArray& x) {
    const size_t N = x.size();  // Get the size of the input array
    const unsigned int bits = log2(N);  // Calculate the number of bits needed for bit-reversal

    // Step 1: Bit reversal reordering
    for (unsigned int i = 0; i < N; ++i) {
        unsigned int j = reverse_bits(i, bits);  // Get the bit-reversed index
        if (i < j)
            std::swap(x[i], x[j]);  // Swap elements to reorder the array
    }

    // Step 2 and 3: FFT computation
    for (size_t len = 2; len <= N; len <<= 1) {  // Loop over the size of subproblems (2, 4, 8, ...)
        double angle = -2 * PI / len;  // Calculate the angle for the complex roots of unity
        Complex wlen(cos(angle), sin(angle));  // Precompute the root of unity for the current stage
        for (size_t i = 0; i < N; i += len) {  // Loop over each subproblem
            Complex w(1);  // Initialize the complex rotation factor
            for (size_t j = 0; j < len / 2; ++j) {  // Loop over the elements within each subproblem
                Complex u = x[i + j];  // Extract the first element of the pair
                Complex v = x[i + j + len / 2] * w;  // Apply the rotation factor to the second element
                x[i + j] = u + v;  // Combine the elements (Butterfly operation)
                x[i + j + len / 2] = u - v;  // Store the difference in the second half
                w *= wlen;  // Update the rotation factor for the next pair
            }
        }
    }
}

/*
 * Function: main
 * --------------
 * Entry point of the program. Initializes a sample input array, performs the FFT,
 * and prints the results to the console.
 *
 * Explanation:
 *   This function demonstrates the usage of the FFT function by providing a simple
 *   example where the input array consists of 1s followed by 0s. The output is the
 *   frequency domain representation of the input signal.
 */
int main() {
    const int N = 8;  // Sample size (must be a power of 2 for this implementation)
    Complex test[N] = { 1, 1, 1, 1, 0, 0, 0, 0 };  // Example input array (time domain signal)

    CArray data(test, test + N);  // Initialize the data array with the sample input

    std::cout << "Input data:" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << data[i] << std::endl;  // Print the input data
    }

    fft(data);  // Perform the FFT on the input data

    std::cout << "\nFFT output:" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << data[i] << std::endl;  // Print the FFT output (frequency domain signal)
    }

    return 0;
}
