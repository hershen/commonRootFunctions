

uint Factorial(const uint n) {
    return n <= 1 ? n : Factorial(n-1)*n;
}
