#pragma once

uint Factorial(const uint n) {
    return n > 1 ? Factorial(n-1)*n : 1;
}
