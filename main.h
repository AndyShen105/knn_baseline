/*
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 */

#ifndef MAIN_H
#define MAIN_H

/* Switches */
#define DEBUG 0
#define DUMPOUTPUT 0
#define DUMPINPUT 0
#define TIMER 1
#define BATCH 0
#define CLOCKS_PER_SECOND 1000000.0
#if DEBUG
#define DBG(x) do { x } while (0)
#else
#define DBG(x) do { } while (0)
#endif


inline void check_alloc(void *ptr) {
    if (ptr == nullptr) {
        std::cerr << "check_alloc: Failed to allocate memory" << std::endl;
        std::cin.get(); // Wait before exit
        exit(EXIT_FAILURE);
    }
}

#endif
