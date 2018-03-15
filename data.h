#ifndef DATA_H
#define DATA_H

#include <string>
using std::string;

void delim_to_array(float* data[], string filename, int length, int width, char delimiter);
void csv_to_array(float* data[], string filename, int length, int width);

#endif
