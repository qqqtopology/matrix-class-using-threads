#pragma once
#include <vector> 
struct point
{
	size_t x, y;

	point(); 

	point(size_t x0, size_t y0);
};

class cell
{
public:

	point left, right;

	size_t square;

	cell(); 

	cell(point l0, point r0);

};

std::vector<cell> createCellForThreads(size_t n, size_t m, size_t numOfThreads);
