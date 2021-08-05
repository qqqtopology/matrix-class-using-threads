#pragma once
#include <iostream>
#include <algorithm>
#include "threadMatrix.h"

point::point() { x = 0; y = 0; };

point::point(size_t x0, size_t y0)
{ x = x0; y = y0; }

cell::cell()
{ square = 0; }

cell::cell(point l0, point r0) 
{
	left = l0;
	right = r0; 
	square = (right.x - left.x) * (right.y - left.y);
}

std::vector<point> divideRectangle(const cell& rect)
{
	std::vector<point> res(4);
	res[0] = point(rect.left.x, rect.left.y);
	res[3] = point(rect.right.x, rect.right.y);
	if ((rect.right.x - rect.left.x) > (rect.right.y - rect.left.y))
	{
		res[1] = point((rect.right.x - rect.left.x) / 2 + rect.left.x, rect.right.y);
		res[2] = point((rect.right.x - rect.left.x) / 2 + rect.left.x, rect.left.y);
	}
	else
	{
		res[1] = point(rect.right.x, (rect.right.y - rect.left.y) / 2 + rect.left.y);
		res[2] = point(rect.left.x, (rect.right.y - rect.left.y) / 2 + rect.left.y);
	}
	return res;
}

bool compareRect(const cell& c1, const cell& c2)
{
	if (c1.square < c2.square)
		return true;
	else
		return false;
}

std::vector<cell> createCellForThreads(size_t n, size_t m, size_t numOfThreads)
{
	std::vector <cell> res = { cell(point(0,0), point(n,m)) };
	for (int i = 0; i < numOfThreads - 1; ++i)
	{
		std::vector<point> vec = divideRectangle(res[res.size() - 1]);
		res.pop_back();
		res.push_back(cell(vec[0], vec[1]));
		res.push_back(cell(vec[2], vec[3]));
		std::sort(res.begin(), res.end(), compareRect);
	}
	return res;
}
