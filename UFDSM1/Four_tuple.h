#pragma once
#include <unordered_map>
#include <functional>

struct  Four_tuple {
	int r1 = -2;
	int c1 = -2;
	int r2 = -2;
	int c2 = -2;
	Four_tuple(int r1_, int c1_, int r2_, int c2_) :r1(r1_), c1(c1_), r2(r2_), c2(c2_) {};
	Four_tuple() {};
	inline const Four_tuple& operator()(int r1_, int c1_, int r2_, int c2_)
	{
		r1 = r1_;
		c1 = c1_;
		r2 = r2_;
		c2 = c2_;
		return *this;
	}
	inline bool operator==(const Four_tuple & a) const
	{
		return r1 == a.r1 && c1 == a.c1 && r2 == a.r2 && c2 == a.c2;
	}

};
struct hash_name {
	size_t operator()(const Four_tuple& ft) const
	{
		return std::hash<int>()(ft.r1) ^ std::hash<int>()(ft.c1) ^ std::hash<int>()(ft.r2) ^ std::hash<int>()(ft.c2);
	};
};



