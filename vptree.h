// C++ implementation of Vantage point trees
// We store the tree as a vector

#include <vector>
#include <cstdlib>
#include <algorithm>

#include <iostream>

using namespace std;

// Default class implementing a projective distance function
//  ie. for a distance function d:T*T->double, a projector
//  for p in T implements the function d_p:T->double
//  This is useful in cases when the distance function in
//  question is defined on a transformation of the data
//  that is actually stored, since p would need to be
//  transformed for each comparison otherwise. If this is
//  the case, you should implement your own projector
template<class T, double (*d)(const T &, const T &)>
class projector {
	private:
		const T & p;
	public:
		projector(const T & elt) : p(elt) { }
		double operator()(const T & q) const { return d(p, q); }
};

template< class T, double (*d)(const T &, const T &), class projector = projector<T, d> >
class VPTree {
	private:
		vector<T> data;
		// i = vantage point, low = i + 1, high = highest element, median occurs at (i + high) / 2
		//  left contains items closer to the vantage point, right contains median and those farther
		struct node {
			unsigned int i;
			unsigned int low, high;
			double median;
			node * left;
			node * right;

			 node() { i = 0; low = high = 0; median = 0; left = right = NULL; }
			~node() { delete left; delete right; }
		} * root;

		// Defines a comparison function given a projector d_p:T->double
		//  It amounts to currying the function bool comp(d_p, a, b)
		struct comp {
			projector & dist;
			comp(projector & proj) : dist(proj) { }
			bool operator()(const T & a, const T & b) {
				return dist(a) < dist(b);
			}
		};

		// Build a tree node given its lowest and highest index
		node * build_node(unsigned int low, unsigned int high) {
			node * ret = NULL;
			unsigned int tmp;

		  if (low == high) {
				// One element in the subtree
				ret         = new node;
				ret->i      = low;
				ret->low    = low;
				ret->high   = high;
				ret->median = 0;
				ret->left   = ret->right = NULL;

			} else if (high - low == 1) {
				// Two elements in subtree

				// Pick a random element
				tmp = (rand() % (high + 1 - low)) + low;
				swap(data[low], data[tmp]);

				ret         = new node;
				ret->i      = low;
				ret->low    = high;
				ret->high   = high;
				ret->median = d(data[low], data[high]);
				ret->left   = NULL;
				ret->right  = build_node(high, high);

			} else if (low < high) {
				// Pick a random element
				tmp = (rand() % (high + 1 - low)) + low;
				swap(data[low], data[tmp]);

				// After, data[tmp] is the median distance from data[low]
				//  and everything left of data[tmp] is closer to data[low]
				tmp = (low + high) / 2;
				projector proj(data[low]);
				nth_element( data.begin() + low + 1, 
										 data.begin() + tmp, 
										 data.begin() + high, 
										 comp(proj));

				ret         = new node;
				ret->i      = low;
				ret->low    = low + 1;
				ret->high   = high;
				ret->median = d(data[low], data[tmp]);
				ret->left   = build_node(ret->low, tmp - 1);
				ret->right  = build_node(tmp, ret->high);

			}
			return ret;
		}

		// -1 indicates nothing closer than epsilon
		int NNsearch(const node * tree, const projector & dq, double * epsilon) {
			int ret = -1;
			int tmp;
			double dist;
			if (tree == NULL) {
				return ret;
			}

			dist = dq(data[tree->i]);

			if (dist < *epsilon) {
				*epsilon = dist;
				ret = tree->i;
			}

			if (dist < tree->median) {
				tmp = NNsearch(tree->left, dq, epsilon);
				if (tmp != -1) ret = tmp;
				if (dist + *epsilon >= tree->median) {
					tmp = NNsearch(tree->right, dq, epsilon);
					if (tmp != -1) ret = tmp;
				}
			} else {
				tmp = NNsearch(tree->right, dq, epsilon);
				if (tmp != -1) ret = tmp;
				if (dist - *epsilon < tree->median) {
					tmp = NNsearch(tree->left, dq, epsilon);
					if (tmp != -1) ret = tmp;
				}
			}

			return ret;
		}


	public:
		typedef typename vector<T>::iterator iterator;
		typedef typename vector<T>::const_iterator const_iterator;


		 VPTree() { root = NULL; }
		~VPTree() { delete root; data.clear(); }

		iterator begin() { return data.begin(); }
		iterator end()   { return data.end();   }
		const_iterator begin() const { return data.begin(); }
		const_iterator end()   const { return data.end();   }
		int size() const { return data.size(); }
		
		// Build a tree given a dataset
		// copies vector contents -- try to not use
		void build_tree(vector<T> & dataset) {
			delete root;
			data = dataset;
			root = build_node(0, data.size() - 1);
		}

		template<class InputIterator>
		void build_tree(InputIterator begin, InputIterator end, int size) {
			data.clear();
			delete root;

			data.reserve(size);
			for (; begin != end; begin++) data.push_back(*begin);
			root = build_node(0, data.size() - 1);
		}

		// Return the nearest neighbour and place the distance from q in epsilon
		const T * nearest_neighbour(const projector & dq, double * epsilon) {
			int index = NNsearch(root, dq, epsilon);
			if (index == -1) {
				return NULL;
			} else {
				return &data[index];
			}
		}

		const T * nearest_neighbour(const T & q, double * epsilon) {
			return nearest_neighbour(projector(q), epsilon);
		}
};
